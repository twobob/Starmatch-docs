const canvas = document.getElementById('space-canvas');
const ctx = canvas.getContext('2d');
const slider = document.getElementById('time-slider');
const scaleSlider = document.getElementById('scale');
const dateLabel = document.getElementById('date-label');
const distanceLabel = document.getElementById('distance-label');
const legendEl = document.getElementById('legend');

const COLORS = [
  '#9ed5ff',
  '#ffc89e',
  '#a3ffba',
  '#ff9ed4',
  '#ffd17e',
  '#9fa7ff',
  '#9efff9',
  '#d1a4ff',
  '#f7ff9e',
  '#ff6666'
];

const formatDate = (jd) => {
  const unixEpochJD = 2440587.5;
  const ms = (jd - unixEpochJD) * 86400000;
  const date = new Date(ms);
  return date.toISOString().split('T')[0];
};

let data;

const DATASET_SOURCES = {
  file: [
    {
      type: 'json',
      url: new URL('./de200_demo_positions.json', import.meta.url)
    }
  ],
  http: [
    {
      type: 'integrated',
      headerUrl: new URL('../data/header.200', import.meta.url)
    },
    {
      type: 'json',
      url: new URL('../data/de200_demo_positions.json', import.meta.url)
    },
    {
      type: 'script',
      url: new URL('../data/de200_demo_positions.js', import.meta.url),
      globals: ['de200_demo_positions', 'DE200_DEMO_POSITIONS', 'demoPositions']
    }
  ]
};

async function loadDatasetFromEphemeris(headerUrl) {
  const response = await fetch(headerUrl);
  if (!response.ok) {
    throw new Error(`Request for ${headerUrl} failed with status ${response.status}`);
  }
  const text = await response.text();
  const { constants } = parseEphemerisHeader(text);
  return integrateDemoSamples(constants);
}

function getDatasetSources() {
  if (window.location.protocol === 'file:') {
    return DATASET_SOURCES.file;
  }
  return DATASET_SOURCES.http;
}

function parseEphemerisHeader(text) {
  const lines = text.split(/\r?\n/);
  let index = 0;

  function findGroup(label) {
    while (index < lines.length && lines[index].trim() !== label) {
      index += 1;
    }
    if (index >= lines.length) {
      throw new Error(`Missing ${label} in ephemeris header`);
    }
  }

  findGroup('GROUP   1030');
  index += 1; // consume group label
  if (lines[index] && lines[index].trim() !== '') {
    throw new Error('Expected blank line after GROUP 1030 header');
  }
  index += 1;
  const [startJD, endJD, stepDays] = lines[index].trim().split(/\s+/).map(Number);
  index += 1;

  findGroup('GROUP   1040');
  index += 1;
  if (lines[index] && lines[index].trim() !== '') {
    throw new Error('Expected blank line after GROUP 1040 header');
  }
  index += 1;
  const constantCount = Number.parseInt(lines[index].trim(), 10);
  index += 1;
  const constantNames = [];
  while (constantNames.length < constantCount && index < lines.length) {
    const parts = lines[index].trim().split(/\s+/).filter(Boolean);
    constantNames.push(...parts);
    index += 1;
  }

  findGroup('GROUP   1041');
  index += 1;
  if (lines[index] && lines[index].trim() !== '') {
    throw new Error('Expected blank line after GROUP 1041 header');
  }
  index += 1;
  const repeatedCount = Number.parseInt(lines[index].trim(), 10);
  if (repeatedCount !== constantCount) {
    throw new Error('Header constant count mismatch');
  }
  index += 1;
  const constantValues = [];
  while (constantValues.length < constantCount && index < lines.length) {
    const line = lines[index].trim();
    index += 1;
    if (!line) {
      continue;
    }
    const parts = line
      .split(/\s+/)
      .filter(Boolean)
      .map((token) => Number.parseFloat(token.replace(/D/i, 'E')));
    constantValues.push(...parts);
  }

  if (constantNames.length !== constantValues.length) {
    throw new Error('Header constants could not be parsed correctly');
  }

  const constants = {};
  for (let i = 0; i < constantCount; i += 1) {
    constants[constantNames[i]] = constantValues[i];
  }

  return {
    constants,
    startJD,
    endJD,
    stepDays,
  };
}

function buildBodyStates(constants) {
  const gmLookup = {
    Mercury: constants.GM1,
    Venus: constants.GM2,
    'Earth-Moon Barycenter': constants.GMB,
    Mars: constants.GM4,
    Jupiter: constants.GM5,
    Saturn: constants.GM6,
    Uranus: constants.GM7,
    Neptune: constants.GM8,
    Pluto: constants.GM9,
    Sun: constants.GMS,
  };

  const bodies = [
    ['Mercury', '1'],
    ['Venus', '2'],
    ['Earth-Moon Barycenter', 'B'],
    ['Mars', '4'],
    ['Jupiter', '5'],
    ['Saturn', '6'],
    ['Uranus', '7'],
    ['Neptune', '8'],
    ['Pluto', '9'],
    ['Sun', 'S'],
  ];

  return bodies.map(([name, token]) => ({
    name,
    position: [
      constants[`X${token}`],
      constants[`Y${token}`],
      constants[`Z${token}`],
    ],
    velocity: [
      constants[`XD${token}`],
      constants[`YD${token}`],
      constants[`ZD${token}`],
    ],
    gm: gmLookup[name],
  }));
}

function integrateDemoSamples(constants) {
  const bodies = buildBodyStates(constants);
  const n = bodies.length;

  const positions = bodies.map((body) => Float64Array.from(body.position));
  const velocities = bodies.map((body) => Float64Array.from(body.velocity));
  const masses = bodies.map((body) => body.gm);

  const stepDays = 1.0;
  const totalSteps = 720;
  const outputStride = 15;
  const samples = [];

  const auKm = constants.AU;
  const jd0 = constants.JDEPOC;

  function computeAccelerationsAt(posList) {
    const result = posList.map(() => new Float64Array(3));
    for (let i = 0; i < n; i += 1) {
      const pi = posList[i];
      const acc = result[i];
      for (let j = 0; j < n; j += 1) {
        if (i === j) {
          continue;
        }
        const pj = posList[j];
        const dx = pj[0] - pi[0];
        const dy = pj[1] - pi[1];
        const dz = pj[2] - pi[2];
        const r2 = dx * dx + dy * dy + dz * dz;
        if (r2 === 0) {
          continue;
        }
        const invR3 = masses[j] / (r2 * Math.sqrt(r2));
        acc[0] += dx * invR3;
        acc[1] += dy * invR3;
        acc[2] += dz * invR3;
      }
    }
    return result;
  }

  let currentAccelerations = computeAccelerationsAt(positions);

  for (let stepIndex = 0; stepIndex <= totalSteps; stepIndex += 1) {
    if (stepIndex % outputStride === 0) {
      const jd = jd0 + stepIndex * stepDays;
      const frame = {
        julian_date: jd,
        positions_km: {},
      };
      for (let i = 0; i < n; i += 1) {
        const coords = positions[i];
        frame.positions_km[bodies[i].name] = [
          coords[0] * auKm,
          coords[1] * auKm,
          coords[2] * auKm,
        ];
      }
      samples.push(frame);
    }

    const nextPositions = positions.map((pos, i) => {
      const vel = velocities[i];
      const acc = currentAccelerations[i];
      return Float64Array.of(
        pos[0] + vel[0] * stepDays + 0.5 * acc[0] * stepDays * stepDays,
        pos[1] + vel[1] * stepDays + 0.5 * acc[1] * stepDays * stepDays,
        pos[2] + vel[2] * stepDays + 0.5 * acc[2] * stepDays * stepDays,
      );
    });

    const nextAccelerations = computeAccelerationsAt(nextPositions);

    for (let i = 0; i < n; i += 1) {
      const vel = velocities[i];
      const acc = currentAccelerations[i];
      const nextAcc = nextAccelerations[i];
      vel[0] += 0.5 * (acc[0] + nextAcc[0]) * stepDays;
      vel[1] += 0.5 * (acc[1] + nextAcc[1]) * stepDays;
      vel[2] += 0.5 * (acc[2] + nextAcc[2]) * stepDays;
    }

    for (let i = 0; i < n; i += 1) {
      const pos = positions[i];
      const next = nextPositions[i];
      pos[0] = next[0];
      pos[1] = next[1];
      pos[2] = next[2];
    }

    currentAccelerations = nextAccelerations;
  }

  return {
    metadata: {
      description: 'Newtonian integration seeded by DE200 constants (browser)',
      start_julian_date: jd0,
      step_days: stepDays,
      output_stride_days: outputStride,
      au_km: auKm,
    },
    bodies: bodies.map((body) => body.name),
    samples,
  };
}

async function fetchJsonDataset(url) {
  const response = await fetch(url);
  if (!response.ok) {
    throw new Error(`Request for ${url} failed with status ${response.status}`);
  }

  const text = await response.text();
  try {
    return JSON.parse(text);
  } catch (error) {
    throw new Error(`Response from ${url} was not valid JSON: ${error.message}`);
  }
}

function loadScriptDataset({ url, globals }) {
  return new Promise((resolve, reject) => {
    const script = document.createElement('script');
    script.src = url.href;
    script.async = true;

    script.onload = () => {
      script.remove();
      for (const name of globals) {
        const value = window[name];
        if (value) {
          resolve(value);
          return;
        }
      }
      reject(new Error(`Loaded ${url} but none of the expected globals (${globals.join(', ')}) were defined`));
    };

    script.onerror = () => {
      script.remove();
      reject(new Error(`Failed to load ephemeris dataset script ${url}`));
    };

    document.head.appendChild(script);
  });
}

async function loadDatasetFromSource(source) {
  if (source.type === 'json') {
    return fetchJsonDataset(source.url);
  }
  if (source.type === 'script') {
    return loadScriptDataset(source);
  }
  if (source.type === 'integrated') {
    return loadDatasetFromEphemeris(source.headerUrl);
  }
  throw new Error(`Unsupported dataset source type: ${source.type}`);
}

async function loadData() {
  const sources = getDatasetSources();
  let lastError;

  for (const source of sources) {
    try {
      data = await loadDatasetFromSource(source);
      break;
    } catch (error) {
      const label = source.url?.href || source.headerUrl?.href || source.type;
      console.error(`Failed to load ephemeris data from ${label}`, error);
      lastError = error;
    }
  }

  if (!data) {
    throw lastError || new Error('No dataset source succeeded');
  }

  slider.max = data.samples.length - 1;
  buildLegend();
  updateScene();
}

function buildLegend() {
  legendEl.innerHTML = '';
  data.bodies.forEach((body, index) => {
    const row = document.createElement('div');
    const badge = document.createElement('span');
    badge.className = 'badge';
    badge.style.background = COLORS[index % COLORS.length];
    const label = document.createElement('span');
    label.textContent = body;
    row.appendChild(badge);
    row.appendChild(label);
    legendEl.appendChild(row);
  });
}

function clearCanvas() {
  ctx.clearRect(0, 0, canvas.width, canvas.height);
  ctx.fillStyle = '#05070f';
  ctx.fillRect(0, 0, canvas.width, canvas.height);
}

function drawAxes(scaleAU) {
  const midX = canvas.width / 2;
  const midY = canvas.height / 2;
  ctx.strokeStyle = 'rgba(94, 197, 255, 0.15)';
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(0, midY);
  ctx.lineTo(canvas.width, midY);
  ctx.moveTo(midX, 0);
  ctx.lineTo(midX, canvas.height);
  ctx.stroke();

  ctx.fillStyle = 'rgba(158, 206, 255, 0.5)';
  ctx.font = '12px "Fira Code", monospace';
  ctx.textAlign = 'center';
  ctx.fillText(`±${scaleAU.toFixed(1)} AU`, midX, canvas.height - 10);
}

function drawBodies(sample, scaleAU) {
  const midX = canvas.width / 2;
  const midY = canvas.height / 2;
  const auKm = data.metadata.au_km;
  const pxPerAu = (canvas.width / 2 - 40) / scaleAU;

  data.bodies.forEach((body, index) => {
    const [x, y] = sample.positions_km[body];
    const xAu = x / auKm;
    const yAu = y / auKm;
    const cx = midX + xAu * pxPerAu;
    const cy = midY - yAu * pxPerAu;

    const radius = Math.max(3, 6 - Math.log(index + 1));
    const color = COLORS[index % COLORS.length];

    const gradient = ctx.createRadialGradient(cx, cy, 1, cx, cy, radius * 2);
    gradient.addColorStop(0, color);
    gradient.addColorStop(1, 'rgba(10, 14, 30, 0.2)');

    ctx.beginPath();
    ctx.fillStyle = gradient;
    ctx.arc(cx, cy, radius, 0, Math.PI * 2);
    ctx.fill();

    if (body === 'Sun') {
      ctx.shadowColor = 'rgba(255, 210, 120, 0.7)';
      ctx.shadowBlur = 24;
      ctx.beginPath();
      ctx.fillStyle = 'rgba(255, 214, 120, 0.8)';
      ctx.arc(cx, cy, radius + 2, 0, Math.PI * 2);
      ctx.fill();
      ctx.shadowBlur = 0;
    }
  });
}

function updateScene() {
  if (!data) return;
  const index = Number(slider.value);
  const sample = data.samples[index];
  const scaleAU = Number(scaleSlider.value);

  clearCanvas();
  drawAxes(scaleAU);
  drawBodies(sample, scaleAU);

  const jd = sample.julian_date;
  dateLabel.textContent = `Julian Date ${jd.toFixed(2)} (${formatDate(jd)})`;
  const farthest = Math.max(...Object.values(sample.positions_km).map(([x, y, z]) => Math.sqrt(x * x + y * y + z * z)));
  distanceLabel.textContent = `Farthest body ≈ ${(farthest / data.metadata.au_km).toFixed(2)} AU`;
}

slider.addEventListener('input', updateScene);
scaleSlider.addEventListener('input', updateScene);

loadData().catch((error) => {
  console.error('Failed to load ephemeris demo data', error);
  dateLabel.textContent = 'Failed to load data';
  distanceLabel.textContent = error.message;
});
