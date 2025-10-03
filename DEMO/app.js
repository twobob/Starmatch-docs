const canvas = document.getElementById('space-canvas');
const ctx = canvas.getContext('2d');
const slider = document.getElementById('time-slider');
const scaleSlider = document.getElementById('scale');
const dateLabel = document.getElementById('date-label');
const distanceLabel = document.getElementById('distance-label');
const legendEl = document.getElementById('legend');

const currentScript = document.currentScript;
const scriptBaseUrl = currentScript ? currentScript.src : window.location.href;
const resolveRelativeUrl = (path) => new URL(path, scriptBaseUrl);

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

let data = window.DE200_DEMO_POSITIONS || null;

const DATASET_SOURCES = {
  file: [
    {
      type: 'json',
      url: resolveRelativeUrl('./de200_demo_positions.json')
    }
  ],
  http: [
    {
      type: 'ephemeris',
      headerUrl: resolveRelativeUrl('../data/header.200'),
      ephUrl: resolveRelativeUrl('../data/de200.eph')
    },
    {
      type: 'json',
      url: resolveRelativeUrl('../data/de200_demo_positions.json')
    },
    {
      type: 'script',
      url: resolveRelativeUrl('../data/de200_demo_positions.js'),
      globals: ['de200_demo_positions', 'DE200_DEMO_POSITIONS', 'demoPositions']
    }
  ]
};

async function computeEphemerisUsage(buffer) {
  const bytes = new Uint8Array(buffer);
  let checksum32 = 0;
  for (let i = 0; i < bytes.length; i += 1) {
    checksum32 = (checksum32 + bytes[i]) >>> 0;
  }

  let sha256 = null;
  if (globalThis.crypto && globalThis.crypto.subtle) {
    try {
      const digest = await globalThis.crypto.subtle.digest('SHA-256', buffer);
      sha256 = Array.from(new Uint8Array(digest))
        .map((byte) => byte.toString(16).padStart(2, '0'))
        .join('');
    } catch (error) {
      console.warn('Failed to compute SHA-256 for ephemeris payload', error);
    }
  }

  return {
    byteLength: bytes.length,
    checksum32: `0x${checksum32.toString(16).padStart(8, '0')}`,
    sha256
  };
}

async function loadDatasetFromEphemeris({ headerUrl, ephUrl }) {
  const [headerResponse, ephResponse] = await Promise.all([
    fetch(headerUrl),
    fetch(ephUrl)
  ]);

  if (!headerResponse.ok) {
    throw new Error(`Request for ${headerUrl} failed with status ${headerResponse.status}`);
  }
  if (!ephResponse.ok) {
    throw new Error(`Request for ${ephUrl} failed with status ${ephResponse.status}`);
  }

  const [headerText, ephBuffer] = await Promise.all([
    headerResponse.text(),
    ephResponse.arrayBuffer()
  ]);

  const { constants, startJD, endJD } = parseEphemerisHeader(headerText);
  const dataset = integrateDemoSamples(constants, { startJD, endJD });
  const usage = await computeEphemerisUsage(ephBuffer);
  dataset.metadata.ephemeris_asset = {
    header_url: headerUrl.href,
    eph_url: ephUrl.href,
    ...usage
  };
  return dataset;
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

function integrateDemoSamples(constants, options = {}) {
  const bodies = buildBodyStates(constants);
  const n = bodies.length;

  const positions = bodies.map((body) => Float64Array.from(body.position));
  const velocities = bodies.map((body) => Float64Array.from(body.velocity));
  const masses = bodies.map((body) => body.gm);

  const stepDays = 1.0;
  const outputStride = 15;
  const samples = [];

  const auKm = constants.AU;
  const jd0 = constants.JDEPOC;

  const defaultTotalSteps = 720;
  const hasHeaderRange = Number.isFinite(options.startJD) && Number.isFinite(options.endJD);
  let rangeStartJD = hasHeaderRange ? options.startJD : jd0;
  let rangeEndJD = hasHeaderRange ? options.endJD : jd0 + defaultTotalSteps * stepDays;

  if (!Number.isFinite(rangeStartJD)) {
    rangeStartJD = jd0;
  }
  if (!Number.isFinite(rangeEndJD)) {
    rangeEndJD = rangeStartJD;
  }
  if (rangeEndJD < rangeStartJD) {
    const tmp = rangeStartJD;
    rangeStartJD = rangeEndJD;
    rangeEndJD = tmp;
  }

  function computeAccelerationsAt(posList, target) {
    let result;
    if (target) {
      result = target;
      for (let i = 0; i < result.length; i += 1) {
        const acc = result[i];
        acc[0] = 0;
        acc[1] = 0;
        acc[2] = 0;
      }
    } else {
      result = posList.map(() => new Float64Array(3));
    }
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

  const currentAccelerations = computeAccelerationsAt(positions);
  const nextAccelerations = positions.map(() => new Float64Array(3));
  const nextPositions = positions.map(() => new Float64Array(3));

  let currentJD = jd0;

  function advanceState(dt) {
    const dt2 = dt * dt;
    for (let i = 0; i < n; i += 1) {
      const pos = positions[i];
      const vel = velocities[i];
      const acc = currentAccelerations[i];
      const next = nextPositions[i];
      next[0] = pos[0] + vel[0] * dt + 0.5 * acc[0] * dt2;
      next[1] = pos[1] + vel[1] * dt + 0.5 * acc[1] * dt2;
      next[2] = pos[2] + vel[2] * dt + 0.5 * acc[2] * dt2;
    }

    computeAccelerationsAt(nextPositions, nextAccelerations);

    for (let i = 0; i < n; i += 1) {
      const vel = velocities[i];
      const acc = currentAccelerations[i];
      const nextAcc = nextAccelerations[i];
      vel[0] += 0.5 * (acc[0] + nextAcc[0]) * dt;
      vel[1] += 0.5 * (acc[1] + nextAcc[1]) * dt;
      vel[2] += 0.5 * (acc[2] + nextAcc[2]) * dt;

      const pos = positions[i];
      const next = nextPositions[i];
      pos[0] = next[0];
      pos[1] = next[1];
      pos[2] = next[2];

      const accTarget = currentAccelerations[i];
      const srcAcc = nextAccelerations[i];
      accTarget[0] = srcAcc[0];
      accTarget[1] = srcAcc[1];
      accTarget[2] = srcAcc[2];
    }

    currentJD += dt;
  }

  function moveToJulianDate(targetJD) {
    if (!Number.isFinite(targetJD) || targetJD === currentJD) {
      return;
    }
    let remaining = targetJD - currentJD;
    const direction = Math.sign(remaining) || 1;
    const step = stepDays * direction;
    const steps = Math.floor(Math.abs(remaining) / stepDays);
    for (let i = 0; i < steps; i += 1) {
      advanceState(step);
    }
    remaining = targetJD - currentJD;
    if (Math.abs(remaining) > 1e-9) {
      advanceState(remaining);
    }
  }

  function addSample(jd) {
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

  moveToJulianDate(rangeStartJD);
  addSample(currentJD);

  const totalSpan = rangeEndJD - rangeStartJD;
  if (totalSpan <= 0) {
    return {
      metadata: {
        description: 'Newtonian integration seeded by DE200 constants (browser)',
        start_julian_date: samples[0].julian_date,
        end_julian_date: samples[0].julian_date,
        step_days: stepDays,
        output_stride_days: outputStride,
        au_km: auKm,
      },
      bodies: bodies.map((body) => body.name),
      samples,
    };
  }

  const totalSteps = Math.floor(totalSpan / stepDays);
  const remainderDays = totalSpan - totalSteps * stepDays;

  let stepsSinceSample = 0;
  for (let stepIndex = 0; stepIndex < totalSteps; stepIndex += 1) {
    advanceState(stepDays);
    stepsSinceSample += 1;
    if (stepsSinceSample >= outputStride) {
      addSample(currentJD);
      stepsSinceSample = 0;
    }
  }

  if (remainderDays > 1e-9) {
    advanceState(remainderDays);
    stepsSinceSample += remainderDays / stepDays;
  }

  const lastSample = samples[samples.length - 1];
  if (!lastSample || Math.abs(lastSample.julian_date - currentJD) > 1e-9) {
    addSample(currentJD);
  }

  return {
    metadata: {
      description: 'Newtonian integration seeded by DE200 constants (browser)',
      start_julian_date: samples[0].julian_date,
      end_julian_date: samples[samples.length - 1].julian_date,
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
  if (source.type === 'ephemeris') {
    return loadDatasetFromEphemeris(source);
  }
  throw new Error(`Unsupported dataset source type: ${source.type}`);
}

async function loadData() {
  if (data) {
    slider.max = data.samples.length - 1;
    buildLegend();
    updateScene();
    return;
  }

  const sources = getDatasetSources();
  let lastError = null;
  for (const source of sources) {
    try {
      data = await loadDatasetFromSource(source);
      slider.max = data.samples.length - 1;
      buildLegend();
      updateScene();
      return;
    } catch (error) {
      lastError = error;
      console.warn(`Failed to load dataset from ${describeSource(source)}`, error);
    }
  }

  if (lastError) {
    throw lastError;
  }
  throw new Error('No dataset sources available');
}

function describeSource(source) {
  if (source.type === 'json') {
    return `JSON ${source.url}`;
  }
  if (source.type === 'script') {
    return `script ${source.url}`;
  }
  if (source.type === 'ephemeris') {
    return `ephemeris ${source.ephUrl}`;
  }
  return source.type;
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
