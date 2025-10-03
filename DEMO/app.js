const canvas = document.getElementById('space-canvas');
const ctx = canvas.getContext('2d');
const slider = document.getElementById('time-slider');
const scaleSlider = document.getElementById('scale');
const speedSlider = document.getElementById('speed-slider');
const dateLabel = document.getElementById('date-label');
const distanceLabel = document.getElementById('distance-label');
const speedLabel = document.getElementById('speed-label');
const legendEl = document.getElementById('legend');

const btnStart = document.getElementById('btn-start');
const btnReverse = document.getElementById('btn-reverse');
const btnStop = document.getElementById('btn-stop');
const btnPlay = document.getElementById('btn-play');
const btnEnd = document.getElementById('btn-end');

let hoveredBody = null;
let bodyPositions = []; // Store body positions for hover detection
let playbackState = 'stopped'; // 'stopped', 'playing', 'reversing'
let playbackInterval = null;
let playbackSpeed = 1; // days per second
let virtualTime = 0; // Continuous time for smooth interpolation

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
    }
    // Removed JSON and script fallbacks - force ephemeris computation only
    // {
    //   type: 'json',
    //   url: resolveRelativeUrl('../data/de200_demo_positions.json')
    // },
    // {
    //   type: 'script',
    //   url: resolveRelativeUrl('../data/de200_demo_positions.js'),
    //   globals: ['de200_demo_positions', 'DE200_DEMO_POSITIONS', 'demoPositions']
    // }
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
  const protocol = window.location.protocol;
  
  if (protocol === 'file:') {
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
  const repeatedCountLine = lines[index].trim();
  const repeatedCount = Number.parseInt(repeatedCountLine, 10);
  if (repeatedCount !== constantCount) {
    throw new Error('Header constant count mismatch');
  }
  index += 1; // Skip the count line itself
  
  const constantValues = [];
  while (constantValues.length < constantCount && index < lines.length) {
    const line = lines[index].trim();
    const currentLineNumber = index + 1; // 1-indexed for readability
    index += 1;
    
    // Stop if we hit a GROUP marker or empty line followed by GROUP
    if (line.startsWith('GROUP') || !line) {
      if (!line && index < lines.length && lines[index].trim().startsWith('GROUP')) {
        break;
      }
      if (!line) {
        continue;
      }
      break;
    }
    
    const parts = line
      .split(/\s+/)
      .filter(Boolean)
      .map((token) => Number.parseFloat(token.replace(/D/i, 'E')));
    
    if (constantValues.length + parts.length > constantCount) {
      // Only add enough values to reach the expected count
      const needed = constantCount - constantValues.length;
      constantValues.push(...parts.slice(0, needed));
      break;
    } else {
      constantValues.push(...parts);
    }
  }
  
  if (constantNames.length !== constantValues.length) {
    throw new Error(`Header constants could not be parsed correctly: ${constantNames.length} names vs ${constantValues.length} values (expected ${constantCount})`);
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
  const outputStride = 1; // Output every day for smooth animation
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
      console.log(`✓ Loaded ${data.samples.length} samples from ${data.bodies.length} bodies (${describeSource(source)})`);
      slider.max = data.samples.length - 1;
      buildLegend();
      updateScene();
      return;
    } catch (error) {
      lastError = error;
      console.warn(`Failed to load from ${describeSource(source)}:`, error.message);
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

  // Center the view on the Sun
  const centerBody = 'Sun';
  const centerPos = sample.positions_km[centerBody];
  const centerX = centerPos ? centerPos[0] : 0;
  const centerY = centerPos ? centerPos[1] : 0;

  // Clear the body positions array for hover detection
  bodyPositions = [];

  data.bodies.forEach((body, index) => {
    const [x, y] = sample.positions_km[body];
    // Subtract center position to make the view relative to the center body
    const xAu = (x - centerX) / auKm;
    const yAu = (y - centerY) / auKm;
    const cx = midX + xAu * pxPerAu;
    const cy = midY - yAu * pxPerAu;

    const radius = Math.max(3, 6 - Math.log(index + 1));
    const color = COLORS[index % COLORS.length];

    // Store position for hover detection
    bodyPositions.push({ name: body, cx, cy, radius: radius + 3 });

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

  // Draw hover label if a body is hovered
  if (hoveredBody) {
    const pos = bodyPositions.find(p => p.name === hoveredBody);
    if (pos) {
      ctx.fillStyle = 'rgba(0, 0, 0, 0.8)';
      ctx.strokeStyle = 'rgba(255, 255, 255, 0.9)';
      ctx.lineWidth = 1;
      ctx.font = '14px "Fira Code", monospace';
      ctx.textAlign = 'center';
      
      const text = hoveredBody;
      const metrics = ctx.measureText(text);
      const padding = 6;
      const labelWidth = metrics.width + padding * 2;
      const labelHeight = 20;
      const labelX = pos.cx - labelWidth / 2;
      const labelY = pos.cy - pos.radius - labelHeight - 5;
      
      // Draw label background
      ctx.fillRect(labelX, labelY, labelWidth, labelHeight);
      ctx.strokeRect(labelX, labelY, labelWidth, labelHeight);
      
      // Draw label text
      ctx.fillStyle = 'rgba(255, 255, 255, 0.95)';
      ctx.fillText(text, pos.cx, labelY + 14);
    }
  }
}

function updateScene() {
  if (!data) return;
  
  // Use the nearest actual sample - never interpolate orbital positions!
  const timeValue = playbackState !== 'stopped' ? virtualTime : Number(slider.value);
  const index = Math.round(timeValue);
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

// Speed slider with logarithmic scale
// 0-100 maps to realtime (very slow) to 100 years/sec
function updateSpeedFromSlider() {
  const sliderValue = Number(speedSlider.value);
  
  // Logarithmic mapping: 
  // 0 = 0.5 days/sec (one sample ~every second if samples are 15 days apart)
  // 100 = 36500 days/sec (100 years/sec)
  const minSpeed = 0.5; // 0.5 days per second
  const maxSpeed = 36500; // 100 years per second
  const logMin = Math.log(minSpeed);
  const logMax = Math.log(maxSpeed);
  const scale = (logMax - logMin) / 100;
  
  playbackSpeed = Math.exp(logMin + scale * sliderValue);
  
  // Update label
  if (playbackSpeed < 1) {
    speedLabel.textContent = `${playbackSpeed.toFixed(2)} day/sec`;
  } else if (playbackSpeed < 365) {
    speedLabel.textContent = `${playbackSpeed.toFixed(1)} days/sec`;
  } else {
    speedLabel.textContent = `${(playbackSpeed / 365).toFixed(1)} years/sec`;
  }
}

speedSlider.addEventListener('input', updateSpeedFromSlider);
updateSpeedFromSlider(); // Initialize

// Playback controls
function stopPlayback() {
  if (playbackInterval) {
    clearInterval(playbackInterval);
    playbackInterval = null;
  }
  playbackState = 'stopped';
  updateButtonStates();
}

function startPlayback(direction) {
  stopPlayback();
  playbackState = direction === 'forward' ? 'playing' : 'reversing';
  updateButtonStates();
  
  // Initialize virtual time to current slider position
  virtualTime = Number(slider.value);
  
  const fps = 30; // Target 30 fps
  const msPerFrame = 1000 / fps;
  
  playbackInterval = setInterval(() => {
    if (!data) return;
    
    // Calculate how many samples to advance based on playback speed
    const daysPerSample = data.metadata.output_stride_days || 15;
    const daysPerFrame = playbackSpeed / fps;
    const samplesPerFrame = daysPerFrame / daysPerSample;
    
    // Update virtual time continuously
    if (direction === 'forward') {
      virtualTime += samplesPerFrame;
    } else {
      virtualTime -= samplesPerFrame;
    }
    
    // Clamp and stop at boundaries
    if (virtualTime >= slider.max) {
      virtualTime = slider.max;
      slider.value = slider.max;
      stopPlayback();
    } else if (virtualTime <= 0) {
      virtualTime = 0;
      slider.value = 0;
      stopPlayback();
    } else {
      // Update slider to reflect current position (fractional values allowed internally)
      slider.value = Math.round(virtualTime);
    }
    
    updateScene();
  }, msPerFrame);
}

function updateButtonStates() {
  btnPlay.classList.toggle('active', playbackState === 'playing');
  btnReverse.classList.toggle('active', playbackState === 'reversing');
  btnStop.classList.toggle('active', playbackState === 'stopped');
}

btnStart.addEventListener('click', () => {
  stopPlayback();
  slider.value = 0;
  updateScene();
});

btnEnd.addEventListener('click', () => {
  stopPlayback();
  slider.value = slider.max;
  updateScene();
});

btnPlay.addEventListener('click', () => {
  if (playbackState === 'playing') {
    stopPlayback();
  } else {
    startPlayback('forward');
  }
});

btnReverse.addEventListener('click', () => {
  if (playbackState === 'reversing') {
    stopPlayback();
  } else {
    startPlayback('reverse');
  }
});

btnStop.addEventListener('click', () => {
  stopPlayback();
});

// Add mouse move listener for hover detection
canvas.addEventListener('mousemove', (event) => {
  const rect = canvas.getBoundingClientRect();
  const mouseX = event.clientX - rect.left;
  const mouseY = event.clientY - rect.top;
  
  let foundBody = null;
  for (const pos of bodyPositions) {
    const dx = mouseX - pos.cx;
    const dy = mouseY - pos.cy;
    const distance = Math.sqrt(dx * dx + dy * dy);
    
    if (distance <= pos.radius) {
      foundBody = pos.name;
      break;
    }
  }
  
  if (foundBody !== hoveredBody) {
    hoveredBody = foundBody;
    canvas.style.cursor = hoveredBody ? 'pointer' : 'default';
    updateScene();
  }
});

// Clear hover when mouse leaves canvas
canvas.addEventListener('mouseleave', () => {
  if (hoveredBody) {
    hoveredBody = null;
    canvas.style.cursor = 'default';
    updateScene();
  }
});

loadData().catch((error) => {
  console.error('Failed to load ephemeris demo data', error);
  dateLabel.textContent = 'Failed to load data';
  distanceLabel.textContent = error.message;
});

// Initialize button states
updateButtonStates();
