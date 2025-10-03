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

function getDatasetSources() {
  if (window.location.protocol === 'file:') {
    return DATASET_SOURCES.file;
  }
  return DATASET_SOURCES.http;
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
      console.error(`Failed to load ephemeris data from ${source.url}`, error);
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
