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

async function loadData() {
  const response = await fetch('../data/de200_demo_positions.json');
  data = await response.json();
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
});
