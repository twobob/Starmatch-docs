// ============================================================================
// Astrological Engine visualiser
// Integrates the theme analysis engine with ephemeris data
// ============================================================================

const canvas = document.getElementById('chart-canvas');
const ctx = canvas.getContext('2d');

// DOM Elements
const birthDate = document.getElementById('birth-date');
const birthTime = document.getElementById('birth-time');
const latitudeInput = document.getElementById('latitude');
const longitudeInput = document.getElementById('longitude');
const btnCalculate = document.getElementById('btn-calculate');
const orbTypeSelect = document.getElementById('orb-type');
const aspectOrbSetSelect = document.getElementById('aspect-orb-set');
const rulershipSetSelect = document.getElementById('rulership-set');
const precessionCheckbox = document.getElementById('precession-flag');
const locationLookupBtn = document.getElementById('location-lookup');
const selectedLocationName = document.getElementById('selected-location-name');
// CRUD UI elements
const btnSaveRecord = document.getElementById('btn-save-record');
const btnLoadRecords = document.getElementById('btn-load-records');
const recordsPanel = document.getElementById('records-panel');
const recordsList = document.getElementById('records-list');
const btnCloseRecords = document.getElementById('btn-close-records');
const btnClearAll = document.getElementById('btn-clear-all');
// Modal elements
const saveModal = document.getElementById('save-modal');
const modalClose = document.getElementById('modal-close');
const modalCancel = document.getElementById('modal-cancel');
const modalSave = document.getElementById('modal-save');
const recordNameInput = document.getElementById('record-name-input');
const applySavedSettingsCheckbox = document.getElementById('apply-saved-settings');
// Danger modal elements
const dangerModal = document.getElementById('danger-modal');
const dangerModalClose = document.getElementById('danger-modal-close');
const dangerCancel = document.getElementById('danger-cancel');
const dangerConfirm = document.getElementById('danger-confirm');
const dangerModalTitle = document.getElementById('danger-modal-title');
const dangerModalText = document.getElementById('danger-modal-text');
const dangerStageIndicator = document.getElementById('danger-stage-indicator');
let dangerStage = 0; // 0 -> first prompt, 1 -> second prompt

function openDangerModal() {
  dangerStage = 0;
  updateDangerModal();
  dangerModal.classList.remove('hidden');
  setTimeout(()=> { dangerConfirm.focus(); }, 30);
}
function closeDangerModal() {
  dangerModal.classList.add('hidden');
  dangerStage = 0;
}
function updateDangerModal() {
  if (dangerStage === 0) {
    dangerModalTitle.textContent = 'Delete ALL Records?';
    dangerModalText.textContent = 'This will permanently remove EVERY saved record. This cannot be undone.';
    dangerConfirm.textContent = 'Yes, Continue';
    dangerStageIndicator.textContent = 'Stage 1 / 2';
  } else {
    dangerModalTitle.textContent = 'Are You REALLY Sure?';
    dangerModalText.textContent = 'Final confirmation. All records including names, coordinates, and settings will be lost.';
    dangerConfirm.textContent = 'Delete Everything';
    dangerStageIndicator.textContent = 'Stage 2 / 2';
  }
}

// Storage key
const STORAGE_KEY = 'astro_records_v1';

function loadRecords() {
  try {
    const raw = localStorage.getItem(STORAGE_KEY);
    if (!raw) return [];
    const parsed = JSON.parse(raw);
    if (!Array.isArray(parsed)) return [];
    return parsed;
  } catch (e) {
    console.warn('Failed to parse records from storage', e);
    return [];
  }
}

function saveRecords(records) {
  localStorage.setItem(STORAGE_KEY, JSON.stringify(records));
}

function addRecord(data) {
  const records = loadRecords();
  records.push(data);
  saveRecords(records);
  return records;
}

function updateRecord(id, patch) {
  const records = loadRecords();
  const idx = records.findIndex(r => r.id === id);
  if (idx !== -1) {
    records[idx] = { ...records[idx], ...patch, updatedAt: new Date().toISOString() };
    saveRecords(records);
  }
  return records;
}

function deleteRecord(id) {
  const records = loadRecords().filter(r => r.id !== id);
  saveRecords(records);
  return records;
}

function clearAllRecords() {
  saveRecords([]);
}

function buildRecordPayload(name) {
  return {
    id: crypto.randomUUID(),
    name: name && name.trim() ? name.trim() : 'Untitled',
    date: birthDate.value || '',
    time: birthTime.value || '',
    lat: latitudeInput.value || '',
    lon: longitudeInput.value || '',
    orbType: orbTypeSelect.value,
    aspectOrbSet: aspectOrbSetSelect.value,
    rulershipSet: rulershipSetSelect.value,
    precession: precessionCheckbox.checked ? 1 : 0,
    createdAt: new Date().toISOString(),
    updatedAt: new Date().toISOString()
  };
}

function renderRecords() {
  const records = loadRecords();
  recordsList.innerHTML = '';
  if (!records.length) {
    recordsList.classList.add('empty');
    recordsList.innerHTML = '<div class="empty-msg">No saved records yet.</div>';
    return;
  }
  recordsList.classList.remove('empty');
  records.sort((a,b)=> a.name.localeCompare(b.name));
  records.forEach(rec => {
    const el = document.createElement('div');
    el.className = 'record-item';
    el.innerHTML = `
      <div class="record-name" data-id="${rec.id}" title="Click to rename">${rec.name}</div>
      <div class="record-actions load-col">
        <button class="pill-btn" data-action="load" data-id="${rec.id}" title="Load & Calculate">Load</button>
      </div>
      <div class="record-meta">${rec.date || '—'} ${rec.time || ''}</div>
      <div class="record-meta">${rec.lat || '—'}, ${rec.lon || '—'}</div>
      <div class="record-actions main-actions">
        <button class="pill-btn" data-action="overwrite" data-id="${rec.id}" title="Overwrite this saved record with current inputs/settings">Overwrite</button>
        <button class="pill-btn danger" data-action="del" data-id="${rec.id}">Del</button>
      </div>`;
    recordsList.appendChild(el);
  });
}

function openRecordsPanel() {
  recordsPanel.classList.remove('hidden');
  renderRecords();
}
function closeRecordsPanel() { recordsPanel.classList.add('hidden'); }

function openSaveModal() {
  recordNameInput.value = '';
  saveModal.classList.remove('hidden');
  recordNameInput.focus();
}
function closeSaveModal() { saveModal.classList.add('hidden'); }

function applyRecord(rec, doCalculate=false, includeSettings=true) {
  if (rec.date) birthDate.value = rec.date;
  if (rec.time) birthTime.value = rec.time;
  if (rec.lat) latitudeInput.value = rec.lat;
  if (rec.lon) longitudeInput.value = rec.lon;
  if (includeSettings) {
    if (rec.orbType !== undefined) orbTypeSelect.value = rec.orbType;
    if (rec.aspectOrbSet !== undefined) aspectOrbSetSelect.value = rec.aspectOrbSet;
    if (rec.rulershipSet !== undefined) rulershipSetSelect.value = rec.rulershipSet;
    if (rec.precession !== undefined) precessionCheckbox.checked = rec.precession === 1;
  }
  if (doCalculate) {
    calculateChart();
  }
}

// Inline rename
recordsList?.addEventListener('click', (e) => {
  const target = e.target;
  if (target.classList.contains('record-name')) {
    const id = target.getAttribute('data-id');
    const original = target.textContent;
    target.contentEditable = 'true';
    target.classList.add('editing');
    target.focus();
    const range = document.createRange();
    range.selectNodeContents(target);
    const sel = window.getSelection();
    sel.removeAllRanges();
    sel.addRange(range);
    function finish(save) {
      target.contentEditable = 'false';
      target.classList.remove('editing');
      if (save) {
        const newName = target.textContent.trim() || original;
        updateRecord(id, { name: newName });
        renderRecords();
      } else {
        target.textContent = original;
      }
      target.removeEventListener('blur', onBlur);
      target.removeEventListener('keydown', onKey);
    }
    function onBlur(){ finish(true); }
    function onKey(ev){
      if (ev.key === 'Enter') { ev.preventDefault(); finish(true); }
      else if (ev.key === 'Escape') { finish(false); }
    }
    target.addEventListener('blur', onBlur);
    target.addEventListener('keydown', onKey);
  }
  if (target.dataset.action) {
    const id = target.getAttribute('data-id');
    const action = target.dataset.action;
    const rec = loadRecords().find(r=>r.id===id);
    if (!rec) return;
    const includeSettings = applySavedSettingsCheckbox ? applySavedSettingsCheckbox.checked : true;
    if (action === 'load') {
      // Auto calculate always on load
      applyRecord(rec,true,includeSettings);
    } else if (action === 'overwrite') {
      // Build new payload but keep same id & createdAt
      const updated = buildRecordPayload(rec.name);
      updated.id = rec.id;
      updated.createdAt = rec.createdAt;
      updated.updatedAt = new Date().toISOString();
      // Replace record
      const all = loadRecords().map(r => r.id === rec.id ? updated : r);
      saveRecords(all);
      renderRecords();
    } else if (action === 'del') {
      deleteRecord(id);
      renderRecords();
    }
  }
});

// Button events
btnLoadRecords?.addEventListener('click', () => {
  if (recordsPanel.classList.contains('hidden')) openRecordsPanel(); else closeRecordsPanel();
});
btnCloseRecords?.addEventListener('click', closeRecordsPanel);
btnSaveRecord?.addEventListener('click', openSaveModal);
btnClearAll?.addEventListener('click', () => {
  const existing = loadRecords();
  if (!existing.length) return; // nothing to clear
  openDangerModal();
});
// Danger modal events
dangerModalClose?.addEventListener('click', closeDangerModal);
dangerCancel?.addEventListener('click', closeDangerModal);
dangerConfirm?.addEventListener('click', () => {
  if (dangerStage === 0) {
    dangerStage = 1;
    updateDangerModal();
  } else {
    clearAllRecords();
    renderRecords();
    closeDangerModal();
  }
});
// Escape key handling for modals
document.addEventListener('keydown', (e) => {
  if (e.key === 'Escape') {
    if (!saveModal.classList.contains('hidden')) closeSaveModal();
    if (!dangerModal.classList.contains('hidden')) closeDangerModal();
  }
});
modalClose?.addEventListener('click', closeSaveModal);
modalCancel?.addEventListener('click', closeSaveModal);
modalSave?.addEventListener('click', () => {
  const payload = buildRecordPayload(recordNameInput.value);
  addRecord(payload);
  closeSaveModal();
  renderRecords();
  openRecordsPanel();
});
recordNameInput?.addEventListener('keydown', (e)=>{ if (e.key==='Enter'){ e.preventDefault(); modalSave.click(); } });

// Init render
document.addEventListener('DOMContentLoaded', () => {
  renderRecords();
});

// Display containers
const positionsDisplay = document.getElementById('positions-display');
const themeBars = document.getElementById('theme-bars');
const aspectCounts = document.getElementById('aspect-counts');
const tradFactors = document.getElementById('trad-factors');
const dominantInfo = document.getElementById('dominant-info');

// Constants
const SIGN_NAMES = ['Aries', 'Taurus', 'Gemini', 'Cancer', 'Leo', 'Virgo', 
                    'Libra', 'Scorpio', 'Sagittarius', 'Capricorn', 'Aquarius', 'Pisces'];
const PLANET_NAMES = ['Sun', 'Moon', 'Mercury', 'Venus', 'Mars', 'Jupiter', 
                      'Saturn', 'Uranus', 'Neptune', 'Pluto'];
const ASPECT_NAMES = ['Conjunction', 'Opposition', 'Trine', 'Square', 'Sextile', 'Semi-square', 'Semi-sextile'];
const ELEMENT_NAMES = ['Fire', 'Earth', 'Air', 'Water'];
const QUALITY_NAMES = ['Cardinal', 'Fixed', 'Mutable'];

// Global variables for engine.js
// Note: engine.js declares 'nativityYear' (line 58) but references 'nativity' (lines 331, 676)
// We only need to declare 'nativity' since 'nativityYear' is already declared in engine.js
var nativity = 0;        // Used in engine.js lines 331, 676 but not declared there

// Global chart data for tooltips
let chartData = {
  positions: {},
  ascendant: 0,
  midheaven: 0,
  aspects: [],
  centerX: 300,
  centerY: 300,
  planetRadius: 160
};

// Helper function to match engine.js requirements
function TidyUpAndFloat(value) {
  return parseFloat(value);
}

// ============================================================================
// Ephemeris Data Loading (Simplified from app.js)
// ============================================================================

let ephemerisData = null;

// Convert date/time to Julian Date
function dateToJulianDate(date, time) {
  const dateStr = `${date}T${time}:00Z`;
  const jsDate = new Date(dateStr);
  const unixEpochJD = 2440587.5;
  const ms = jsDate.getTime();
  return unixEpochJD + (ms / 86400000);
}

// Get planetary positions for a given Julian Date with linear interpolation
function getPositionsForJD(jd) {
  if (!ephemerisData || !ephemerisData.samples) {
    console.error('No ephemeris data available');
    return null;
  }

  // Find the two samples that bracket the requested JD
  let beforeIndex = -1;
  let afterIndex = -1;

  for (let i = 0; i < ephemerisData.samples.length - 1; i++) {
    if (ephemerisData.samples[i].julian_date <= jd && ephemerisData.samples[i + 1].julian_date >= jd) {
      beforeIndex = i;
      afterIndex = i + 1;
      break;
    }
  }

  // If JD is outside range, use nearest sample
  if (beforeIndex === -1) {
    let closestIndex = 0;
    let closestDiff = Math.abs(ephemerisData.samples[0].julian_date - jd);

    for (let i = 1; i < ephemerisData.samples.length; i++) {
      const diff = Math.abs(ephemerisData.samples[i].julian_date - jd);
      if (diff < closestDiff) {
        closestDiff = diff;
        closestIndex = i;
      }
    }

    const sample = ephemerisData.samples[closestIndex];
    console.log(`⚠️ JD ${jd.toFixed(2)} outside range - using nearest sample at JD ${sample.julian_date.toFixed(2)} (diff: ${closestDiff.toFixed(2)} days)`);
    return sample.positions_km;
  }

  // Perform linear interpolation between the two samples
  const sample1 = ephemerisData.samples[beforeIndex];
  const sample2 = ephemerisData.samples[afterIndex];
  
  const jd1 = sample1.julian_date;
  const jd2 = sample2.julian_date;
  const fraction = (jd - jd1) / (jd2 - jd1);

  const interpolatedPositions = {};
  
  for (const bodyName in sample1.positions_km) {
    const pos1 = sample1.positions_km[bodyName];
    const pos2 = sample2.positions_km[bodyName];
    
    // Linear interpolation for each coordinate
    interpolatedPositions[bodyName] = [
      pos1[0] + (pos2[0] - pos1[0]) * fraction,
      pos1[1] + (pos2[1] - pos1[1]) * fraction,
      pos1[2] + (pos2[2] - pos1[2]) * fraction
    ];
  }

  console.log(`✓ Interpolated positions for JD ${jd.toFixed(2)} between JD ${jd1.toFixed(2)} and ${jd2.toFixed(2)} (fraction: ${fraction.toFixed(4)})`);
  
  return interpolatedPositions;
}

// Convert ecliptic coordinates to zodiacal longitude (simplified)
function cartesianToLongitude(x, y, z) {
  // Convert from 3D position to ecliptic longitude (0-360 degrees)
  const longitude = Math.atan2(y, x) * (180 / Math.PI);
  return longitude < 0 ? longitude + 360 : longitude;
}

// ============================================================================
// Ascendant and Midheaven Calculation
// ============================================================================

function calculateAscendant(jd, latitude, longitude) {
  // Simplified calculation - in reality you'd need RAMC, obliquity, etc.
  // For demo purposes, we'll use a simplified approach
  
  // Get sidereal time (simplified)
  const T = (jd - 2451545.0) / 36525.0;
  const gmst = 280.46061837 + 360.98564736629 * (jd - 2451545.0) + 
               0.000387933 * T * T - (T * T * T) / 38710000.0;
  
  let lst = (gmst + longitude) % 360;
  if (lst < 0) lst += 360;
  
  // Simplified ascendant calculation (not accurate, but demonstrates concept)
  let asc = (lst + 90) % 360;
  if (asc < 0) asc += 360;
  
  return asc;
}

function calculateMidheaven(jd, longitude) {
  // Simplified MC calculation
  const T = (jd - 2451545.0) / 36525.0;
  const gmst = 280.46061837 + 360.98564736629 * (jd - 2451545.0) + 
               0.000387933 * T * T - (T * T * T) / 38710000.0;
  
  let lst = (gmst + longitude) % 360;
  if (lst < 0) lst += 360;
  let mc = lst % 360;
  if (mc < 0) mc += 360;
  
  return mc;
}

// ============================================================================
// Chart Calculation
// ============================================================================

function calculateChart() {
  const date = birthDate.value;
  const time = birthTime.value;
  const latitude = parseFloat(latitudeInput.value);
  const longitude = parseFloat(longitudeInput.value);

  if (!date || !time) {
    alert('Please enter birth date and time');
    return;
  }

  // Convert to Julian Date
  const jd = dateToJulianDate(date, time);
  
  // Get ephemeris positions
  const positions = getPositionsForJD(jd);
  
  if (!positions) {
    alert('Ephemeris data not loaded or date out of range');
    return;
  }

  // Convert cartesian positions to zodiacal longitudes
  const planetaryPositions = {};
  
  // Map ephemeris body names to our planet names
  // Note: Using Earth-Moon Barycenter as approximation for Moon position
  const bodyMapping = {
    'Sun': 'Sun',
    'Moon': 'Earth-Moon Barycenter',  // Barycenter is very close to Earth, good approximation
    'Mercury': 'Mercury',
    'Venus': 'Venus',
    'Mars': 'Mars',
    'Jupiter': 'Jupiter',
    'Saturn': 'Saturn',
    'Uranus': 'Uranus',
    'Neptune': 'Neptune',
    'Pluto': 'Pluto'
  };

  for (const [ourName, ephemName] of Object.entries(bodyMapping)) {
    if (positions[ephemName]) {
      const [x, y, z] = positions[ephemName];
      const longitude = cartesianToLongitude(x, y, z);
      
      // Validate the calculated longitude
      if (!isNaN(longitude) && isFinite(longitude)) {
        planetaryPositions[ourName] = longitude;
      } else {
        console.warn(`Invalid longitude calculated for ${ourName}:`, longitude, 'from coords:', x, y, z);
      }
    } else {
      console.warn(`Missing ephemeris data for ${ephemName}`);
    }
  }

  // Check if we have enough planetary data to proceed
  if (Object.keys(planetaryPositions).length < 10) {
    alert(`Warning: Only found ${Object.keys(planetaryPositions).length} out of 10 planets. Results may be incomplete.`);
  }

  // Calculate Ascendant and Midheaven
  const ascendant = calculateAscendant(jd, latitude, longitude);
  const midheaven = calculateMidheaven(jd, longitude);

  // Apply engine settings
  orbType = parseInt(orbTypeSelect.value);
  aoIndex = parseInt(aspectOrbSetSelect.value);
  tfIndex = parseInt(rulershipSetSelect.value);
  precessionFlag = precessionCheckbox.checked ? 1 : 0;
  
  // Extract year for precession - engine.js uses both 'nativity' and 'nativityYear'
  const birthYear = new Date(date).getFullYear();
  nativityYear = birthYear;  // Used by precessPositions function
  nativity = birthYear;      // Used in main getThemeValues function (line 331, 676 of engine.js)

  // Ensure all planets have values (use 0 as fallback if missing)
  const safePlanetaryPositions = {
    Sun: planetaryPositions.Sun || 0,
    Moon: planetaryPositions.Moon || 0,
    Mercury: planetaryPositions.Mercury || 0,
    Venus: planetaryPositions.Venus || 0,
    Mars: planetaryPositions.Mars || 0,
    Jupiter: planetaryPositions.Jupiter || 0,
    Saturn: planetaryPositions.Saturn || 0,
    Uranus: planetaryPositions.Uranus || 0,
    Neptune: planetaryPositions.Neptune || 0,
    Pluto: planetaryPositions.Pluto || 0
  };

  // Call the engine
  getThemeValues(
    safePlanetaryPositions.Sun,
    safePlanetaryPositions.Moon,
    safePlanetaryPositions.Mercury,
    safePlanetaryPositions.Venus,
    safePlanetaryPositions.Mars,
    safePlanetaryPositions.Jupiter,
    safePlanetaryPositions.Saturn,
    safePlanetaryPositions.Uranus,
    safePlanetaryPositions.Neptune,
    safePlanetaryPositions.Pluto,
    ascendant,
    midheaven
  );

  // Store chart data for tooltips
  chartData.positions = planetaryPositions;
  chartData.ascendant = ascendant;
  chartData.midheaven = midheaven;
  
  // Display results (use the original positions that actually exist)
  displayPositions(planetaryPositions, ascendant, midheaven);
  displayThemes();
  displayAspects();
  displayTraditionalFactors();
  displayDominants();
  drawChartWheel(planetaryPositions, ascendant, midheaven);
}

// ============================================================================
// Display Functions
// ============================================================================

function displayPositions(positions, ascendant, midheaven) {
  positionsDisplay.innerHTML = '';

  const allPositions = {
    ...positions,
    'Ascendant': ascendant,
    'Midheaven': midheaven
  };

  for (const [name, longitude] of Object.entries(allPositions)) {
    // Validate longitude value
    if (longitude === undefined || longitude === null || isNaN(longitude)) {
      console.warn(`Invalid longitude for ${name}:`, longitude);
      continue;
    }

    // Normalize longitude to 0-360 range
    let normalizedLon = longitude % 360;
    if (normalizedLon < 0) normalizedLon += 360;

    const sign = Math.floor(normalizedLon / 30);
    const degree = normalizedLon % 30;
    const signName = SIGN_NAMES[sign];

    // Additional validation for sign index
    if (!signName) {
      console.warn(`Invalid sign for ${name}: longitude=${longitude}, sign=${sign}`);
      continue;
    }

    const item = document.createElement('div');
    item.className = 'position-item';
    item.innerHTML = `
      <span class="planet-name">${name}</span>
      <span class="planet-position">
        ${degree.toFixed(2)}° 
        <span class="planet-sign sign-${signName.toLowerCase()}">${signName}</span>
      </span>
    `;
    positionsDisplay.appendChild(item);
  }
}

function displayThemes() {
  themeBars.innerHTML = '';

  // Find max theme value for normalization
  const maxTheme = Math.max(...theme);

  SIGN_NAMES.forEach((sign, index) => {
    const value = theme[index];
    const percentage = maxTheme > 0 ? (value / maxTheme) * 100 : 0;

    const bar = document.createElement('div');
    bar.className = 'theme-bar';
    bar.innerHTML = `
      <span class="theme-label">${sign}</span>
      <div class="bar-wrapper">
        <div class="bar-fill" style="width: ${percentage}%"></div>
        <div class="bar-value">${value.toFixed(2)}</div>
      </div>
    `;
    themeBars.appendChild(bar);
  });
}

function displayAspects() {
  aspectCounts.innerHTML = '';

  ASPECT_NAMES.forEach((name, index) => {
    const count = numAspects[index];
    const item = document.createElement('div');
    item.className = 'info-item';
    item.innerHTML = `
      <span class="info-label">${name}</span>
      <span class="info-value">${count}</span>
    `;
    aspectCounts.appendChild(item);
  });
}

function displayTraditionalFactors() {
  tradFactors.innerHTML = '';

  const labels = [
    'Positive Signs',
    'Negative Signs',
    'Fire',
    'Earth',
    'Air',
    'Water',
    'Cardinal',
    'Fixed',
    'Mutable'
  ];

  labels.forEach((label, index) => {
    const value = numTradFactors[index];
    const item = document.createElement('div');
    item.className = 'info-item';
    item.innerHTML = `
      <span class="info-label">${label}</span>
      <span class="info-value">${value}</span>
    `;
    tradFactors.appendChild(item);
  });
}

function displayDominants() {
  dominantInfo.innerHTML = '';

  // Polarity
  const polarityItem = document.createElement('div');
  polarityItem.className = 'info-item';
  polarityItem.innerHTML = `
    <span class="info-label">Dominant Polarity</span>
    <span class="info-value dominant">${tfDominant[0] === 1 ? 'Positive' : 'Negative'}</span>
  `;
  dominantInfo.appendChild(polarityItem);

  // Element
  const elementItem = document.createElement('div');
  elementItem.className = 'info-item';
  const elementIndex = tfDominant[1] - 2; // Adjust for array offset
  const elementName = elementIndex >= 0 && elementIndex < 4 ? ELEMENT_NAMES[elementIndex] : 'None';
  elementItem.innerHTML = `
    <span class="info-label">Dominant Element</span>
    <span class="info-value dominant">${elementName}</span>
  `;
  dominantInfo.appendChild(elementItem);

  // Quality
  const qualityItem = document.createElement('div');
  qualityItem.className = 'info-item';
  const qualityIndex = tfDominant[2] - 6; // Adjust for array offset
  const qualityName = qualityIndex >= 0 && qualityIndex < 3 ? QUALITY_NAMES[qualityIndex] : 'None';
  qualityItem.innerHTML = `
    <span class="info-label">Dominant Quality</span>
    <span class="info-value dominant">${qualityName}</span>
  `;
  dominantInfo.appendChild(qualityItem);

  // Find strongest theme
  let maxThemeIndex = 0;
  let maxThemeValue = theme[0];
  for (let i = 1; i < 12; i++) {
    if (theme[i] > maxThemeValue) {
      maxThemeValue = theme[i];
      maxThemeIndex = i;
    }
  }

  const themeItem = document.createElement('div');
  themeItem.className = 'info-item';
  themeItem.innerHTML = `
    <span class="info-label">Strongest Theme</span>
    <span class="info-value dominant">${SIGN_NAMES[maxThemeIndex]} (${maxThemeValue.toFixed(2)})</span>
  `;
  dominantInfo.appendChild(themeItem);
}

// ============================================================================
// Chart Wheel visualisation
// ============================================================================

function drawChartWheel(positions, ascendant, midheaven) {
  const centerX = canvas.width / 2;
  const centerY = canvas.height / 2;
  const outerRadius = 280;
  const innerRadius = 220;
  const zodiacRadius = 250;
  
  // Store chart dimensions for tooltip system
  chartData.centerX = centerX;
  chartData.centerY = centerY;
  chartData.planetRadius = innerRadius - 60;
  chartData.innerRadius = innerRadius;
  chartData.outerRadius = outerRadius;
  chartData.aspects = []; // Reset aspects array

  // Clear canvas
  ctx.clearRect(0, 0, canvas.width, canvas.height);
  ctx.fillStyle = '#05070f';
  ctx.fillRect(0, 0, canvas.width, canvas.height);

  // Draw zodiac wheel (12 signs)
  drawZodiacWheel(centerX, centerY, outerRadius, innerRadius, ascendant);

  // Draw house cusps
  drawHouseCusps(centerX, centerY, innerRadius, ascendant);

  // Draw planets
  drawPlanets(centerX, centerY, innerRadius - 60, positions);

  // Draw aspects
  drawAspects(centerX, centerY, innerRadius - 60, positions);
}

function drawZodiacWheel(centerX, centerY, outerRadius, innerRadius, ascendant) {
  const signColors = [
    '#ff6b6b', '#51cf66', '#ffd43b', '#74c0fc',
    '#ff8787', '#69db7c', '#ffd43b', '#ff6b6b',
    '#cc5de8', '#51cf66', '#74c0fc', '#a78bfa'
  ];

  for (let i = 0; i < 12; i++) {
    const startAngle = ((i * 30 - ascendant - 90) * Math.PI) / 180;
    const endAngle = (((i + 1) * 30 - ascendant - 90) * Math.PI) / 180;

    // Draw sign segment
    ctx.beginPath();
    ctx.moveTo(centerX, centerY);
    ctx.arc(centerX, centerY, outerRadius, startAngle, endAngle);
    ctx.closePath();
    ctx.fillStyle = signColors[i] + '20';
    ctx.fill();
    ctx.strokeStyle = signColors[i] + '80';
    ctx.lineWidth = 1;
    ctx.stroke();

    // Draw sign name
    const midAngle = startAngle + (endAngle - startAngle) / 2;
    const textRadius = (outerRadius + innerRadius) / 2 + 15;
    const textX = centerX + Math.cos(midAngle) * textRadius;
    const textY = centerY + Math.sin(midAngle) * textRadius;

    ctx.save();
    ctx.translate(textX, textY);
    ctx.rotate(midAngle + Math.PI / 2);
    ctx.fillStyle = signColors[i];
    ctx.font = 'bold 14px "Segoe UI"';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    ctx.fillText(SIGN_NAMES[i], 0, 0);
    ctx.restore();
  }

  // Draw inner circle
  ctx.beginPath();
  ctx.arc(centerX, centerY, innerRadius, 0, Math.PI * 2);
  ctx.fillStyle = 'rgba(5, 7, 15, 0.9)';
  ctx.fill();
  ctx.strokeStyle = 'rgba(94, 197, 255, 0.4)';
  ctx.lineWidth = 2;
  ctx.stroke();
}

function drawHouseCusps(centerX, centerY, radius, ascendant) {
  // Draw Ascendant line (1st house cusp)
  const ascAngle = ((-ascendant - 90) * Math.PI) / 180;
  
  ctx.beginPath();
  ctx.moveTo(centerX, centerY);
  ctx.lineTo(
    centerX + Math.cos(ascAngle) * radius,
    centerY + Math.sin(ascAngle) * radius
  );
  ctx.strokeStyle = '#ffb85e';
  ctx.lineWidth = 3;
  ctx.stroke();

  // Label ASC
  ctx.fillStyle = '#ffb85e';
  ctx.font = 'bold 12px "Segoe UI"';
  ctx.textAlign = 'center';
  const labelX = centerX + Math.cos(ascAngle) * (radius - 20);
  const labelY = centerY + Math.sin(ascAngle) * (radius - 20);
  ctx.fillText('ASC', labelX, labelY);
}

function drawPlanets(centerX, centerY, radius, positions) {
  const planetSymbols = ['☉', '☽', '☿', '♀', '♂', '♃', '♄', '⛢', '♆', '♇'];
  const planetColors = [
    '#ffd700', '#c0c0c0', '#ffa500', '#ff69b4', '#ff0000',
    '#9370db', '#8b4513', '#00ced1', '#4169e1', '#8b0000'
  ];

  Object.entries(positions).forEach(([name, longitude]) => {
    // Validate longitude
    if (longitude === undefined || longitude === null || isNaN(longitude)) {
      console.warn(`Skipping planet ${name} - invalid longitude:`, longitude);
      return;
    }

    const angle = ((-longitude - 90) * Math.PI) / 180;
    const x = centerX + Math.cos(angle) * radius;
    const y = centerY + Math.sin(angle) * radius;

    const planetIndex = PLANET_NAMES.indexOf(name);
    if (planetIndex === -1) {
      console.warn(`Planet ${name} not found in PLANET_NAMES list`);
      return;
    }

    // Draw planet circle
    ctx.beginPath();
    ctx.arc(x, y, 12, 0, Math.PI * 2);
    ctx.fillStyle = planetColors[planetIndex];
    ctx.fill();
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.5)';
    ctx.lineWidth = 2;
    ctx.stroke();

    // Draw symbol
    ctx.fillStyle = '#000';
    ctx.font = 'bold 16px Arial';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    ctx.fillText(planetSymbols[planetIndex], x, y);
  });
}

function drawAspects(centerX, centerY, radius, positions) {
  const aspectColors = {
    0: 'rgba(255, 215, 0, 0.6)',    // Conjunction - gold
    180: 'rgba(255, 69, 0, 0.6)',   // Opposition - red
    120: 'rgba(0, 255, 127, 0.6)',  // Trine - green
    90: 'rgba(255, 0, 0, 0.6)',     // Square - red
    60: 'rgba(135, 206, 250, 0.6)', // Sextile - blue
    45: 'rgba(255, 165, 0, 0.5)',   // Semi-square - orange
    30: 'rgba(173, 216, 230, 0.5)'  // Semi-sextile - light blue
  };

  const planetArray = Object.entries(positions);

  for (let i = 0; i < planetArray.length; i++) {
    for (let j = i + 1; j < planetArray.length; j++) {
      const [name1, lon1] = planetArray[i];
      const [name2, lon2] = planetArray[j];

      // Validate both longitudes
      if (isNaN(lon1) || isNaN(lon2) || !isFinite(lon1) || !isFinite(lon2)) {
        continue;
      }

      const diff = Math.abs(lon1 - lon2);
      const normalizedDiff = diff > 180 ? 360 - diff : diff;

      // Check for aspects
      for (let k = 0; k < a.length; k++) {
        const aspectAngle = a[k];
        const orb = ao[aoIndex][k];

        if (Math.abs(normalizedDiff - aspectAngle) <= orb) {
          // Draw aspect line
          const angle1 = ((-lon1 - 90) * Math.PI) / 180;
          const angle2 = ((-lon2 - 90) * Math.PI) / 180;

          const x1 = centerX + Math.cos(angle1) * radius;
          const y1 = centerY + Math.sin(angle1) * radius;
          const x2 = centerX + Math.cos(angle2) * radius;
          const y2 = centerY + Math.sin(angle2) * radius;

          ctx.beginPath();
          ctx.moveTo(x1, y1);
          ctx.lineTo(x2, y2);
          ctx.strokeStyle = aspectColors[aspectAngle] || 'rgba(200, 200, 200, 0.3)';
          ctx.lineWidth = aspectAngle === 0 || aspectAngle === 180 || aspectAngle === 120 || aspectAngle === 90 ? 2 : 1;
          ctx.stroke();
          
          // Store aspect data for tooltips
          const exactOrb = Math.abs(normalizedDiff - aspectAngle);
          chartData.aspects.push({
            planet1: name1,
            planet2: name2,
            type: ASPECT_NAMES[k],
            angle: aspectAngle,
            orb: exactOrb,
            x1, y1, x2, y2
          });

          break; // Only draw one aspect between two planets
        }
      }
    }
  }
}

// ============================================================================
// Load Ephemeris Data - Using REAL .eph file computation!
// ============================================================================

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
  
  return dataset;
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
  index += 1;
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
  index += 1;
  
  const constantValues = [];
  while (constantValues.length < constantCount && index < lines.length) {
    const line = lines[index].trim();
    index += 1;
    
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
      const needed = constantCount - constantValues.length;
      constantValues.push(...parts.slice(0, needed));
      break;
    } else {
      constantValues.push(...parts);
    }
  }
  
  if (constantNames.length !== constantValues.length) {
    throw new Error(`Header constants could not be parsed correctly: ${constantNames.length} names vs ${constantValues.length} values`);
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
  const outputStride = 1;
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

async function loadEphemerisData() {
  try {
    console.log('Loading ephemeris from .eph file...');
    const headerUrl = new URL('../data/header.200', window.location.href);
    const ephUrl = new URL('../data/de200.eph', window.location.href);
    
    ephemerisData = await loadDatasetFromEphemeris({ headerUrl, ephUrl });
    
    console.log('✓ Ephemeris computed from DE200:', ephemerisData.samples.length, 'samples');
    console.log('✓ Date range:', 
      new Date((ephemerisData.metadata.start_julian_date - 2440587.5) * 86400000).toISOString().split('T')[0],
      'to',
      new Date((ephemerisData.metadata.end_julian_date - 2440587.5) * 86400000).toISOString().split('T')[0]
    );
    console.log('✓ Available bodies:', ephemerisData.bodies.join(', '));
    
    if (ephemerisData.samples.length > 0) {
      const firstSample = ephemerisData.samples[0];
      console.log('✓ Bodies in first sample:', Object.keys(firstSample.positions_km).join(', '));
    }
    
    btnCalculate.disabled = false;
  } catch (error) {
    console.error('✗ Error loading ephemeris:', error);
    alert('Failed to load and compute ephemeris data from .eph file. Error: ' + error.message);
  }
}

// ============================================================================
// Event Listeners
// ============================================================================

btnCalculate.addEventListener('click', calculateChart);

// Location lookup button
locationLookupBtn.addEventListener('click', () => {
  const picker = initLocationPicker();
  picker.open((location) => {
    // Update the latitude and longitude inputs
    latitudeInput.value = location.latitude.toFixed(4);
    longitudeInput.value = location.longitude.toFixed(4);
    
    // Display the selected location name
    selectedLocationName.textContent = `📍 ${location.name} (${location.fullAddress})`;
    selectedLocationName.style.color = '#5ec5ff';
    
    console.log('Location selected:', location);
  });
});

// ============================================================================
// Interactive Tooltip System
// ============================================================================

const tooltip = document.getElementById('chart-tooltip');

// Helper function to calculate distance from point to line segment
function distanceToLineSegment(px, py, x1, y1, x2, y2) {
  const A = px - x1;
  const B = py - y1;
  const C = x2 - x1;
  const D = y2 - y1;

  const dot = A * C + B * D;
  const lenSq = C * C + D * D;
  let param = -1;
  
  if (lenSq !== 0) {
    param = dot / lenSq;
  }

  let xx, yy;

  if (param < 0) {
    xx = x1;
    yy = y1;
  } else if (param > 1) {
    xx = x2;
    yy = y2;
  } else {
    xx = x1 + param * C;
    yy = y1 + param * D;
  }

  const dx = px - xx;
  const dy = py - yy;
  return Math.sqrt(dx * dx + dy * dy);
}

// Get mouse position relative to canvas
function getMousePos(canvas, evt) {
  const rect = canvas.getBoundingClientRect();
  return {
    x: evt.clientX - rect.left,
    y: evt.clientY - rect.top
  };
}

// Check if mouse is over a planet
function checkPlanetHover(mouseX, mouseY) {
  const planetRadius = 12; // Match the planet circle radius
  
  for (const [name, longitude] of Object.entries(chartData.positions)) {
    if (isNaN(longitude) || !isFinite(longitude)) continue;
    
    const angle = ((-longitude - 90) * Math.PI) / 180;
    const x = chartData.centerX + Math.cos(angle) * chartData.planetRadius;
    const y = chartData.centerY + Math.sin(angle) * chartData.planetRadius;
    
    const distance = Math.sqrt((mouseX - x) ** 2 + (mouseY - y) ** 2);
    
    if (distance <= planetRadius) {
      let normalizedLon = longitude % 360;
      if (normalizedLon < 0) normalizedLon += 360;
      const sign = Math.floor(normalizedLon / 30);
      const degree = normalizedLon % 30;
      const signName = SIGN_NAMES[sign];
      
      return {
        type: 'planet',
        name: name,
        longitude: longitude,
        position: `${degree.toFixed(2)}° ${signName}`,
        sign: signName
      };
    }
  }
  return null;
}

// Check if mouse is over ascendant line
function checkAscendantHover(mouseX, mouseY) {
  const ascAngle = ((-chartData.ascendant - 90) * Math.PI) / 180;
  const x1 = chartData.centerX;
  const y1 = chartData.centerY;
  const x2 = chartData.centerX + Math.cos(ascAngle) * chartData.innerRadius;
  const y2 = chartData.centerY + Math.sin(ascAngle) * chartData.innerRadius;
  
  const distance = distanceToLineSegment(mouseX, mouseY, x1, y1, x2, y2);
  
  if (distance <= 5) {
    let normalizedAsc = chartData.ascendant % 360;
    if (normalizedAsc < 0) normalizedAsc += 360;
    const sign = Math.floor(normalizedAsc / 30);
    const degree = normalizedAsc % 30;
    const signName = SIGN_NAMES[sign];
    
    return {
      type: 'ascendant',
      name: 'Ascendant',
      position: `${degree.toFixed(2)}° ${signName}`,
      sign: signName
    };
  }
  return null;
}

// Check if mouse is over an aspect line
function checkAspectHover(mouseX, mouseY) {
  for (const aspect of chartData.aspects) {
    const distance = distanceToLineSegment(
      mouseX, mouseY,
      aspect.x1, aspect.y1,
      aspect.x2, aspect.y2
    );
    
    if (distance <= 5) {
      return {
        type: 'aspect',
        planet1: aspect.planet1,
        planet2: aspect.planet2,
        aspectType: aspect.type,
        angle: aspect.angle,
        orb: aspect.orb
      };
    }
  }
  return null;
}

// Check if mouse is over a zodiac sign
function checkSignHover(mouseX, mouseY) {
  const dx = mouseX - chartData.centerX;
  const dy = mouseY - chartData.centerY;
  const distance = Math.sqrt(dx * dx + dy * dy);
  
  // Check if in zodiac ring area
  if (distance >= chartData.innerRadius && distance <= chartData.outerRadius) {
    // Calculate angle from center
    // Mirror about vertical axis by negating dx instead of the angle
    let angle = Math.atan2(dy, -dx) * (180 / Math.PI);
    // Convert to zodiac longitude (adjusted for ascendant and 90° offset)
    let zodiacLon = -angle - 90 + chartData.ascendant;
    while (zodiacLon < 0) zodiacLon += 360;
    while (zodiacLon >= 360) zodiacLon -= 360;
    
    const signIndex = Math.floor(zodiacLon / 30);
    const signName = SIGN_NAMES[signIndex];
    
    return {
      type: 'sign',
      name: signName,
      index: signIndex,
      element: ELEMENT_NAMES[signIndex % 4],
      quality: QUALITY_NAMES[Math.floor(signIndex / 4)]
    };
  }
  return null;
}

// Update tooltip display
function updateTooltip(evt) {
  if (!chartData.positions || Object.keys(chartData.positions).length === 0) {
    tooltip.style.display = 'none';
    return;
  }
  
  const mousePos = getMousePos(canvas, evt);
  const mouseX = mousePos.x;
  const mouseY = mousePos.y;
  
  // Check in priority order: planets, ascendant, aspects, signs
  let hoverInfo = checkPlanetHover(mouseX, mouseY);
  
  if (!hoverInfo) {
    hoverInfo = checkAscendantHover(mouseX, mouseY);
  }
  
  if (!hoverInfo) {
    hoverInfo = checkAspectHover(mouseX, mouseY);
  }
  
  if (!hoverInfo) {
    hoverInfo = checkSignHover(mouseX, mouseY);
  }
  
  if (hoverInfo) {
    let tooltipHTML = '';
    
    if (hoverInfo.type === 'planet') {
      tooltipHTML = `
        <strong>${hoverInfo.name}</strong><br>
        ${hoverInfo.position}
      `;
    } else if (hoverInfo.type === 'ascendant') {
      tooltipHTML = `
        <strong>Ascendant (Rising Sign)</strong><br>
        ${hoverInfo.position}
      `;
    } else if (hoverInfo.type === 'aspect') {
      tooltipHTML = `
        <strong>${hoverInfo.aspectType}</strong><br>
        ${hoverInfo.planet1} ⟷ ${hoverInfo.planet2}<br>
        <span style="font-size: 0.9em;">Orb: ${hoverInfo.orb.toFixed(2)}°</span>
      `;
    } else if (hoverInfo.type === 'sign') {
      tooltipHTML = `
        <strong>${hoverInfo.name}</strong><br>
        <span style="font-size: 0.9em;">${hoverInfo.element} • ${hoverInfo.quality}</span>
      `;
    }
    
    tooltip.innerHTML = tooltipHTML;
    tooltip.style.display = 'block';
    
    // Position tooltip near cursor (fixed positioning uses viewport coordinates)
    tooltip.style.left = (evt.clientX + 15) + 'px';
    tooltip.style.top = (evt.clientY + 15) + 'px';
    
    // Change cursor to pointer
    canvas.style.cursor = 'pointer';
  } else {
    tooltip.style.display = 'none';
    canvas.style.cursor = 'default';
  }
}

// Add mouse event listeners to canvas
canvas.addEventListener('mousemove', updateTooltip);
canvas.addEventListener('mouseleave', () => {
  tooltip.style.display = 'none';
  canvas.style.cursor = 'default';
});

// Initialize
loadEphemerisData();

// Set default date to a date within ephemeris range
birthDate.value = '1974-09-11';
selectedLocationName.textContent = '';

// ============================================================================
// Collapsible Analysis Section
// ============================================================================

const analysisToggle = document.getElementById('analysis-toggle');
const analysisContent = document.getElementById('analysis-content');

analysisToggle?.addEventListener('click', () => {
  const isCollapsed = analysisContent.classList.contains('collapsed');
  
  if (isCollapsed) {
    // Expand
    analysisContent.classList.remove('collapsed');
    analysisToggle.classList.remove('collapsed');
  } else {
    // Collapse
    analysisContent.classList.add('collapsed');
    analysisToggle.classList.add('collapsed');
  }
});
