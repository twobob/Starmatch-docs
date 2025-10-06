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
const dateRangeHint = document.getElementById('date-range-hint');
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
    
    // Create name element
    const nameEl = document.createElement('div');
    nameEl.className = 'record-name';
    nameEl.setAttribute('data-id', rec.id);
    nameEl.setAttribute('title', 'Click to rename');
    nameEl.textContent = rec.name;
    
    // Create actions row (all buttons together)
    const actionsRow = document.createElement('div');
    actionsRow.className = 'record-actions-row';
    actionsRow.innerHTML = `
      <button class="pill-btn" data-action="load" data-id="${rec.id}" title="Load & Calculate">Load</button>
      <button class="pill-btn" data-action="overwrite" data-id="${rec.id}" title="Overwrite this saved record with current inputs/settings">Overwrite</button>
      <button class="pill-btn danger" data-action="del" data-id="${rec.id}">Del</button>`;
    
    // Create metadata row
    const metaRow = document.createElement('div');
    metaRow.className = 'record-meta-row';
    metaRow.innerHTML = `
      <span class="record-meta">${rec.date || '—'} ${rec.time || ''}</span>
      <span class="record-meta">${rec.lat || '—'}, ${rec.lon || '—'}</span>`;
    
    // Append all elements in correct order
    el.appendChild(nameEl);
    el.appendChild(actionsRow);
    el.appendChild(metaRow);
    
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
  const bodyMapping = {
    'Sun': 'Sun',
    'Moon': 'Moon',  // DE200 has separate Moon ephemeris
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
      let longitude;
      
      // Special case for Sun: In heliocentric coords, Sun is at origin (0,0,0)
      // So Sun's apparent position from Earth is Earth's position + 180°
      if (ourName === 'Sun') {
        // Use Earth-Moon Barycenter as proxy for Earth's position
        const earthPos = positions['Earth-Moon Barycenter'];
        if (earthPos) {
          // Sun as seen from Earth is opposite direction from Earth's heliocentric position
          longitude = cartesianToLongitude(-earthPos[0], -earthPos[1], -earthPos[2]);
        } else {
          console.warn('Earth-Moon Barycenter position not found for Sun calculation');
          continue;
        }
      } else {
        longitude = cartesianToLongitude(x, y, z);
      }
      
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

  // Load app-chebyshev.js module functions
  const { parseEphemerisBuffer, generateSamplesFromEphemeris } = window;
  if (!parseEphemerisBuffer || !generateSamplesFromEphemeris) {
    throw new Error('app-chebyshev.js not loaded');
  }

  const { constants } = parseEphemerisHeader(headerText);
  const ephemeris = parseEphemerisBuffer(ephBuffer);
  const dataset = generateSamplesFromEphemeris(ephemeris, constants);
  
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
  const protocol = window.location.protocol;
  
  // Determine which data source to use based on protocol
  if (protocol === 'file:') {
    // When running from file://, load the pre-computed data from a script
    // (fetch doesn't work with file:// due to CORS restrictions)
    try {
      console.log('Running from file:// protocol - loading pre-computed data via script...');
      
      const scriptUrl = new URL('../data/de200_demo_positions.js', window.location.href);
      
      ephemerisData = await new Promise((resolve, reject) => {
        const script = document.createElement('script');
        script.src = scriptUrl.href;
        script.async = true;

        script.onload = () => {
          script.remove();
          // Check for the global variable that the script should define
          const globals = ['de200_demo_positions', 'DE200_DEMO_POSITIONS', 'demoPositions'];
          for (const name of globals) {
            const value = window[name];
            if (value) {
              console.log(`✓ Loaded ephemeris data from global variable: ${name}`);
              resolve(value);
              return;
            }
          }
          reject(new Error(`Script loaded but none of the expected globals (${globals.join(', ')}) were defined`));
        };

        script.onerror = () => {
          script.remove();
          reject(new Error(`Failed to load ephemeris dataset script from ${scriptUrl}`));
        };

        document.head.appendChild(script);
      });
      
      console.log('✓ Loaded ephemeris data:', ephemerisData.samples.length, 'samples');
      
      // Get date range from metadata or samples
      let startJD, endJD;
      if (ephemerisData.metadata && ephemerisData.metadata.start_julian_date && ephemerisData.metadata.end_julian_date) {
        startJD = ephemerisData.metadata.start_julian_date;
        endJD = ephemerisData.metadata.end_julian_date;
      } else if (ephemerisData.samples && ephemerisData.samples.length > 0) {
        // Fallback: get range from actual samples
        startJD = ephemerisData.samples[0].julian_date;
        endJD = ephemerisData.samples[ephemerisData.samples.length - 1].julian_date;
        console.log('ℹ️ Using date range from samples (metadata not available)');
      } else {
        throw new Error('Cannot determine date range from ephemeris data');
      }
      
      const startDate = new Date((startJD - 2440587.5) * 86400000);
      const endDate = new Date((endJD - 2440587.5) * 86400000);
      const startDateStr = startDate.toISOString().split('T')[0];
      const endDateStr = endDate.toISOString().split('T')[0];
      
      console.log('✓ Date range:', startDateStr, 'to', endDateStr);
      console.log('✓ Available bodies:', ephemerisData.bodies.join(', '));
      
      if (ephemerisData.samples.length > 0) {
        const firstSample = ephemerisData.samples[0];
        console.log('✓ Bodies in first sample:', Object.keys(firstSample.positions_km).join(', '));
      }
      
      // Set default date to middle of available range
      const midJD = (startJD + endJD) / 2;
      const midDate = new Date((midJD - 2440587.5) * 86400000);
      birthDate.value = midDate.toISOString().split('T')[0];
      
      // Set date input constraints
      birthDate.min = startDateStr;
      birthDate.max = endDateStr;
      birthDate.title = `Available range: ${startDateStr} to ${endDateStr}`;
      
      // Show date range hint
      if (dateRangeHint) {
        dateRangeHint.textContent = `📅 Available: ${startDateStr} to ${endDateStr}`;
        dateRangeHint.style.display = 'block';
      }
      
      btnCalculate.disabled = false;
    } catch (error) {
      console.error('✗ Error loading ephemeris data:', error);
      alert('Failed to load ephemeris data from script file. Error: ' + error.message);
    }
  } else {
    // When running from http://, compute from .eph file
    try {
      console.log('Running from http:// protocol - computing ephemeris from .eph file...');
      const urlParams = new URLSearchParams(location.search);
      const useReduced = urlParams.get('reduced') !== '0'; // Use reduced by default
      
      const headerUrl = new URL(useReduced ? '../data/reduced_header.200' : '../data/header.200', window.location.href);
      const ephUrl = new URL(useReduced ? '../data/reduced_de200.eph' : '../data/de200.eph', window.location.href);
      
      ephemerisData = await loadDatasetFromEphemeris({ headerUrl, ephUrl });
      
      console.log('✓ Ephemeris computed from DE200:', ephemerisData.samples.length, 'samples');
      
      const startDate = new Date((ephemerisData.metadata.start_julian_date - 2440587.5) * 86400000);
      const endDate = new Date((ephemerisData.metadata.end_julian_date - 2440587.5) * 86400000);
      const startDateStr = startDate.toISOString().split('T')[0];
      const endDateStr = endDate.toISOString().split('T')[0];
      
      console.log('✓ Date range:', startDateStr, 'to', endDateStr);
      console.log('✓ Available bodies:', ephemerisData.bodies.join(', '));
      
      if (ephemerisData.samples.length > 0) {
        const firstSample = ephemerisData.samples[0];
        console.log('✓ Bodies in first sample:', Object.keys(firstSample.positions_km).join(', '));
        
        // Data quality check: Display positions for a few bodies
        console.log('\n📊 DATA QUALITY CHECK - First Sample:');
        const AU_TO_KM = 149597870.7;
        const bodiesToCheck = ['Sun', 'Mercury', 'Venus', 'Earth-Moon Barycenter', 'Mars', 'Jupiter'];
        bodiesToCheck.forEach(body => {
          if (firstSample.positions_km[body]) {
            const pos_km = firstSample.positions_km[body];
            // Positions are arrays [x, y, z] not objects
            const x_au = pos_km[0] / AU_TO_KM;
            const y_au = pos_km[1] / AU_TO_KM;
            const z_au = pos_km[2] / AU_TO_KM;
            const dist_au = Math.sqrt(x_au*x_au + y_au*y_au + z_au*z_au);
            console.log(`  ${body}: distance from origin = ${dist_au.toFixed(3)} AU`);
          }
        });
        
        // Check a middle sample too
        const midIndex = Math.floor(ephemerisData.samples.length / 2);
        const midSample = ephemerisData.samples[midIndex];
        const midJD = midSample.julian_date;
        const midDateObj = new Date((midJD - 2440587.5) * 86400000);
        console.log(`\n📊 DATA QUALITY CHECK - Middle Sample (${midDateObj.toISOString().split('T')[0]}):`);
        bodiesToCheck.forEach(body => {
          if (midSample.positions_km[body]) {
            const pos_km = midSample.positions_km[body];
            // Positions are arrays [x, y, z] not objects
            const x_au = pos_km[0] / AU_TO_KM;
            const y_au = pos_km[1] / AU_TO_KM;
            const z_au = pos_km[2] / AU_TO_KM;
            const dist_au = Math.sqrt(x_au*x_au + y_au*y_au + z_au*z_au);
            console.log(`  ${body}: distance from origin = ${dist_au.toFixed(3)} AU`);
          }
        });
        console.log('\n');
      }
      
      // Set default date to middle of available range
      const midJD = (ephemerisData.metadata.start_julian_date + ephemerisData.metadata.end_julian_date) / 2;
      const midDate = new Date((midJD - 2440587.5) * 86400000);
      birthDate.value = midDate.toISOString().split('T')[0];
      
      // Set date input constraints
      birthDate.min = startDateStr;
      birthDate.max = endDateStr;
      birthDate.title = `Available range: ${startDateStr} to ${endDateStr}`;
      
      // Show date range hint
      if (dateRangeHint) {
        dateRangeHint.textContent = `📅 Available: ${startDateStr} to ${endDateStr}`;
        dateRangeHint.style.display = 'block';
      }
      
      btnCalculate.disabled = false;
      
      // Show date range hint
      if (dateRangeHint) {
        dateRangeHint.textContent = `📅 Available: ${startDateStr} to ${endDateStr}`;
        dateRangeHint.style.display = 'block';
      }
      
      btnCalculate.disabled = false;
    } catch (error) {
      console.error('✗ Error loading ephemeris:', error);
      alert('Failed to load and compute ephemeris data from .eph file. Error: ' + error.message);
    }
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
  const scaleX = canvas.width / rect.width;
  const scaleY = canvas.height / rect.height;
  return {
    x: (evt.clientX - rect.left) * scaleX,
    y: (evt.clientY - rect.top) * scaleY
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

// Note: Default date is set in loadEphemerisData() based on available data range
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

// ============================================================================
// Toast Notification System
// ============================================================================

function showToast(message, type = 'info', duration = 4000) {
  const container = document.getElementById('toast-container');
  if (!container) return;
  
  const toast = document.createElement('div');
  toast.className = `toast ${type}`;
  
  const icons = {
    success: '✓',
    error: '✕',
    warning: '⚠',
    info: 'ℹ'
  };
  
  toast.innerHTML = `
    <div class="toast-icon">${icons[type] || icons.info}</div>
    <div class="toast-content">
      <div class="toast-message">${message}</div>
    </div>
    <button class="toast-close" aria-label="Close">×</button>
  `;
  
  container.appendChild(toast);
  
  const closeBtn = toast.querySelector('.toast-close');
  closeBtn.addEventListener('click', () => {
    toast.remove();
  });
  
  if (duration > 0) {
    setTimeout(() => {
      toast.style.opacity = '0';
      setTimeout(() => toast.remove(), 300);
    }, duration);
  }
}

// ============================================================================
// Starmatch Mode (Comparison)
// ============================================================================

const btnChartMode = document.getElementById('btn-chart-mode');
const btnStarmatchMode = document.getElementById('btn-starmatch-mode');
const starmatchSection = document.getElementById('starmatch-section');
const chartInputControls = document.querySelector('.input-controls');
const chartVisualisation = document.querySelector('.chart-visualisation');
const resultsContainer = document.querySelector('.results-container');
const analysisDetails = document.querySelector('.analysis-details');

const subjectSelect = document.getElementById('subject-select');
const targetSelect = document.getElementById('target-select');
const btnLoadSubject = document.getElementById('btn-load-subject');
const btnLoadTarget = document.getElementById('btn-load-target');
const btnCompare = document.getElementById('btn-compare');
const subjectInfo = document.getElementById('subject-info');
const targetInfo = document.getElementById('target-info');
const comparisonResults = document.getElementById('comparison-results');
const comparisonOutput = document.getElementById('comparison-output');

let currentSubject = null;
let currentTarget = null;

// Mode switching
function switchToChartMode() {
  btnChartMode.classList.add('active');
  btnStarmatchMode.classList.remove('active');
  
  starmatchSection.classList.add('hidden');
  chartInputControls.style.display = 'grid';
  chartVisualisation.style.display = 'block';
  resultsContainer.style.display = 'grid';
  analysisDetails.style.display = 'block';
}

function switchToStarmatchMode() {
  const records = loadRecords();
  
  if (records.length < 2) {
    showToast('Please create at least 2 records before using Starmatch mode.', 'warning', 5000);
    return;
  }
  
  btnStarmatchMode.classList.add('active');
  btnChartMode.classList.remove('active');
  
  starmatchSection.classList.remove('hidden');
  chartInputControls.style.display = 'none';
  chartVisualisation.style.display = 'none';
  resultsContainer.style.display = 'none';
  analysisDetails.style.display = 'none';
  
  populateComparisonSelects();
}

btnChartMode?.addEventListener('click', switchToChartMode);
btnStarmatchMode?.addEventListener('click', switchToStarmatchMode);

// Populate dropdowns with saved records
function populateComparisonSelects() {
  const records = loadRecords();
  
  subjectSelect.innerHTML = '<option value="">-- Select Subject --</option>';
  targetSelect.innerHTML = '<option value="">-- Select Target --</option>';
  
  records.forEach(rec => {
    const optionSubject = document.createElement('option');
    optionSubject.value = rec.id;
    optionSubject.textContent = rec.name;
    subjectSelect.appendChild(optionSubject);
    
    const optionTarget = document.createElement('option');
    optionTarget.value = rec.id;
    optionTarget.textContent = rec.name;
    targetSelect.appendChild(optionTarget);
  });
}

// Load subject
function loadSubjectForComparison() {
  const selectedId = subjectSelect.value;
  if (!selectedId) {
    showToast('Please select a subject from the dropdown.', 'warning');
    return;
  }
  
  const records = loadRecords();
  const record = records.find(r => r.id === selectedId);
  
  if (!record) {
    showToast('Subject record not found.', 'error');
    return;
  }
  
  // Check if same as target
  if (currentTarget && currentTarget.id === record.id) {
    showToast('Subject and Target cannot be the same person.', 'error');
    return;
  }
  
  currentSubject = record;
  displayPersonInfo(record, subjectInfo);
  updateCompareButton();
  showToast(`Loaded subject: ${record.name}`, 'success', 2000);
}

// Load target
function loadTargetForComparison() {
  const selectedId = targetSelect.value;
  if (!selectedId) {
    showToast('Please select a target from the dropdown.', 'warning');
    return;
  }
  
  const records = loadRecords();
  const record = records.find(r => r.id === selectedId);
  
  if (!record) {
    showToast('Target record not found.', 'error');
    return;
  }
  
  // Check if same as subject
  if (currentSubject && currentSubject.id === record.id) {
    showToast('Subject and Target cannot be the same person.', 'error');
    return;
  }
  
  currentTarget = record;
  displayPersonInfo(record, targetInfo);
  updateCompareButton();
  showToast(`Loaded target: ${record.name}`, 'success', 2000);
}

// Display person info
function displayPersonInfo(record, container) {
  container.innerHTML = `
    <strong>${record.name}</strong><br>
    Date: ${record.date || 'N/A'}<br>
    Time: ${record.time || 'N/A'}<br>
    Lat: ${record.lat || 'N/A'}, Lon: ${record.lon || 'N/A'}<br>
    <em style="opacity: 0.7; font-size: 0.8rem;">Settings: Orb ${record.orbType || 0}, Asp ${record.aspectOrbSet || 0}, Rule ${record.rulershipSet || 0}</em>
  `;
}

// Update compare button state
function updateCompareButton() {
  if (currentSubject && currentTarget) {
    btnCompare.disabled = false;
  } else {
    btnCompare.disabled = true;
  }
}

// Perform comparison
function performComparison() {
  if (!currentSubject || !currentTarget) {
    showToast('Please load both Subject and Target.', 'warning');
    return;
  }
  
  //showToast('Calculating comparison...', 'info', 2000);
  
  // Calculate charts for both
  try {
    // Subject chart
    const subjectJD = dateToJulianDate(currentSubject.date, currentSubject.time);
    const subjectPositions = getPositionsForJD(subjectJD);
    const subjectPlanetaryPositions = extractPlanetaryPositions(subjectPositions);
    const subjectAsc = calculateAscendant(subjectJD, parseFloat(currentSubject.lat), parseFloat(currentSubject.lon));
    const subjectMC = calculateMidheaven(subjectJD, parseFloat(currentSubject.lon));
    
    // Target chart
    const targetJD = dateToJulianDate(currentTarget.date, currentTarget.time);
    const targetPositions = getPositionsForJD(targetJD);
    const targetPlanetaryPositions = extractPlanetaryPositions(targetPositions);
    const targetAsc = calculateAscendant(targetJD, parseFloat(currentTarget.lat), parseFloat(currentTarget.lon));
    const targetMC = calculateMidheaven(targetJD, parseFloat(currentTarget.lon));
    
    // Apply settings for subject
    orbType = parseInt(currentSubject.orbType || 0);
    aoIndex = parseInt(currentSubject.aspectOrbSet || 0);
    tfIndex = parseInt(currentSubject.rulershipSet || 0);
    precessionFlag = currentSubject.precession || 0;
    const subjectYear = new Date(currentSubject.date).getFullYear();
    nativityYear = subjectYear;
    nativity = subjectYear;
    
    // Calculate subject themes
    getThemeValues(
      subjectPlanetaryPositions.Sun,
      subjectPlanetaryPositions.Moon,
      subjectPlanetaryPositions.Mercury,
      subjectPlanetaryPositions.Venus,
      subjectPlanetaryPositions.Mars,
      subjectPlanetaryPositions.Jupiter,
      subjectPlanetaryPositions.Saturn,
      subjectPlanetaryPositions.Uranus,
      subjectPlanetaryPositions.Neptune,
      subjectPlanetaryPositions.Pluto,
      subjectAsc,
      subjectMC
    );
    
    const subjectThemes = [...theme]; // Copy theme values
    
    // Apply settings for target
    orbType = parseInt(currentTarget.orbType || 0);
    aoIndex = parseInt(currentTarget.aspectOrbSet || 0);
    tfIndex = parseInt(currentTarget.rulershipSet || 0);
    precessionFlag = currentTarget.precession || 0;
    const targetYear = new Date(currentTarget.date).getFullYear();
    nativityYear = targetYear;
    nativity = targetYear;
    
    // Calculate target themes
    getThemeValues(
      targetPlanetaryPositions.Sun,
      targetPlanetaryPositions.Moon,
      targetPlanetaryPositions.Mercury,
      targetPlanetaryPositions.Venus,
      targetPlanetaryPositions.Mars,
      targetPlanetaryPositions.Jupiter,
      targetPlanetaryPositions.Saturn,
      targetPlanetaryPositions.Uranus,
      targetPlanetaryPositions.Neptune,
      targetPlanetaryPositions.Pluto,
      targetAsc,
      targetMC
    );
    
    const targetThemes = [...theme]; // Copy theme values
    
    // Display comparison results with planetary positions
    displayComparisonResults(
      subjectThemes, 
      targetThemes, 
      subjectPlanetaryPositions, 
      targetPlanetaryPositions,
      subjectAsc,
      targetAsc
    );
    
    comparisonResults.classList.remove('hidden');
    //showToast('Comparison complete!', 'success', 3000);
    
  } catch (error) {
    console.error('Comparison error:', error);
    showToast('Error calculating comparison: ' + error.message, 'error', 5000);
  }
}

// Extract planetary positions helper
function extractPlanetaryPositions(positions) {
  const bodyMapping = {
    'Sun': 'Sun',
    'Moon': 'Moon',  // DE200 has separate Moon ephemeris
    'Mercury': 'Mercury',
    'Venus': 'Venus',
    'Mars': 'Mars',
    'Jupiter': 'Jupiter',
    'Saturn': 'Saturn',
    'Uranus': 'Uranus',
    'Neptune': 'Neptune',
    'Pluto': 'Pluto'
  };
  
  const planetaryPositions = {};
  
  for (const [ourName, ephemName] of Object.entries(bodyMapping)) {
    if (positions[ephemName]) {
      const [x, y, z] = positions[ephemName];
      let longitude;
      
      // Special case for Sun: In heliocentric coords, Sun is at origin (0,0,0)
      // So Sun's apparent position from Earth is Earth's position + 180°
      if (ourName === 'Sun') {
        // Use Earth-Moon Barycenter as proxy for Earth's position
        const earthPos = positions['Earth-Moon Barycenter'];
        if (earthPos) {
          // Sun as seen from Earth is opposite direction from Earth's heliocentric position
          longitude = cartesianToLongitude(-earthPos[0], -earthPos[1], -earthPos[2]);
        } else {
          console.warn('Earth-Moon Barycenter position not found for Sun calculation');
          longitude = 0;
        }
      } else {
        longitude = cartesianToLongitude(x, y, z);
      }
      
      planetaryPositions[ourName] = longitude;
    } else {
      planetaryPositions[ourName] = 0;
    }
  }
  
  return planetaryPositions;
}

// Calculate xProfile value (similarity-complementarity spectrum)
// Returns value from -1 (complementarity/inverted) to +1 (similarity/same shape)
// Values near 0 indicate balanced relationships (most significant/lasting)
function calculateXProfileValue(subjectThemes, targetThemes) {
  // Calculate correlation coefficient between the two theme arrays
  const n = subjectThemes.length;
  
  // Calculate means
  const meanSubject = subjectThemes.reduce((a, b) => a + b, 0) / n;
  const meanTarget = targetThemes.reduce((a, b) => a + b, 0) / n;
  
  // Calculate correlation
  let numerator = 0;
  let denomSubject = 0;
  let denomTarget = 0;
  
  for (let i = 0; i < n; i++) {
    const diffSubject = subjectThemes[i] - meanSubject;
    const diffTarget = targetThemes[i] - meanTarget;
    numerator += diffSubject * diffTarget;
    denomSubject += diffSubject * diffSubject;
    denomTarget += diffTarget * diffTarget;
  }
  
  const correlation = numerator / Math.sqrt(denomSubject * denomTarget);
  
  // Correlation ranges from -1 to +1
  // +1 = perfect positive correlation (same shape) = similarity
  // -1 = perfect negative correlation (inverted shape) = complementarity
  // 0 = no correlation = equality/balance
  
  return correlation;
}

// Get relationship type interpretation based on xProfile value
function getRelationshipTypeInterpretation(xProfileValue) {
  const absValue = Math.abs(xProfileValue);
  
  if (absValue < 0.2) {
    return {
      type: 'Equality/Balance',
      color: '#51cf66',
      description: 'Optimal for long-lasting, significant relationships. A balanced blend of similarity and complementarity.',
      significance: 'High'
    };
  } else if (xProfileValue > 0.7) {
    return {
      type: 'Strong Similarity',
      color: '#74c0fc',
      description: 'Charts have the same shape. Good relationship potential, though may lack the balance for deepest partnerships.',
      significance: 'Moderate'
    };
  } else if (xProfileValue > 0.4) {
    return {
      type: 'Moderate Similarity',
      color: '#69db7c',
      description: 'Similar energies with some variation. Good compatibility with room for growth.',
      significance: 'Good'
    };
  } else if (xProfileValue < -0.7) {
    return {
      type: 'Strong Complementarity',
      color: '#b85eff',
      description: 'Charts are inverted relative to each other. Complementary energies, though may lack balance for lasting partnerships.',
      significance: 'Moderate'
    };
  } else if (xProfileValue < -0.4) {
    return {
      type: 'Moderate Complementarity',
      color: '#a78bfa',
      description: 'Complementary energies provide contrast and growth opportunities.',
      significance: 'Good'
    };
  } else {
    return {
      type: 'Mixed Balance',
      color: '#ffd43b',
      description: 'A mixture of similar and complementary energies. Approaching ideal balance.',
      significance: 'Good'
    };
  }
}

// Display comparison results
function displayComparisonResults(subjectThemes, targetThemes, subjectPos, targetPos, subjectAsc, targetAsc) {
  const SIGN_NAMES = ['Aries', 'Taurus', 'Gemini', 'Cancer', 'Leo', 'Virgo', 
                      'Libra', 'Scorpio', 'Sagittarius', 'Capricorn', 'Aquarius', 'Pisces'];
  
  // Calculate xProfile value
  const xProfileValue = calculateXProfileValue(subjectThemes, targetThemes);
  const relType = getRelationshipTypeInterpretation(xProfileValue);
  
  let html = '<div class="comparison-grid">';
  
  // xProfile Spectrum Display - wrapped in its own container
  html += `<div class="xprofile-spectrum-container">
    <h4 style="color: var(--accent); margin-top: 0;">xProfile Relationship Spectrum</h4>
    <div style="background: rgba(94,197,255,0.1); padding: 1.5rem; border-radius: 8px; border: 1px solid rgba(94,197,255,0.3);">
      
      <!-- Spectrum Bar -->
      <div style="margin-bottom: 1.5rem;">
        <div style="display: flex; justify-content: space-between; font-size: 0.7rem; color: #8fa8ce; margin-bottom: 0.5rem;">
          <span>Complementarity</span>
          <span>Equality</span>
          <span>Similarity</span>
        </div>
        <div style="position: relative; height: 30px; background: linear-gradient(90deg, #b85eff 0%, #ffd43b 50%, #74c0fc 100%); border-radius: 6px; border: 1px solid rgba(94,197,255,0.3);">
          <!-- Marker -->
          <div style="position: absolute; left: ${((xProfileValue + 1) / 2) * 100}%; top: -5px; transform: translateX(-50%);">
            <div style="width: 0; height: 0; border-left: 8px solid transparent; border-right: 8px solid transparent; border-top: 10px solid white;"></div>
          </div>
          <!-- Value marker line -->
          <div style="position: absolute; left: ${((xProfileValue + 1) / 2) * 100}%; top: 0; bottom: 0; width: 2px; background: white; transform: translateX(-50%);"></div>
        </div>
        <div style="display: flex; justify-content: space-between; font-size: 0.65rem; color: #6a7fa0; margin-top: 0.25rem;">
          <span>-1.0</span>
          <span>0.0</span>
          <span>+1.0</span>
        </div>
      </div>
      
      <!-- xProfile Value -->
      <div style="text-align: center; margin-bottom: 1rem;">
        <div style="font-size: 0.85rem; color: #8fa8ce; margin-bottom: 0.5rem;">xProfile Value</div>
        <div style="font-size: 3rem; font-weight: 700; color: ${relType.color};">${xProfileValue.toFixed(3)}</div>
      </div>
      
      <!-- Relationship Type -->
      <div style="background: rgba(0,0,0,0.3); padding: 1rem; border-radius: 6px; border-left: 4px solid ${relType.color};">
        <div style="font-size: 1.1rem; font-weight: 600; color: ${relType.color}; margin-bottom: 0.5rem;">${relType.type}</div>
        <div style="font-size: 0.85rem; line-height: 1.6; color: #b8d0f0; margin-bottom: 0.75rem;">${relType.description}</div>
        <div style="font-size: 0.75rem; color: #8fa8ce;">
          <strong>Relationship Significance:</strong> ${relType.significance}
        </div>
      </div>
      
      ${Math.abs(xProfileValue) < 0.2 ? 
        '<div style="margin-top: 1rem; padding: 0.75rem; background: rgba(81,207,102,0.15); border-radius: 6px; border: 1px solid rgba(81,207,102,0.3); font-size: 0.8rem; color: #51cf66;">★ Optimal balance for long-lasting partnerships</div>' : 
        Math.abs(xProfileValue) > 0.7 ?
        '<div style="margin-top: 1rem; padding: 0.75rem; background: rgba(255,212,59,0.15); border-radius: 6px; border: 1px solid rgba(255,212,59,0.3); font-size: 0.8rem; color: #ffd43b;">⚠ Extreme values suggest good relationships but less likely for deep partnerships</div>' :
        ''}
    </div>
  </div>`;
  
  // Theme comparison - wrapped in its own container
  html += '<div class="theme-comparison-container">';
  html += '<h4 style="color: var(--accent); margin-top: 0;">Theme-by-Theme Analysis</h4>';
  html += '<div style="display: flex; flex-direction: column; gap: 0.6rem;">';
  
  // Find max value for scaling
  const maxTheme = Math.max(...subjectThemes, ...targetThemes);
  
  for (let i = 0; i < 12; i++) {
    const subjectVal = subjectThemes[i];
    const targetVal = targetThemes[i];
    const subjectPercent = (subjectVal / maxTheme) * 100;
    const targetPercent = (targetVal / maxTheme) * 100;
    const diff = Math.abs(subjectVal - targetVal);
    
    html += `
      <div style="display: flex; align-items: center; gap: 0.5rem; font-size: 0.75rem;">
        <div style="min-width: 70px; color: #b8d0f0; text-align: right; font-weight: 500;">${SIGN_NAMES[i]}</div>
        
        <!-- Subject bar (left side, blue) -->
        <div style="flex: 1; display: flex; justify-content: flex-end; align-items: center; gap: 0.3rem;">
          <div style="font-family: 'Fira Code', monospace; font-size: 0.7rem; color: #74c0fc; min-width: 35px; text-align: right;">${subjectVal.toFixed(1)}</div>
          <div style="width: 100%; height: 20px; background: rgba(10,13,19,0.8); border-radius: 3px; overflow: hidden; border: 1px solid rgba(116,197,252,0.3); position: relative;">
            <div style="position: absolute; right: 0; height: 100%; width: ${subjectPercent}%; background: linear-gradient(90deg, rgba(116,197,252,0.3), #74c0fc); transition: width 0.6s;"></div>
          </div>
        </div>
        
        <!-- Target bar (right side, purple) -->
        <div style="flex: 1; display: flex; align-items: center; gap: 0.3rem;">
          <div style="width: 100%; height: 20px; background: rgba(10,13,19,0.8); border-radius: 3px; overflow: hidden; border: 1px solid rgba(184,94,255,0.3); position: relative;">
            <div style="position: absolute; left: 0; height: 100%; width: ${targetPercent}%; background: linear-gradient(90deg, #b85eff, rgba(184,94,255,0.3)); transition: width 0.6s;"></div>
          </div>
          <div style="font-family: 'Fira Code', monospace; font-size: 0.7rem; color: #b85eff; min-width: 35px;">${targetVal.toFixed(1)}</div>
        </div>
        
        <!-- Difference indicator -->
        <div style="min-width: 40px; text-align: center; font-size: 0.65rem; color: ${diff < 2 ? '#51cf66' : diff < 5 ? '#ffd43b' : '#ff6b6b'}; font-family: 'Fira Code', monospace;">
          Δ${diff.toFixed(1)}
        </div>
      </div>
    `;
  }
  
  html += '</div>';
  html += '<div style="margin-top: 1rem; padding-top: 1rem; border-top: 1px solid rgba(94,197,255,0.15); font-size: 0.7rem; color: #8fa8ce; display: flex; justify-content: space-between; align-items: center;">';
  html += '<div style="display: flex; gap: 1.5rem;">';
  html += '<div><span style="color: #74c0fc;">━━━</span> Subject</div>';
  html += '<div><span style="color: #b85eff;">━━━</span> Target</div>';
  html += '</div>';
  html += '<div style="font-style: italic;">Δ = Difference</div>';
  html += '</div>';
  html += '</div>';
  
  html += '</div>';
  
  // Bottom section with chart visualization
  html += `<div class="comparison-bottom-grid">
    <div style="padding: 1rem; background: rgba(10,13,19,0.6); border-radius: 8px; border: 1px solid rgba(94,197,255,0.15);">
      <div style="font-size: 0.75rem; color: #8fa8ce; line-height: 1.6;">
        <strong style="color: #b8d0f0;">Understanding xProfile Values:</strong><br>
        <span style="color: #74c0fc;">+1.0</span> = Charts have same shape (similarity)<br>
        <span style="color: #ffd43b;">0.0</span> = Perfect balance (ideal for lasting relationships)<br>
        <span style="color: #b85eff;">-1.0</span> = Charts are inverted (complementarity)
      </div>
      <div style="margin-top: 1rem; font-size: 0.7rem; color: #6a7fa0; font-style: italic;">
        Subject: ${currentSubject.name} • Target: ${currentTarget.name}
      </div>
    </div>
    <div style="padding: 1rem; background: rgba(10,13,19,0.6); border-radius: 8px; border: 1px solid rgba(94,197,255,0.15);">
      <h4 style="color: var(--accent); margin-top: 0; margin-bottom: 0.75rem; font-size: 0.9rem;">Chart Overlay</h4>
      <canvas id="comparison-chart-canvas" width="400" height="400" style="width: 100%; max-width: 400px; aspect-ratio: 1/1; display: block; margin: 0 auto;"></canvas>
      <div id="comparison-tooltip" class="chart-tooltip"></div>
      <div style="margin-top: 0.75rem; font-size: 0.7rem; color: #8fa8ce; display: flex; justify-content: center; gap: 1.5rem;">
        <div><span style="color: #74c0fc;">●</span> Subject (${currentSubject.name})</div>
        <div><span style="color: #b85eff;">●</span> Target (${currentTarget.name})</div>
      </div>
    </div>
  </div>`;
  
  comparisonOutput.innerHTML = html;
  
  // Draw the comparison chart after the HTML is rendered
  setTimeout(() => {
    const canvas = document.getElementById('comparison-chart-canvas');
    if (canvas) {
      // Force square aspect ratio
      const rect = canvas.getBoundingClientRect();
      canvas.style.height = rect.width + 'px';
    }
    drawComparisonChart(subjectPos, targetPos, subjectAsc, targetAsc);
    setupComparisonChartTooltips(subjectPos, targetPos, subjectAsc, targetAsc);
  }, 50);
}

// Draw overlaid comparison chart
function drawComparisonChart(subjectPos, targetPos, subjectAsc, targetAsc) {
  const canvas = document.getElementById('comparison-chart-canvas');
  if (!canvas) return;
  
  const ctx = canvas.getContext('2d');
  const centerX = canvas.width / 2;
  const centerY = canvas.height / 2;
  const outerRadius = 180;
  const innerRadius = 140;
  const subjectPlanetRadius = 120;
  const targetPlanetRadius = 95;
  
  // Clear canvas
  ctx.clearRect(0, 0, canvas.width, canvas.height);
  ctx.fillStyle = '#05070f';
  ctx.fillRect(0, 0, canvas.width, canvas.height);
  
  const SIGN_NAMES = ['Aries', 'Taurus', 'Gemini', 'Cancer', 'Leo', 'Virgo', 
                      'Libra', 'Scorpio', 'Sagittarius', 'Capricorn', 'Aquarius', 'Pisces'];
  
  // Draw zodiac wheel for subject (using subject's ascendant)
  const signColors = [
    '#ff6b6b', '#51cf66', '#ffd43b', '#74c0fc',
    '#ff8787', '#69db7c', '#ffd43b', '#ff6b6b',
    '#cc5de8', '#51cf66', '#74c0fc', '#a78bfa'
  ];
  
  for (let i = 0; i < 12; i++) {
    const startAngle = ((i * 30 - subjectAsc - 90) * Math.PI) / 180;
    const endAngle = (((i + 1) * 30 - subjectAsc - 90) * Math.PI) / 180;
    
    // Draw sign segment
    ctx.beginPath();
    ctx.moveTo(centerX, centerY);
    ctx.arc(centerX, centerY, outerRadius, startAngle, endAngle);
    ctx.closePath();
    ctx.fillStyle = signColors[i] + '20';
    ctx.fill();
    ctx.strokeStyle = signColors[i] + '60';
    ctx.lineWidth = 1;
    ctx.stroke();
    
    // Draw sign name
    const midAngle = startAngle + (endAngle - startAngle) / 2;
    const textRadius = (outerRadius + innerRadius) / 2 + 10;
    const textX = centerX + Math.cos(midAngle) * textRadius;
    const textY = centerY + Math.sin(midAngle) * textRadius;
    
    ctx.save();
    ctx.translate(textX, textY);
    ctx.rotate(midAngle + Math.PI / 2);
    ctx.fillStyle = signColors[i];
    ctx.font = 'bold 10px "Segoe UI"';
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
  
  // Draw subject planets (blue) on outer ring
  const PLANET_SYMBOLS = ['☉', '☽', '☿', '♀', '♂', '♃', '♄', '⛢', '♆', '♇'];
  const PLANET_NAMES = ['Sun', 'Moon', 'Mercury', 'Venus', 'Mars', 'Jupiter', 
                        'Saturn', 'Uranus', 'Neptune', 'Pluto'];
  
  Object.entries(subjectPos).forEach(([name, longitude]) => {
    if (longitude === undefined || longitude === null || isNaN(longitude)) return;
    
    const angle = ((-longitude - 90) * Math.PI) / 180;
    const x = centerX + Math.cos(angle) * subjectPlanetRadius;
    const y = centerY + Math.sin(angle) * subjectPlanetRadius;
    
    const planetIndex = PLANET_NAMES.indexOf(name);
    if (planetIndex === -1) return;
    
    // Draw planet circle
    ctx.beginPath();
    ctx.arc(x, y, 10, 0, Math.PI * 2);
    ctx.fillStyle = '#74c0fc';
    ctx.fill();
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.5)';
    ctx.lineWidth = 2;
    ctx.stroke();
    
    // Draw planet symbol
    ctx.fillStyle = '#000';
    ctx.font = 'bold 14px Arial';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    ctx.fillText(PLANET_SYMBOLS[planetIndex], x, y);
  });
  
  // Draw target planets (purple) on inner ring
  Object.entries(targetPos).forEach(([name, longitude]) => {
    if (longitude === undefined || longitude === null || isNaN(longitude)) return;
    
    const angle = ((-longitude - 90) * Math.PI) / 180;
    const x = centerX + Math.cos(angle) * targetPlanetRadius;
    const y = centerY + Math.sin(angle) * targetPlanetRadius;
    
    const planetIndex = PLANET_NAMES.indexOf(name);
    if (planetIndex === -1) return;
    
    // Draw planet circle
    ctx.beginPath();
    ctx.arc(x, y, 10, 0, Math.PI * 2);
    ctx.fillStyle = '#b85eff';
    ctx.fill();
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.5)';
    ctx.lineWidth = 2;
    ctx.stroke();
    
    // Draw planet symbol
    ctx.fillStyle = '#000';
    ctx.font = 'bold 14px Arial';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    ctx.fillText(PLANET_SYMBOLS[planetIndex], x, y);
  });
  
  // Draw ascendant lines
  // Subject ascendant (blue) - solid line
  const subjectAscAngle = ((-subjectAsc - 90) * Math.PI) / 180;
  ctx.beginPath();
  ctx.moveTo(centerX, centerY);
  ctx.lineTo(
    centerX + Math.cos(subjectAscAngle) * innerRadius,
    centerY + Math.sin(subjectAscAngle) * innerRadius
  );
  ctx.strokeStyle = '#74c0fc';
  ctx.lineWidth = 3;
  ctx.stroke();
  
  // Label for subject ASC
  ctx.fillStyle = '#74c0fc';
  ctx.font = 'bold 10px "Segoe UI"';
  ctx.textAlign = 'center';
  const subjectLabelX = centerX + Math.cos(subjectAscAngle) * (innerRadius - 15);
  const subjectLabelY = centerY + Math.sin(subjectAscAngle) * (innerRadius - 15);
  ctx.fillText('ASC', subjectLabelX, subjectLabelY);
  
  // Target ascendant (purple) - dashed line
  const targetAscAngle = ((-targetAsc - 90) * Math.PI) / 180;
  ctx.beginPath();
  ctx.moveTo(centerX, centerY);
  ctx.lineTo(
    centerX + Math.cos(targetAscAngle) * innerRadius,
    centerY + Math.sin(targetAscAngle) * innerRadius
  );
  ctx.strokeStyle = '#b85eff';
  ctx.lineWidth = 3;
  ctx.setLineDash([5, 5]);
  ctx.stroke();
  ctx.setLineDash([]);
  
  // Label for target ASC
  ctx.fillStyle = '#b85eff';
  ctx.font = 'bold 10px "Segoe UI"';
  ctx.textAlign = 'center';
  const targetLabelX = centerX + Math.cos(targetAscAngle) * (innerRadius - 30);
  const targetLabelY = centerY + Math.sin(targetAscAngle) * (innerRadius - 30);
  ctx.fillText('ASC', targetLabelX, targetLabelY);
}

// Store event listeners for cleanup
let comparisonTooltipListeners = null;

// Setup tooltips for comparison chart
function setupComparisonChartTooltips(subjectPos, targetPos, subjectAsc, targetAsc) {
  const canvas = document.getElementById('comparison-chart-canvas');
  if (!canvas) {
    console.error('Comparison canvas not found');
    return;
  }
  
  const tooltip = document.getElementById('comparison-tooltip');
  if (!tooltip) {
    console.error('Comparison tooltip element not found');
    return;
  }
  
  // Remove old event listeners if they exist
  if (comparisonTooltipListeners) {
    canvas.removeEventListener('mousemove', comparisonTooltipListeners.mousemove);
    canvas.removeEventListener('mouseleave', comparisonTooltipListeners.mouseleave);
  }
  
  const SIGN_NAMES = ['Aries', 'Taurus', 'Gemini', 'Cancer', 'Leo', 'Virgo', 
                      'Libra', 'Scorpio', 'Sagittarius', 'Capricorn', 'Aquarius', 'Pisces'];
  const PLANET_NAMES = ['Sun', 'Moon', 'Mercury', 'Venus', 'Mars', 'Jupiter', 
                        'Saturn', 'Uranus', 'Neptune', 'Pluto'];
  
  // Store chart data for tooltip calculations
  const chartData = {
    centerX: canvas.width / 2,
    centerY: canvas.height / 2,
    outerRadius: 180,
    innerRadius: 140,
    subjectPlanetRadius: 120,
    targetPlanetRadius: 95,
    subjectPos,
    targetPos,
    subjectAsc,
    targetAsc
  };
  
  function getMousePos(canvas, evt) {
    const rect = canvas.getBoundingClientRect();
    const scaleX = canvas.width / rect.width;
    const scaleY = canvas.height / rect.height;
    return {
      x: (evt.clientX - rect.left) * scaleX,
      y: (evt.clientY - rect.top) * scaleY
    };
  }
  
  function checkPlanetHover(mouseX, mouseY) {
    const planetRadius = 15; // Increased hit detection radius for easier hovering
    
    // Check subject planets
    for (const [name, longitude] of Object.entries(chartData.subjectPos)) {
      if (isNaN(longitude) || !isFinite(longitude)) continue;
      
      const angle = ((-longitude - 90) * Math.PI) / 180;
      const x = chartData.centerX + Math.cos(angle) * chartData.subjectPlanetRadius;
      const y = chartData.centerY + Math.sin(angle) * chartData.subjectPlanetRadius;
      
      const distance = Math.sqrt((mouseX - x) ** 2 + (mouseY - y) ** 2);
      
      if (distance <= planetRadius) {
        let normalizedLon = longitude % 360;
        if (normalizedLon < 0) normalizedLon += 360;
        const sign = Math.floor(normalizedLon / 30);
        const degree = normalizedLon % 30;
        const signName = SIGN_NAMES[sign];
        
        return {
          type: 'subject-planet',
          name: name,
          longitude: longitude,
          position: `${degree.toFixed(2)}° ${signName}`,
          sign: signName,
          person: currentSubject.name
        };
      }
    }
    
    // Check target planets
    for (const [name, longitude] of Object.entries(chartData.targetPos)) {
      if (isNaN(longitude) || !isFinite(longitude)) continue;
      
      const angle = ((-longitude - 90) * Math.PI) / 180;
      const x = chartData.centerX + Math.cos(angle) * chartData.targetPlanetRadius;
      const y = chartData.centerY + Math.sin(angle) * chartData.targetPlanetRadius;
      
      const distance = Math.sqrt((mouseX - x) ** 2 + (mouseY - y) ** 2);
      
      if (distance <= planetRadius) {
        let normalizedLon = longitude % 360;
        if (normalizedLon < 0) normalizedLon += 360;
        const sign = Math.floor(normalizedLon / 30);
        const degree = normalizedLon % 30;
        const signName = SIGN_NAMES[sign];
        
        return {
          type: 'target-planet',
          name: name,
          longitude: longitude,
          position: `${degree.toFixed(2)}° ${signName}`,
          sign: signName,
          person: currentTarget.name
        };
      }
    }
    
    return null;
  }
  
  function checkAscendantHover(mouseX, mouseY) {
    // Check subject ascendant line (blue)
    const subjectAscAngle = ((-chartData.subjectAsc - 90) * Math.PI) / 180;
    const subjectX2 = chartData.centerX + Math.cos(subjectAscAngle) * chartData.innerRadius;
    const subjectY2 = chartData.centerY + Math.sin(subjectAscAngle) * chartData.innerRadius;
    
    const distanceToSubjectAsc = distanceToLineSegment(
      mouseX, mouseY,
      chartData.centerX, chartData.centerY,
      subjectX2, subjectY2
    );
    
    if (distanceToSubjectAsc <= 5) {
      let normalizedLon = chartData.subjectAsc % 360;
      if (normalizedLon < 0) normalizedLon += 360;
      const sign = Math.floor(normalizedLon / 30);
      const degree = normalizedLon % 30;
      const signName = SIGN_NAMES[sign];
      
      return {
        type: 'subject-ascendant',
        person: currentSubject.name,
        position: `${degree.toFixed(2)}° ${signName}`,
        sign: signName
      };
    }
    
    // Check target ascendant line (purple)
    const targetAscAngle = ((-chartData.targetAsc - 90) * Math.PI) / 180;
    const targetX2 = chartData.centerX + Math.cos(targetAscAngle) * chartData.innerRadius;
    const targetY2 = chartData.centerY + Math.sin(targetAscAngle) * chartData.innerRadius;
    
    const distanceToTargetAsc = distanceToLineSegment(
      mouseX, mouseY,
      chartData.centerX, chartData.centerY,
      targetX2, targetY2
    );
    
    if (distanceToTargetAsc <= 5) {
      let normalizedLon = chartData.targetAsc % 360;
      if (normalizedLon < 0) normalizedLon += 360;
      const sign = Math.floor(normalizedLon / 30);
      const degree = normalizedLon % 30;
      const signName = SIGN_NAMES[sign];
      
      return {
        type: 'target-ascendant',
        person: currentTarget.name,
        position: `${degree.toFixed(2)}° ${signName}`,
        sign: signName
      };
    }
    
    return null;
  }
  
  function checkSignHover(mouseX, mouseY) {
    const ELEMENT_NAMES = ['Fire', 'Earth', 'Air', 'Water'];
    const QUALITY_NAMES = ['Cardinal', 'Fixed', 'Mutable'];
    
    const dx = mouseX - chartData.centerX;
    const dy = mouseY - chartData.centerY;
    const distance = Math.sqrt(dx * dx + dy * dy);
    
    // Check if in zodiac ring area
    if (distance >= chartData.innerRadius && distance <= chartData.outerRadius) {
      // Calculate angle from center
      // Mirror about vertical axis by negating dx instead of the angle
      let angle = Math.atan2(dy, -dx) * (180 / Math.PI);
      // Convert to zodiac longitude (adjusted for ascendant and 90° offset)
      let zodiacLon = -angle - 90 + chartData.subjectAsc;
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
  
  function updateTooltip(evt) {
    const mousePos = getMousePos(canvas, evt);
    const mouseX = mousePos.x;
    const mouseY = mousePos.y;
    
    // Check in priority order: planets, ascendants, signs
    let hoverInfo = checkPlanetHover(mouseX, mouseY);
    
    if (!hoverInfo) {
      hoverInfo = checkAscendantHover(mouseX, mouseY);
    }
    
    if (!hoverInfo) {
      hoverInfo = checkSignHover(mouseX, mouseY);
    }
    
    if (hoverInfo) {
      let tooltipHTML = '';
      
      if (hoverInfo.type === 'subject-planet') {
        tooltipHTML = `
          <div style="font-weight: 600; color: #74c0fc; margin-bottom: 0.25rem;">${hoverInfo.name} (${hoverInfo.person})</div>
          <div style="font-size: 0.85rem;">${hoverInfo.position}</div>
        `;
      } else if (hoverInfo.type === 'target-planet') {
        tooltipHTML = `
          <div style="font-weight: 600; color: #b85eff; margin-bottom: 0.25rem;">${hoverInfo.name} (${hoverInfo.person})</div>
          <div style="font-size: 0.85rem;">${hoverInfo.position}</div>
        `;
      } else if (hoverInfo.type === 'subject-ascendant') {
        tooltipHTML = `
          <div style="font-weight: 600; color: #74c0fc; margin-bottom: 0.25rem;">Ascendant (${hoverInfo.person})</div>
          <div style="font-size: 0.85rem;">${hoverInfo.position}</div>
        `;
      } else if (hoverInfo.type === 'target-ascendant') {
        tooltipHTML = `
          <div style="font-weight: 600; color: #b85eff; margin-bottom: 0.25rem;">Ascendant (${hoverInfo.person})</div>
          <div style="font-size: 0.85rem;">${hoverInfo.position}</div>
        `;
      } else if (hoverInfo.type === 'sign') {
        tooltipHTML = `
          <div style="font-weight: 600; color: var(--accent-warm); margin-bottom: 0.25rem;">${hoverInfo.name}</div>
          <div style="font-size: 0.85rem; color: #b8d0f0;">${hoverInfo.element} • ${hoverInfo.quality}</div>
        `;
      }
      
      tooltip.innerHTML = tooltipHTML;
      tooltip.style.display = 'block';
      
      // Position tooltip near cursor (fixed positioning uses viewport coordinates)
      tooltip.style.left = (evt.clientX + 15) + 'px';
      tooltip.style.top = (evt.clientY + 15) + 'px';
      
      canvas.style.cursor = 'pointer';
    } else {
      tooltip.style.display = 'none';
      canvas.style.cursor = 'default';
    }
  }
  
  // Define event listeners
  const mouseleaveHandler = () => {
    tooltip.style.display = 'none';
    canvas.style.cursor = 'default';
  };
  
  // Store references for cleanup
  comparisonTooltipListeners = {
    mousemove: updateTooltip,
    mouseleave: mouseleaveHandler
  };
  
  // Add event listeners
  canvas.addEventListener('mousemove', updateTooltip);
  canvas.addEventListener('mouseleave', mouseleaveHandler);
}

// Helper function to draw planets for comparison (no longer needed but keeping for compatibility)
function drawComparisonPlanets(ctx, centerX, centerY, radius, positions, color, size) {
  Object.entries(positions).forEach(([name, longitude]) => {
    if (longitude === undefined || longitude === null || isNaN(longitude)) return;
    
    const angle = ((-longitude - 90) * Math.PI) / 180;
    const x = centerX + Math.cos(angle) * radius;
    const y = centerY + Math.sin(angle) * radius;
    
    // Draw planet circle
    ctx.beginPath();
    ctx.arc(x, y, size, 0, Math.PI * 2);
    ctx.fillStyle = color;
    ctx.fill();
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.3)';
    ctx.lineWidth = 1;
    ctx.stroke();
  });
}

// Event listeners
btnLoadSubject?.addEventListener('click', loadSubjectForComparison);
btnLoadTarget?.addEventListener('click', loadTargetForComparison);
btnCompare?.addEventListener('click', performComparison);

// Prevent same person selection
subjectSelect?.addEventListener('change', () => {
  if (currentTarget && subjectSelect.value === currentTarget.id) {
    showToast('Subject and Target cannot be the same person.', 'warning');
    subjectSelect.value = '';
  }
});

targetSelect?.addEventListener('change', () => {
  if (currentSubject && targetSelect.value === currentSubject.id) {
    showToast('Subject and Target cannot be the same person.', 'warning');
    targetSelect.value = '';
  }
});
