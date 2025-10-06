// Chebyshev polynomial interpolation for DE200 ephemeris
// This replaces N-body integration with actual JPL coefficient interpolation
window.APP_CHEBYSHEV_VERSION = '3.0-heliocentric';

/**
 * Parse ephemeris data from an ArrayBuffer (used when buffer is already loaded)
 * @param {ArrayBuffer} buffer - The ephemeris file data
 * @returns {Object} - Ephemeris records and metadata
 */
function parseEphemerisBuffer(buffer) {
  const view = new DataView(buffer);
  
  const RECORD_SIZE = 1652 * 8; // KSIZE * sizeof(double) = 13216 bytes
  const records = [];
  
  // Skip header record (first 13216 bytes)
  let offset = RECORD_SIZE;
  
  while (offset + RECORD_SIZE <= buffer.byteLength) {
    // Read start and end JD (little-endian doubles)
    const startJD = view.getFloat64(offset, true);
    const endJD = view.getFloat64(offset + 8, true);
    
    // Validate JD values
    if (startJD < 2000000 || startJD > 3000000 || endJD < 2000000 || endJD > 3000000) {
      break; // Invalid data, stop reading
    }
    
    // Read all coefficients (1650 doubles after the 2 JD values)
    const coefficients = new Float64Array(1650);
    for (let i = 0; i < 1650; i++) {
      coefficients[i] = view.getFloat64(offset + 16 + (i * 8), true);
    }
    
    records.push({
      startJD,
      endJD,
      coefficients
    });
    
    offset += RECORD_SIZE;
  }
  
  console.log(`Parsed ${records.length} ephemeris records from buffer`);
  if (records.length > 0) {
    console.log(`Date range: JD ${records[0].startJD} to ${records[records.length-1].endJD}`);
  }
  
  return {
    records,
    recordSize: RECORD_SIZE,
    startJD: records[0].startJD,
    endJD: records[records.length - 1].endJD
  };
}

/**
 * Load and parse the binary ephemeris file
 * @param {string} ephUrl - URL to the .eph file
 * @param {Object} header - Parsed header with KSIZE info
 * @returns {Object} - Ephemeris records and metadata
 */
async function loadEphemerisData(ephUrl, header) {
  const response = await fetch(ephUrl);
  if (!response.ok) {
    throw new Error(`Failed to load ephemeris: ${response.status}`);
  }
  
  const buffer = await response.arrayBuffer();
  return parseEphemerisBuffer(buffer);
}

// Export for use in other modules
window.parseEphemerisBuffer = parseEphemerisBuffer;
window.loadEphemerisData = loadEphemerisData;

/**
 * Chebyshev polynomial evaluation
 * @param {number} x - Normalized time in [-1, 1]
 * @param {Float64Array} coeffs - Chebyshev coefficients
 * @returns {number} - Interpolated value
 */
function evaluateChebyshev(x, coeffs) {
  const n = coeffs.length;
  if (n === 0) return 0;
  if (n === 1) return coeffs[0];
  
  // Clenshaw algorithm for stable evaluation
  let b_k2 = 0; // b_{k+2}
  let b_k1 = 0; // b_{k+1}
  
  for (let k = n - 1; k >= 1; k--) {
    const b_k = 2 * x * b_k1 - b_k2 + coeffs[k];
    b_k2 = b_k1;
    b_k1 = b_k;
  }
  
  return x * b_k1 - b_k2 + coeffs[0];
}

/**
 * Find the record containing a specific JD
 * @param {Array} records - Array of ephemeris records
 * @param {number} jd - Julian date to find
 * @returns {Object|null} - Record or null if not found
 */
function findRecordForJD(records, jd) {
  // Binary search for efficiency
  let left = 0;
  let right = records.length - 1;
  
  while (left <= right) {
    const mid = Math.floor((left + right) / 2);
    const record = records[mid];
    
    if (jd < record.startJD) {
      right = mid - 1;
    } else if (jd > record.endJD) {
      left = mid + 1;
    } else {
      return record; // Found it!
    }
  }
  
  return null; // JD not in any record
}

/**
 * Interpolate position and velocity for a body at a specific time
 * @param {Object} record - Ephemeris record
 * @param {number} jd - Julian date
 * @param {number} bodyOffset - Offset in coefficients for this body (1-indexed per JPL convention)
 * @param {number} numCoeffs - Number of coefficients STORED per component (for layout)
 * @param {number} subintervals - Number of subintervals within this record
 * @param {number} usableCoeffs - Number of coefficients to ACTUALLY USE (defaults to numCoeffs)
 * @returns {Object} - {position: [x,y,z], velocity: [vx,vy,vz]} in AU and AU/day
 */
function interpolateBody(record, jd, bodyOffset, numCoeffs, subintervals, usableCoeffs = numCoeffs) {
  // The offset is 1-indexed position in the RECORD (including 2 JD values at start)
  // Coefficient array is 0-indexed and does NOT include the JDs
  // So: position 3 in record = index 0 in coefficients array
  // Therefore: array index = bodyOffset - 3
  const baseOffset = bodyOffset - 3;
  
  // Determine which subinterval this JD falls into
  const recordSpan = record.endJD - record.startJD;
  const subintervalSpan = recordSpan / subintervals;
  const normalizedTime = (jd - record.startJD) / recordSpan;
  const subintervalIndex = Math.min(Math.floor(normalizedTime * subintervals), subintervals - 1);
  
  // Normalize time within this specific subinterval to [-1, 1]
  const subintervalStart = record.startJD + (subintervalIndex * subintervalSpan);
  const subintervalEnd = subintervalStart + subintervalSpan;
  const t = ((jd - subintervalStart) / subintervalSpan) * 2 - 1;
  
  const position = [0, 0, 0];
  const velocity = [0, 0, 0];
  
  // Each subinterval has 3 components (x, y, z), each with numCoeffs coefficients
  // Layout PER SUBINTERVAL: [X coeffs, Y coeffs, Z coeffs]
  const coeffsPerSubinterval = numCoeffs * 3;
  
  for (let component = 0; component < 3; component++) {
    // Offset for this subinterval, then offset for this component within the subinterval
    const subintervalCoeffOffset = baseOffset + (subintervalIndex * coeffsPerSubinterval);
    const componentOffset = subintervalCoeffOffset + (component * numCoeffs);
    
    // Read stored coefficients but only use usableCoeffs for Chebyshev (handles DE200 padding)
    const allCoeffs = record.coefficients.slice(componentOffset, componentOffset + numCoeffs);
    const coeffs = allCoeffs.slice(0, usableCoeffs);
    
    // Evaluate Chebyshev polynomial for position
    position[component] = evaluateChebyshev(t, coeffs);
    
    // For velocity, compute derivative of Chebyshev polynomial
    // Using Chebyshev derivative recurrence: T'_n(x) = n * U_{n-1}(x)
    // For simplicity, use finite differences (more accurate derivative would use proper recurrence)
    const dt = 0.0001; // Small normalized time step
    const t_plus = Math.min(1.0, t + dt);
    const t_minus = Math.max(-1.0, t - dt);
    const val_plus = evaluateChebyshev(t_plus, coeffs);
    const val_minus = evaluateChebyshev(t_minus, coeffs);
    
    // dPosition/dt in normalized subinterval time, convert to days
    const dpos_dnorm = (val_plus - val_minus) / (2 * dt);
    velocity[component] = dpos_dnorm * (2 / subintervalSpan); // Convert from normalized to JD units (AU/day)
  }
  
  return { position, velocity };
}

/**
 * Generate samples from ephemeris coefficients
 * @param {Object} ephemeris - Loaded ephemeris data
 * @param {Object} constants - Header constants (for AU, body names, etc.)
 * @param {Object} options - Options like outputStride
 * @returns {Object} - Dataset compatible with existing visualization
 */
function generateSamplesFromEphemeris(ephemeris, constants, options = {}) {
  const auKm = constants.AU;
  
  // DE200 IPT (Index Pointer Table) from GROUP 1050 in header.200
  // Format: [offset_start, num_coeffs, num_subintervals] for each body
  // Bodies: Mercury, Venus, Earth-Moon Bary, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto, Moon, Sun, Nutations
  // Offsets:  3, 147, 183, 273, 303, 330, 354, 378, 396, 414, 702, 747
  // Num coeffs per component: 12, 12, 15, 10, 9, 8, 8, 6, 6, 12, 15, 10
  // Subintervals: 4, 1, 2, 1, 1, 1, 1, 1, 1, 8, 1, 4
  
  const bodies = [
    { name: 'Mercury', offset: 3, coeffs: 12, subintervals: 4 },
    { name: 'Venus', offset: 147, coeffs: 12, subintervals: 1 },
    { name: 'Earth-Moon Barycenter', offset: 183, coeffs: 15, subintervals: 2, usableCoeffs: 13 },  // Layout has 15, but only first 13 are Chebyshev coeffs  // Restored to 15 per header.200
    { name: 'Mars', offset: 273, coeffs: 10, subintervals: 1 },
    { name: 'Jupiter', offset: 303, coeffs: 9, subintervals: 1 },
    { name: 'Saturn', offset: 330, coeffs: 8, subintervals: 1 },
    { name: 'Uranus', offset: 354, coeffs: 8, subintervals: 1 },
    { name: 'Neptune', offset: 378, coeffs: 6, subintervals: 1 },
    { name: 'Pluto', offset: 396, coeffs: 6, subintervals: 1 },
    { name: 'Moon', offset: 414, coeffs: 12, subintervals: 8 },
    { name: 'Sun', offset: 702, coeffs: 15, subintervals: 1, usableCoeffs: 13 }  // Layout has 15, but only first 13 are Chebyshev coeffs  // Restored to 15 per header.200
    // Nutations at offset 747 with 10 coeffs, 4 subintervals - not needed for planetary positions
  ];
  
  const samples = [];
  
  // Generate data for the full ephemeris range
  const startJD = ephemeris.startJD;
  const endJD = ephemeris.endJD;
  const outputStride = 5; // Sample every 5 days for smooth orbits
  
  let sampleCount = 0;
  let skippedCount = 0;
  
  for (let jd = startJD; jd <= endJD; jd += outputStride) {
    const record = findRecordForJD(ephemeris.records, jd);
    if (!record) {
      // Skip dates outside ephemeris range silently
      skippedCount++;
      continue;
    }
    
    sampleCount++;
    const positions_km = {};
    
    // Get all positions from ephemeris - RAW barycentric (NO CONVERSION)
    for (const body of bodies) {
      const { position } = interpolateBody(
        record, 
        jd, 
        body.offset, 
        body.coeffs, 
        body.subintervals,
        body.usableCoeffs  // Use usableCoeffs if specified (for DE200 padding)
      );
      positions_km[body.name] = position;
    }
    
    // CONVERT TO HELIOCENTRIC: Subtract Sun's position from all bodies
    const sunPos = positions_km['Sun'];
    const heliocentric_km = {};
    
    for (const [bodyName, pos] of Object.entries(positions_km)) {
      heliocentric_km[bodyName] = [
        pos[0] - sunPos[0],
        pos[1] - sunPos[1],
        pos[2] - sunPos[2]
      ];
    }
    
    samples.push({
      julian_date: jd,
      positions_km: heliocentric_km
    });
  }
  
  return {
    metadata: {
      description: 'Chebyshev interpolation from DE200 ephemeris coefficients (HELIOCENTRIC)',
      start_julian_date: startJD,
      end_julian_date: endJD,
      output_stride_days: outputStride,
      au_km: auKm,
      source: 'DE200 ephemeris file',
      coordinate_system: 'heliocentric',
      version: '3.0-heliocentric'
    },
    bodies: bodies.map(b => b.name),
    samples
  };
}

// Export for use in main app
if (typeof module !== 'undefined' && module.exports) {
  module.exports = {
    loadEphemerisData,
    generateSamplesFromEphemeris,
    evaluateChebyshev,
    findRecordForJD,
    interpolateBody
  };
}
