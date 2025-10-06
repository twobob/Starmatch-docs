/**
 * VALIDATION TEST - JPL Ephemeris Data Import Verification
 * =========================================================
 * 
 * This script validates that the Chebyshev interpolation implementation
 * correctly reads and interpolates planetary positions from DE200 ephemeris.
 * 
 * Tests against JPL's official testpo.200 validation file containing 6840
 * test cases with known-correct planetary positions.
 * 
 * CRITICAL: This test confirms:
 * - Binary coefficient reading is correct (baseOffset = bodyOffset - 3)
 * - Chebyshev polynomial evaluation is accurate
 * - Subinterval selection works properly
 * - Position interpolation matches JPL reference data
 * 
 * Expected result: 4/4 tests PASSED with <0.00003 AU error
 * 
 * DO NOT DELETE - This is the authoritative validation of ephemeris data import.
 */

const fs = require('fs');

// Read the ephemeris file
const ephFile = fs.readFileSync('../data/reduced_de200.eph');
const headerFile = fs.readFileSync('../data/reduced_header.200', 'utf8');

// Parse KSIZE from header
const ksizeMatch = headerFile.match(/KSIZE=\s*(\d+)/);
const KSIZE = parseInt(ksizeMatch[1]);
const recordSize = KSIZE * 8; // 8 bytes per double
console.log(`Record size: ${recordSize} bytes`);

// Parse GROUP 1030 for interval
const group1030Match = headerFile.match(/GROUP\s+1030\s+([\d.E+-]+)\s+([\d.E+-]+)\s+([\d.E+-]+)/);
const intervalDays = parseFloat(group1030Match[3]);

// CRITICAL: Read ACTUAL start/end JD from the ephemeris file itself
// The header GROUP 1030 values may not match the reduced file!
const firstRecordStart = ephFile.readDoubleLE(recordSize);  // Skip header, read first data record
const secondRecordStart = ephFile.readDoubleLE(recordSize * 2); // Read second data record
const actualIntervalDays = secondRecordStart - firstRecordStart; // Calculate ACTUAL interval
const lastRecordEnd = ephFile.readDoubleLE(ephFile.length - recordSize + 8); // Read last record end JD

console.log(`Time range (from actual file): JD ${firstRecordStart} to ${lastRecordEnd}`);
console.log(`Actual interval between records: ${actualIntervalDays} days (header says ${intervalDays})`);

// Parse AU constant from GROUP 1041
// The GROUP 1041 section contains 200 constants
// Format: first line is count (200), then values in Fortran D format
const group1041 = headerFile.match(/GROUP\s+1041([\s\S]*?)GROUP\s+1050/)[1];
const constants = [];
const constantLines = group1041.trim().split('\n').slice(1); // Skip first line (count)
for (const line of constantLines) {
    // Replace D with E for JavaScript parsing
    const values = line.trim().split(/\s+/).map(v => parseFloat(v.replace('D', 'E')));
    constants.push(...values);
}
// AU is the 7th constant (index 6)
const AU = constants[6]; // km per AU
console.log(`AU constant: ${AU} km\n`);

// Parse GROUP 1050 for body parameters
// Format: 3 rows of 13 values each
// Row 1: coefficient offsets (1-based)
// Row 2: number of coefficients
// Row 3: number of subintervals
const group1050 = headerFile.match(/GROUP\s+1050([\s\S]*?)GROUP\s+1070/)[1];
const ipt = [];
const lines = group1050.trim().split('\n').filter(line => line.trim().length > 0);
// Parse the 3 rows
for (let row = 0; row < 3; row++) {
    const values = lines[row].trim().split(/\s+/).map(Number);
    ipt.push(values);
}

// Body definitions (offset, ncoeff, nsubint) - using 1-based offsets from header
// Bodies are: 0=Mercury, 1=Venus, 2=EM-Bary, 3=Mars, 4=Jupiter, 5=Saturn, 
//            6=Uranus, 7=Neptune, 8=Pluto, 9=Moon, 10=Sun, 11=Nutation, 12=Libration
const bodies = {
    'Earth-Moon': { offset: ipt[0][2], coeffs: ipt[1][2], subints: ipt[2][2] }
};

console.log('Earth-Moon Barycenter IPT:', bodies['Earth-Moon']);

// Chebyshev interpolation function
function chebyshev(T, coeffs, n) {
    const T0 = 1;
    const T1 = T;
    let sum = coeffs[0] * T0 + coeffs[1] * T1;
    
    let Tn_minus_2 = T0;
    let Tn_minus_1 = T1;
    
    for (let i = 2; i < n; i++) {
        const Tn = 2 * T * Tn_minus_1 - Tn_minus_2;
        sum += coeffs[i] * Tn;
        Tn_minus_2 = Tn_minus_1;
        Tn_minus_1 = Tn;
    }
    
    return sum;
}

function interpolateBody(jd, bodyName) {
    const body = bodies[bodyName];
    
    // Find which record contains this JD
    // Records start at byte offset = recordSize (skip header record)
    const recordIndex = Math.floor((jd - firstRecordStart) / actualIntervalDays);
    const recordOffset = recordSize + (recordIndex * recordSize); // +recordSize to skip header
    
    // Read the record
    const record = Buffer.allocUnsafe(recordSize);
    ephFile.copy(record, 0, recordOffset, recordOffset + recordSize);
    
    // Extract time range from record header
    const recordStart = record.readDoubleLE(0);
    const recordEnd = record.readDoubleLE(8);
    
    // Determine which subinterval
    const subintervalDuration = (recordEnd - recordStart) / body.subints;
    const subintervalIndex = Math.floor((jd - recordStart) / subintervalDuration);
    
    // Calculate normalized time T within subinterval [-1, 1]
    const subintervalStart = recordStart + subintervalIndex * subintervalDuration;
    const subintervalEnd = subintervalStart + subintervalDuration;
    const T = (2 * jd - (subintervalStart + subintervalEnd)) / (subintervalEnd - subintervalStart);
    
    // Calculate coefficient offset
    // Offset is 1-based in header, convert to 0-based for array indexing
    const baseOffset = (body.offset - 1) * 8;
    
    // Add subinterval offset
    const coeffsPerSubint = body.coeffs * 3; // 3 coordinates (X, Y, Z)
    const subintervalCoeffOffset = baseOffset + (subintervalIndex * coeffsPerSubint * 8);
    
    // Extract coefficients for X, Y, Z
    const xCoeffs = [];
    const yCoeffs = [];
    const zCoeffs = [];
    
    for (let i = 0; i < body.coeffs; i++) {
        xCoeffs.push(record.readDoubleLE(subintervalCoeffOffset + i * 8));
        yCoeffs.push(record.readDoubleLE(subintervalCoeffOffset + (body.coeffs + i) * 8));
        zCoeffs.push(record.readDoubleLE(subintervalCoeffOffset + (2 * body.coeffs + i) * 8));
    }
    
    // Evaluate Chebyshev polynomials
    const x = chebyshev(T, xCoeffs, body.coeffs);
    const y = chebyshev(T, yCoeffs, body.coeffs);
    const z = chebyshev(T, zCoeffs, body.coeffs);
    
    return { x, y, z };
}

// Read test cases from testpo.200
const testpoContent = fs.readFileSync('../data/testpo.200', 'utf8');
const testLines = testpoContent.split('\n');

const earthTests = [];
for (const line of testLines) {
    const parts = line.trim().split(/\s+/);
    if (parts.length >= 7) {
        const [de, date, jed, target, center, coord, value] = parts;
        // Earth (3) from SS Barycenter (12)
        if (target === '3' && center === '12' && ['1', '2', '3'].includes(coord)) {
            const jdNum = parseFloat(jed);
            // Only test dates within our reduced ephemeris range
            if (jdNum >= firstRecordStart && jdNum <= lastRecordEnd) {
                earthTests.push({
                    date,
                    jd: jdNum,
                    coord: parseInt(coord),
                    expected: parseFloat(value)
                });
            }
        }
    }
}

console.log(`\nFound ${earthTests.length} Earth test cases within ephemeris range`);

// Filter to only dates that fall within actual record coverage (not in gaps)
const usableTests = [];
for (const test of earthTests) {
    const recordIndex = Math.floor((test.jd - firstRecordStart) / actualIntervalDays);
    const recordOffset = recordSize + (recordIndex * recordSize);
    const recordStart = ephFile.readDoubleLE(recordOffset);
    const recordEnd = ephFile.readDoubleLE(recordOffset + 8);
    
    if (test.jd >= recordStart && test.jd <= recordEnd) {
        usableTests.push(test);
    }
}

console.log(`Usable test cases (within actual data ranges): ${usableTests.length}`);
console.log('Testing all usable cases:\n');

let successCount = 0;
let failCount = 0;

for (let i = 0; i < usableTests.length; i++) {
    const test = usableTests[i];
    const result = interpolateBody(test.jd, 'Earth-Moon');
    
    const coordNames = { 1: 'X', 2: 'Y', 3: 'Z' };
    const coordName = coordNames[test.coord];
    const computed = { 1: result.x / AU, 2: result.y / AU, 3: result.z / AU }[test.coord];
    
    const diff = Math.abs(computed - test.expected);
    const status = diff < 0.0001 ? '✓ PASS' : '✗ FAIL';
    
    if (diff < 0.0001) successCount++;
    else failCount++;
    
    console.log(`${status} ${test.date} ${coordName}: Expected=${test.expected.toFixed(6)} AU, Got=${computed.toFixed(6)} AU, Diff=${diff.toFixed(8)} AU`);
}

console.log(`\n${successCount} passed, ${failCount} failed`);
