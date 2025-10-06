# Ephemeris Data Validation

## Validation Test

**File:** `test_against_testpo.js`

This is the **authoritative validation** that confirms our ephemeris data import is correct.

### What it does:
- Tests Chebyshev interpolation against JPL's official test data (testpo.200)
- Validates 4 test cases for Earth positions across different time periods
- Confirms coefficient reading, polynomial evaluation, and interpolation accuracy

### How to run:
```bash
cd DEMO
node test_against_testpo.js
```

### Expected output:
```
✓ Test 1 PASSED: Earth 1904-01-01 Z = -0.265604 AU (error: 0.00000785 AU)
✓ Test 2 PASSED: Earth 1925-01-01 Z = 0.386123 AU (error: 0.00000193 AU)
✓ Test 3 PASSED: Earth 1952-01-01 Y = -0.611072 AU (error: 0.00002584 AU)
✓ Test 4 PASSED: Earth 1996-01-01 X = 0.766646 AU (error: 0.00000748 AU)

4 passed, 0 failed
```

**All errors < 0.00003 AU confirms the implementation is correct!** ✅

## Visual Integration Tests

### `index.html` - Orbital Visualization Test
**Purpose:** Visual validation that planetary positions form circular orbits at correct distances.

**What it tests:**
- Binary ephemeris file loading and parsing
- Chebyshev interpolation across time range
- Heliocentric coordinate conversion (Sun at origin)
- Orbital mechanics accuracy

**How to test:**
1. Open `http://localhost:8000/DEMO/index.html` in browser
2. Check browser console for data quality output
3. Verify the visualization shows:
   - ✅ Sun at center (0,0)
   - ✅ Circular orbits (not diagonal lines or ellipses)
   - ✅ Correct distances: Mercury ~0.4 AU, Earth ~1.0 AU, Mars ~1.5 AU, Jupiter ~5.2 AU

**Console output should show:**
```
  DATA QUALITY CHECK - First Sample:
  Sun: distance from origin = 0.000 AU
  Mercury: distance from origin = 0.464 AU
  Venus: distance from origin = 0.720 AU
  Earth-Moon Barycenter: distance from origin = 0.984 AU
  Mars: distance from origin = 1.621 AU
  Jupiter: distance from origin = 5.338 AU
```

**File selection:**
- Default: Uses `reduced_de200.eph` (1900-2035, 9.7 MB)
- Add `?reduced=0` to URL to use full `de200.eph` (1599-2169, 40.9 MB)

### `engine.html` - Astrological Chart Calculation Test
**Purpose:** End-to-end validation that planetary positions are correctly converted to chart data.

**What it tests:**
- Date/time to Julian Date conversion
- Position interpolation for specific moments
- Ecliptic longitude calculation from Cartesian coordinates
- **Critical:** Sun position from Earth's perspective (inverted heliocentric)
- **Critical:** Moon position from Moon ephemeris 

**How to test:**
1. Open `http://localhost:8000/DEMO/engine.html` in browser
2. Check console shows data quality validation on load
3. Calculate charts for different birth dates and verify:
   - ✅ Sun changes ~1° per day (moves through zodiac signs by birth date)
   - ✅ Moon changes ~13° per day (different even for same-day births at different times)
   - ✅ Planetary positions vary realistically


**Starmatch mode:** Tests comparison calculations use corrected Moon/Sun positions for both subject and target charts.

## Integration Tests

The live application (`engine.html`) also validates data quality on load by displaying:
- Planetary distances from the Sun at two sample dates
- Expected values: Sun ~0 AU, Earth ~1.0 AU, Jupiter ~5.2 AU, etc.

This confirms the data pipeline from binary ephemeris → Chebyshev interpolation → chart calculations is working correctly.
