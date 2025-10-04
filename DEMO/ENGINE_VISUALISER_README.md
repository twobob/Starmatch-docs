# Astrological Engine visualiser

This is an interactive HTML visualisation tool for the astrological theme analysis engine (`engine.js`). It combines the ephemeris data system from the DEMO folder with the astrological calculations to create a complete chart analysis tool.

## Files Created

1. **engine-visualiser.html** - Main HTML page with UI controls and display areas
2. **engine-visualiser.css** - Styling that matches the DEMO aesthetic
3. **engine-visualiser.js** - JavaScript integration layer

## Features

### Input Controls
- **Birth Data Input**: Date, time (UTC), latitude, and longitude
- **Engine Settings**:
  - Orb Type (Aspect Orbs vs Planet Orbs)
  - Aspect Orb Sets (5 different configurations)
  - Rulership System (Ancient vs Modern)
  - Precession Correction toggle

### visualisations

1. **Planetary Positions Display**
   - Shows all 10 planets plus Ascendant and Midheaven
   - Displays degrees and zodiac signs
   - Color-coded by sign

2. **Theme Values Chart**
   - Bar chart showing the strength of each of the 12 zodiac themes
   - Normalized values for easy comparison
   - Based on the engine's thematic analysis algorithm

3. **Aspect Analysis**
   - Counts of each aspect type (Conjunction, Opposition, Trine, Square, Sextile, etc.)
   - Derived from engine calculations

4. **Traditional Factors**
   - Polarity distribution (Positive/Negative)
   - Element distribution (Fire, Earth, Air, Water)
   - Quality distribution (Cardinal, Fixed, Mutable)

5. **Dominant Themes**
   - Identifies dominant polarity, element, and quality
   - Shows strongest zodiac theme

6. **Chart Wheel visualisation**
   - Interactive chart wheel showing:
     - 12 zodiac signs (color-coded)
     - Planetary positions
     - Ascendant line
     - Aspect lines between planets

## How It Works

### Ephemeris Integration
The visualiser uses the same ephemeris loading system as `app.js`:
- Loads DE200 ephemeris data from `../data/de200_demo_positions.json`
- Converts birth date/time to Julian Date
- Finds closest ephemeris sample
- Converts 3D cartesian coordinates to zodiacal longitudes

### Engine Integration
The visualiser calls the `getThemeValues()` function from `engine.js` with:
- 10 planetary positions (Sun through Pluto)
- Calculated Ascendant
- Calculated Midheaven

The engine then calculates:
- Theme values for each of the 12 zodiac signs
- Aspect counts
- Traditional factor distributions
- Dominant elements

### Calculations Performed

#### Ascendant & Midheaven
Uses simplified GMST (Greenwich Mean Sidereal Time) calculations:
- Converts Julian Date to sidereal time
- Adjusts for observer longitude
- Calculates Ascendant and Midheaven positions

**Note**: These are simplified calculations for demonstration purposes. Production astrology software would use more precise algorithms.

## Usage

1. Open `engine-visualiser.html` in a web browser
2. Ensure the ephemeris data file exists at `../data/de200_demo_positions.json`
3. Enter birth data (date, time, location)
4. Adjust engine settings if desired
5. Click "Calculate Chart"
6. View the comprehensive analysis results

## Dependencies

- `../engine.js` - The astrological theme analysis engine
- `../data/de200_demo_positions.json` - Ephemeris data
- Modern web browser with Canvas support

## Technical Details

### Coordinate Systems
- **Ecliptic Longitude**: 0-360 degrees measured from the vernal equinox
- **Zodiac Signs**: 12 signs of 30 degrees each
- **Houses**: Calculated from the Ascendant

### Aspect Detection
Uses the engine's built-in aspect detection with configurable orbs:
- Conjunction (0°)
- Opposition (180°)
- Trine (120°)
- Square (90°)
- Sextile (60°)
- Semi-square (45°)
- Semi-sextile (30°)

### Theme Calculation Algorithm
The engine uses a point-based system that considers:
1. Planets in houses and signs
2. Sign rulers and their aspects
3. Exaltations and debilities
4. Dominant elements, qualities, and polarities
5. Mutual receptions
6. Aspect patterns

## Customization

### Colors
Sign colors are defined in CSS using class names like `.sign-aries`, `.sign-taurus`, etc.

### Orb Sets
Five different orb configurations are available in the engine settings.

### Rulership Systems
Switch between ancient (traditional) and modern rulerships for outer planets.

## Known Limitations

1. **Simplified Ascendant/MC**: The house cusp calculations are simplified approximations
2. **Date Range**: Limited to the range of available ephemeris data
3. **Time Zone**: Input time must be in UTC
4. **Aspects to Angles**: Aspects to Ascendant and Midheaven use simplified orbs

## Future Enhancements

Potential improvements:
- More accurate house system calculations (Placidus, Koch, etc.)
- Transit and progression calculations
- Synastry chart comparison
- Additional aspect patterns (Grand Trines, T-Squares, etc.)
- Export chart data
- Print-friendly chart output

## Credits

- **Engine Algorithm**: Will 18
- **Ephemeris Integration**: Based on DEMO/app.js
- **visualisation**: Custom implementation
- **Ephemeris Data**: JPL DE200

## License

Follows the same terms as the original engine.js - "in lieu of copyright" as stated in the engine source.
