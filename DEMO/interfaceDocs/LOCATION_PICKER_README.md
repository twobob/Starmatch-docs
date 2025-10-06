# Location Picker Feature

## Overview
The location picker allows users to search for any location worldwide and automatically populate the latitude and longitude fields needed for accurate astrological chart calculations.

## Features

### 🔍 Location Search
- **Search by city name**: "London", "New York", "Tokyo"
- **Search with country**: "Paris, France", "Sydney, Australia"
- **Search hospitals**: "New York Hospital", "St Mary's Hospital London"
- **Search addresses**: Full street addresses for precise locations

### 🌍 Geocoding Service
Uses **Nominatim** (OpenStreetMap's free geocoding API):
- ✅ No API key required
- ✅ Worldwide coverage
- ✅ Major cities and hospitals included
- ✅ Returns precise latitude/longitude coordinates
- ✅ Up to 10 results per search

### 📍 Usage

1. **Click the 🔍 button** next to the Latitude/Longitude fields
2. **Enter a location** in the search box
3. **Click Search** or press Enter
4. **Select a location** from the results
5. **Click "Use Selected Location"** to apply

The latitude and longitude fields will be automatically populated with the selected location's coordinates.

## Technical Details

### Components

#### `location-picker.js`
- `LocationPicker` class handles all search and UI logic
- Communicates with Nominatim API
- Manages modal state and user interactions
- Formats and displays search results

#### `location-picker.css`
- Modern dark theme matching the visualiser
- Responsive modal design
- Smooth animations and transitions
- Accessible UI elements

### API Usage

**Endpoint**: `https://nominatim.openstreetmap.org/search`

**Query Parameters**:
- `q`: Search query
- `format`: json
- `addressdetails`: 1 (include full address)
- `limit`: 10 (max results)
- `accept-language`: en

**Rate Limiting**: 
- Nominatim has a usage policy of max 1 request per second
- Please be respectful of the free service

### Response Format

```javascript
{
  name: "London",
  fullAddress: "London, Greater London, England, United Kingdom",
  latitude: 51.5074,
  longitude: -0.1278,
  type: "city",
  importance: 0.8
}
```

## Examples

### Major Cities
- "New York, USA" → 40.7128°, -74.0060°
- "Tokyo, Japan" → 35.6762°, 139.6503°
- "London, UK" → 51.5074°, -0.1278°
- "Paris, France" → 48.8566°, 2.3522°
- "Sydney, Australia" → -33.8688°, 151.2093°

### Hospitals (Common Birth Locations)
- "St Thomas' Hospital London"
- "Mount Sinai Hospital New York"
- "Toronto General Hospital"
- "Royal Women's Hospital Melbourne"
- "Lenox Hill Hospital"

### Specific Addresses
- "1600 Pennsylvania Avenue, Washington DC"
- "10 Downing Street, London"
- "Buckingham Palace"

## Privacy & Data

- **No data is stored** - searches are performed client-side
- **No personal information collected**
- Uses public OpenStreetMap data
- Coordinates are only used locally for chart calculations

## Future Enhancements

Potential additions:
- 🗺️ Interactive map interface (Leaflet.js)
- 📌 Click-to-select on map
- ⭐ Save favorite locations
- 🕰️ Timezone detection from coordinates
- 📱 Geolocation API (use current location)
- 🏥 Hospital database for common birth locations

## Troubleshooting

### "No locations found"
- Try a more general search (e.g., "London" instead of "123 Random Street")
- Include country name for disambiguation
- Check spelling

### Search not working
- Check internet connection (requires online access to Nominatim)
- Check browser console for errors
- Ensure you're not exceeding rate limits (1 req/sec)

### Wrong location selected
- Read the full address carefully before selecting
- Some cities have same names in different countries
- Use country name to disambiguate

## Attribution

Location data © [OpenStreetMap](https://www.openstreetmap.org/copyright) contributors

Geocoding service provided by [Nominatim](https://nominatim.openstreetmap.org/)
