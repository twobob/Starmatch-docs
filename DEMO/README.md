# DE200 Barycentric Demo

This folder contains a self-contained browser demo that visualizes barycentric
positions for the major solar-system bodies. The display loads precomputed
samples generated from the official JPL DE200 ephemeris constants
(`data/header.200`).

A lightweight Newtonian integrator seeded with the DE200 state vectors produces
the positions stored in `data/de200_demo_positions.json` (mirrored to
`DEMO/de200_demo_positions.json` so the demo can be hosted as a standalone
folder). When the demo runs over HTTP/S it will fetch the canonical JSON data
and fall back to the legacy `data/de200_demo_positions.js` bundle if present.
The integration is
intended for qualitative exploration and should not be used for navigation or
precision astronomy.

Because the project infrastructure rejects binary uploads, the canonical
`de200.eph` file is not committed to the repository. Run:

```bash
python DEMO/fetch_ephemeris.py
```

to download verified copies of `de200.eph` and `header.200` directly into the
`data/` directory. Once fetched you can open `index.html` directly in a browser
and explore the solar system snapshot in an entirely client-side experience.
Running `python DEMO/generate_positions.py` will regenerate the dataset and
refresh the mirrored copy used by the demo.
