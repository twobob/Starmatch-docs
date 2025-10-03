# DE200 Barycentric Demo

This folder contains a self-contained browser demo that visualizes barycentric
positions for the major solar-system bodies. The display loads precomputed
samples generated from the official JPL DE200 ephemeris constants
(`data/header.200`).

A lightweight Newtonian integrator seeded with the DE200 state vectors produces
the positions stored in `data/de200_demo_positions.json`. The integration is
intended for qualitative exploration and should not be used for navigation or
precision astronomy.

The repository already includes the canonical `de200.eph` and `header.200`
artifacts inside the top-level `data/` directory. You can therefore open
`index.html` directly in a browser without running any helper scripts.

If you need to refresh the artifacts for verification purposes, run:

```bash
python DEMO/fetch_ephemeris.py
```

This downloads verified copies of the ephemeris files into `data/`.
