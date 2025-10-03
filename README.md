# StarMatch

StarMatch is a celestially computational natal chart research tool. The
documentation and demos ship alongside the repository so you can explore the
project completely offline.

## Quick start

1. Clone or download the repository.
2. Open `DEMO/index.html` in any modern browser to interact with the DE200
   barycentric demo. The required ephemeris artifacts (`data/de200.eph` and
   `data/header.200`) as well as the precomputed sample data are already stored
   in the `data/` directory, so no additional downloads are necessary.
3. For a richer exploration of the documentation, open `index.html` from the
   repository root.

If you prefer to serve the files over HTTP instead of `file://`, launch a quick
static server from the repository root:

```bash
python -m http.server 8000
```

Then navigate to <http://localhost:8000/DEMO/> in your browser.
