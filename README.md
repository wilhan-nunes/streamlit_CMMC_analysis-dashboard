# CMMC Analysis Dashboard

Interactive Streamlit dashboard for exploring **CMMC enrichment** results together with **GNPS2 Feature-Based Molecular Networking (FBMN)** outputs and sample metadata.

The app lets you:

- combine CMMC enrichment annotations with FBMN quantification tables
- merge those results with sample metadata
- compare metabolite abundances across groups with statistical boxplots
- inspect source/origin overlaps with UpSet plots
- explore molecular network neighborhoods
- launch MicrobeMASST searches from enriched features

## What You Need

For a full analysis with your own data you typically need:

- a **CMMC Enrichment Task ID**
- an **FBMN Task ID**
- a metadata table with a required `filename` column

The app can also run with built-in demo data, which is the fastest way to confirm everything is working.

## Metadata Format

Your metadata file must contain:

- `filename`: sample file name matching the FBMN quantification table
- one or more grouping columns, ideally prefixed with `ATTRIBUTE_`

Accepted file types: `.csv`, `.tsv`, `.txt`, `.xlsx`

| filename | ATTRIBUTE_group | ATTRIBUTE_sex | ATTRIBUTE_batch |
|---|---|---|---|
| Sample_01.raw | case | F | 1 |
| Sample_02.raw | control | M | 1 |
| Sample_03.raw | case | F | 2 |

## Local Setup

**Python 3.10 is required.**

```bash
git clone https://github.com/wilhan-nunes/streamlit_CMMC_analysis-dashboard
cd streamlit_CMMC_analysis-dashboard
git submodule update --init --recursive

python3.10 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

Start the app:

```bash
streamlit run app.py
```

By default Streamlit serves on `http://localhost:8501`. To use a different port:

```bash
streamlit run app.py --server.port 8512
```

Task IDs can also be passed via URL query parameters:

```
http://localhost:8501/?cmmc_task_id=<id>&fbmn_task_id=<id>
```

## Docker Setup

```bash
docker compose -f docker-compose.yml -f docker-compose-dev.yml up --build
```

Then open `http://localhost:5000`. Helper Makefile targets are also available:

```bash
make server-compose-interactive   # build + run in foreground
make server-compose               # build + run in background
make server-compose-production    # production stack (no dev overlay)
```

## Usage

1. Launch the app.
2. In the sidebar, either enable **Load Example Data** or enter your own CMMC and FBMN task IDs.
3. Provide metadata — use the table bundled with the FBMN task or upload your own file.
4. Optionally enable **Use uploaded quantification table** and upload a custom quant file.
5. Click **Run Analysis**.

Results are shown in two tabs:

- **Data Explorer** — statistical boxplots and UpSet plots
- **Advanced Visualizations** — molecular network neighborhoods and MicrobeMASST search

### Tips

- If **Include all features in the analysis** is enabled the merged table includes features without CMMC matches, which can significantly increase processing time.
- Boxplots will be empty if metadata `filename` values do not match the FBMN sample names after stripping extensions (`.mzML`, `.mzXML`, `.raw`).
- The fastest first run is the built-in demo dataset ([Quinn et al. 2020](https://doi.org/10.1038/s41586-020-2047-9)): germ-free vs. colonized mice across multiple body sites).

## Repository Layout

- [`app.py`](app.py) — main Streamlit entry point; sidebar, data processing pipeline, tab rendering
- [`utils.py`](utils.py) — GNPS2 API fetching, data preparation, SMILES/RDKit helpers
- [`enhanced_boxplot.py`](enhanced_boxplot.py) — statistical boxplot UI with group selection and significance testing
- [`network_cluster_plotter.py`](network_cluster_plotter.py) — Plotly molecular network visualization
- [`upset_plot.py`](upset_plot.py) — UpSet plot generation by source or origin
- [`microbemass_frame.py`](microbemass_frame.py) — MicrobeMASST search interface
- [`box_plot.py`](box_plot.py) — lower-level boxplot utilities
- [`welcome.py`](welcome.py) — landing page shown before analysis runs
- [`microbe_masst/`](microbe_masst/) — git submodule for MicrobeMASST search backend

## Troubleshooting

- **`Port 8501 is already in use`** — run on another port: `streamlit run app.py --server.port 8512`
- **`No metadata table was found for this FBMN task`** — upload a metadata file manually and confirm it includes a `filename` column
- **`Task ID must be exactly 32 characters long`** — double-check that you copied the full GNPS2 task identifier
- **Boxplots are empty or groups are missing** — confirm that metadata filenames match the FBMN sample names after extensions are stripped

## References

- [CMMC enrichment workflow](https://cmmc.gnps2.org/network_enrichment/)
- [FBMN documentation](https://wang-bioinformatics-lab.github.io/GNPS2_Documentation/fbmn/)
- [CMMC knowledge base](https://cmmc-kb.gnps2.org)
