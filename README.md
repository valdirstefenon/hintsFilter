# hintsFilter GFF

**Annotation filter for BRAKER — transcript-level confidence for cleaner gene sets.**

[![Python](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey.svg)](#installation)

---

## Overview

**hintsFilter GFF** is a Python tool that filters gene annotations produced by **BRAKER** (or any GFF/GTF file) using external evidence (**hints**) and quality criteria. It produces a **confident transcript set** and a **confident gene set**, plus filtered functional annotation files from **EggNOG** and **InterProScan**.

The tool ships with a **Tkinter GUI** so biologists and bioinformaticians can run the pipeline without command-line expertise.

Version **3.3** introduces **transcript-level confidence**, separating evidence from EggNOG, InterProScan, and orphan (no-hit) transcripts, and reporting both transcript-level and gene-level summaries.

---

## Key Features

- **Transcript-level filtering** — confidence is decided per transcript, not per gene.
- **Evidence-aware classification** — separates EggNOG hits, InterProScan hits, and orphans.
- **Quality filters** — protein length, exon count, and BRAKER score (ab initio mode).
- **Hints-aware** — checks hints overlapping each transcript within a user-defined margin.
- **Orphan detection** — flags transcripts that pass quality but lack functional annotation (candidates for lineage-specific genes; require validation).
- **Two-level report** — transcript-level and gene-level summaries that add up correctly.
- **Cross-platform GUI** with soft green tones — Linux, macOS, Windows.
- **UTF-8 / latin-1 safe** file handling.

---

## What's New in 3.3

| Change | Description |
|---|---|
| **Transcript-level confidence** | Confidence is now decided for each mRNA, then propagated to genes. |
| **New outputs** | `eggnog_confident.tsv`, `interpro_confident.tsv`, `orphans_confident.tsv`, `transcripts_confident.txt`, `braker_confident.gff3`. |
| **Two-level report** | Counts are reported at the transcript level and at the gene level, without double counting. |
| **Orphan terminology** | "Orphans" replaced the misleading "species-specific" label. Orphans are transcripts with no EggNOG/InterPro hit that pass quality filters. |
| **Cleaner GFF3 output** | Only confident transcripts (and their sub-features) are written, not entire genes. |
| **Updated defaults** | More lenient defaults suited to BRAKER output (`min_protein_len = 80`, `min_exons = 1`, `min_score = 1.0`). |

---

## How It Works

1. **Load hints** — reads the hints GFF, indexes by chromosome, sorts by position.
2. **Parse GFF** — extracts genes and transcripts. For each transcript, computes:
   - exon count,
   - protein length (sum of CDS lengths ÷ 3),
   - BRAKER score (if present in the score column).
3. **Apply hints to transcripts** — a transcript is *hinted* if at least one hint with score ≥ `min_score` overlaps its coordinates expanded by `margin` bp.
4. **Apply quality filters** — protein length ≥ `min_protein_len`, exon count ≥ `min_exons`, and (only in ab initio mode) BRAKER score ≥ `min_braker_score`.
5. **Classify transcripts** — a transcript is **confident** if:

   ```
   confident = passes_quality  AND  (has_hint OR ab_initio_mode)
   ```

   Confident transcripts are then split by evidence:
   - **EggNOG confident** — confident and has an EggNOG hit.
   - **InterProScan confident** — confident and has an InterProScan hit.
   - **Orphan confident** — confident and has **no** EggNOG and **no** InterProScan hit.

6. **Propagate to genes** — a gene is confident if at least one of its transcripts is confident.
7. **Write outputs** — filtered GFF3, filtered EggNOG/InterProScan files, orphan list, transcript and gene lists, and a two-level report.

---

## Understanding the Logic

### Default mode (with hints)

A transcript is **confident** if:

- it passes all quality filters (`min_protein_len`, `min_exons`), **and**
- at least one hint (score ≥ `min_score`) overlaps its coordinates within `margin` bp.

The BRAKER score is **ignored** in this mode.

### Ab initio mode

Hints are ignored. A transcript is **confident** if:

- protein length ≥ `min_protein_len`, **and**
- exon count ≥ `min_exons`, **and**
- BRAKER score ≥ `min_braker_score`.

Use this when you have no external evidence and want to rely on BRAKER's intrinsic confidence.

### Orphans

Orphans are confident transcripts with **no functional annotation** in EggNOG or InterProScan. They are **candidates** for lineage-specific genes, but they may also reflect:

- fragmented BRAKER models,
- non-coding transcripts or long UTRs mis-annotated as CDS,
- fast-evolving genes with no detectable similarity,
- gaps or limitations in EggNOG/InterProScan databases,
- isoform-specific effects (one isoform annotated, another not).

Treat them as **hypotheses**, not as validated species-specific genes.

---

## Installation

### Prerequisites

- **Python 3.8+** (tested on 3.8, 3.9, 3.10, 3.11, 3.12)
- **Tkinter** — usually bundled with Python. On Debian/Ubuntu:
  ```bash
  sudo apt install python3-tk
  ```
- **Python libraries**: `pandas`, `numpy`

### Install dependencies

```bash
pip install pandas numpy
```

Or with a virtual environment:

```bash
python3 -m venv venv
source venv/bin/activate      # Linux / macOS
# venv\Scripts\activate       # Windows
pip install pandas numpy
```

### Get the script

```bash
git clone https://github.com/valdirstefenon/hintsFilter.git
cd hintsFilter
```

---

## Usage

### GUI mode (recommended)

```bash
python3 hintsFilter.py
```

1. Use **Browse...** to select:
   - BRAKER GFF3/GTF file
   - Hints file (optional in ab initio mode)
   - EggNOG `.emapper.annotations` file
   - InterProScan `.tsv` file
   - Output directory
2. Adjust parameters in the **Filtering Parameters** panel.
3. Click **RUN PIPELINE**.
4. Follow progress in the **Execution Log**.

---

## Parameters

| Parameter | GUI label | Description | Default |
|---|---|---|---|
| `min_score` | Min. hints score | Minimum hint score to be considered. | `1.0` |
| `margin` | Margin (bp) | Flanking region around each transcript to search for hints. | `500` |
| `min_protein_len` | Min. protein (aa) | Minimum protein length to pass quality. | `80` |
| `min_exons` | Min. exons | Minimum number of exons to pass quality. | `1` |
| `min_braker_score` | BRAKER score ≥ | Minimum BRAKER score (only used in ab initio mode). | `0.0` |
| `use_abinitio` | Ab initio mode | Ignore hints; use only quality filters. | `False` |
| `run_analysis` | Run GFF analysis | Run the integrated GFF statistics module. | `True` |

---

## Input Files

| File | Description | Format |
|---|---|---|
| **GFF/GTF** | Gene annotation from BRAKER (or compatible tool). | GFF3 / GTF |
| **Hints** | Evidence file (e.g., `hintsfile.gff`) from BRAKER. | GFF |
| **EggNOG** | Functional annotation from EggNOG-mapper. | `*.emapper.annotations` |
| **InterProScan** | Protein domain annotation from InterProScan. | TSV |

> In ab initio mode, the hints file is optional and the BRAKER score filter becomes active.

---

## Output Files

All outputs are written to the user-specified directory.

| File | Description |
|---|---|
| `braker_confident.gff3` | Filtered GFF3 with only confident genes and confident transcripts (plus their exons, CDS, etc.). |
| `transcripts_confident.txt` | List of confident transcript IDs. |
| `genes_confident.txt` | List of genes with at least one confident transcript. |
| `transcripts_without_hints.txt` | Transcripts with no overlapping hints (diagnostic). |
| `eggnog_confident.tsv` | EggNOG hits restricted to confident transcripts. |
| `interpro_confident.tsv` | InterProScan hits restricted to confident transcripts. |
| `orphans_confident.tsv` | Confident transcripts with no functional annotation (with coordinates, protein length, exon count, score, hint count). |
| `pipeline_report.txt` | Two-level summary (transcript and gene) with parameters and file paths. |
| `filtered_gff_analysis_report.txt` | Detailed statistics of the filtered GFF (if GFF analysis enabled). |
| `filtered_gene_details.csv` | Per-gene table with exon/CDS/transcript counts, start/stop codon presence, length. |

---

## Tips for Best Results

1. **Inspect the GFF first.** Make sure it contains `exon` and `CDS` features. If BRAKER was run with `--prot_seq` or `--rnaseq`, the score column may be empty — that is why the BRAKER score filter is disabled by default.
2. **Use the log statistics.** The pipeline prints min/max/mean of protein length, exon count, and BRAKER score. Use them to choose realistic cutoffs.
3. **Margin tuning.** For protein-based hints, 300–500 bp is safe. For RNA-seq hints, 100–200 bp is usually enough.
4. **Start lenient, then tighten.** Run once with permissive parameters, inspect the report, then tighten thresholds to see how many transcripts are retained.
5. **Orphans are candidates, not conclusions.** Validate orphan transcripts with additional evidence (expression, synteny, BLASTp against nr, orthology) before calling them lineage-specific.

---

## Troubleshooting

| Issue | Likely cause | Solution |
|---|---|---|
| No transcripts kept with strict filters | Protein length or exon cutoffs too high. | Check log statistics; lower thresholds. |
| BRAKER score always 0 | Score column contains `.` (common with `--prot_seq`). | Keep the BRAKER score filter disabled (default). |
| GUI does not open | Tkinter not installed. | Install Tkinter (`sudo apt install python3-tk` on Ubuntu). |
| EggNOG/InterProScan files empty | ID matching failed. | Check that protein IDs in EggNOG/InterProScan match transcript IDs in the GFF (e.g., `g1.t1`). Adjust the regex in `match_transcript_id()` if needed. |
| Too many orphans | Databases lack representation for your lineage, or models are fragmented. | Inspect orphan coordinates and protein lengths; consider stricter quality filters. |

---

## Project Structure

```
hintsFilter/
├── hintsFilter.py        # Main script (GUI + pipeline)
├── README.md             # This file
└── LICENSE               # MIT license
```

---

## License

Distributed under the **MIT License**. See `LICENSE` for details.

---

## Citation

If you use this tool in your research, please cite:

> **hintsFilter GFF — Annotation filter for BRAKER.**
> GitHub: https://github.com/valdirstefenon/hintsFilter

---

## valdir.stefenon@ufsc.br

For questions, bug reports, or feature requests, please open an issue on the GitHub repository.

---

*Enjoy filtering your BRAKER annotations!*
