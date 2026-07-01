# isoquant_annotator

[IsoQuant](https://github.com/ablab/IsoQuant)-based transcript annotation for protist long-read RNA-seq: sets up a species, runs IsoQuant per platform (PacBio + Nanopore), merges the results, and evaluates the annotation with gene/transcript counts and BUSCO. Built for an HPC/SLURM cluster.

Engine repository — long-read assembly + annotation, *guided* (uses a reference annotation when available, falls back to de novo mode otherwise). Used alongside [LyRic_annotator](https://github.com/Janek21/LyRic_annotator) (unguided) and [geneid-training](https://github.com/Janek21/geneid-training) in a larger protist annotation pipeline.

## Overview

`isoquant_execute.sh` chains the whole thing for one species, via SLURM job dependencies:

1. **Setup** (inline in `isoquant_execute.sh`) — stages the genome FASTA, converts the reference GFF to a clean GTF if one exists (or runs de novo otherwise), selects long-read SRA runs, and submits the ENA read-download array job (skipped if the reads were already downloaded by `LyRic_annotator` on the same machine).
2. **`scripts/isoquant_run.sh`** — SLURM array (task 0 = PacBio, task 1 = Nanopore) that runs IsoQuant for whichever platforms have selected runs, once downloads finish.
3. **`evaluation.sh`** — merges the per-platform `transcript_models.gtf` outputs (AGAT), extracts the longest isoform per gene, predicts ORFs (TransDecoder), counts gene/transcript models, and runs BUSCO (taxon-specific + eukaryote lineages).

`turbo.sh` is a SLURM array job that runs `isoquant_execute.sh` across a whole list of species.

## Repository structure

```
isoquant_annotator/
├── isoquant_execute.sh    # set up a species, submit downloads, launch IsoQuant + evaluation
├── evaluation.sh           # merge per-platform output, infer ORFs, count models, run BUSCO
├── turbo.sh                # SLURM array: runs isoquant_execute.sh over a species list
└── scripts/
    ├── select_accessions.py      # picks the best SRA runs per platform x developmental stage
    ├── srr_dw.sh                  # SLURM array job: downloads SRA reads from ENA
    ├── isoquant_run.sh            # SLURM array job: runs IsoQuant (PacBio / Nanopore)
    ├── get_genetic_code.py        # resolves NCBI nuclear genetic code per taxon
    ├── get_busco_db.py            # resolves the BUSCO lineage per taxon
    ├── buscoPlot.py               # plots BUSCO completeness from result JSONs
    ├── make_summary_tables.py     # aggregates summary/ into final TSV reports
    └── check_runs_status.sh       # reports each species' pipeline stage
```

## Requirements

- SLURM cluster
- conda env `isoquant` (IsoQuant) for the assembly step
- conda env `buscomania`: AGAT, gffread, TransDecoder/TD2, BUSCO, pigz
- NCBI Entrez access (email / optional API key) for taxonomy lookups in `get_genetic_code.py` and `get_busco_db.py`
- `wget` for read downloads
- Reference data one level up:
  `../data/species/<species_name>*/GC*/` (genome FASTA, optionally a reference GFF/GFF3, gzipped or plain) and `../data/longread_protists.tsv`

## Usage

```bash
# run the full pipeline for one species (Genus_species[_extra] naming)
bash isoquant_execute.sh <species_name> [master_tsv] [busco_db]

# or run a batch of species as a SLURM array (one line per species in the list file)
sbatch turbo.sh <species_list.txt>

# check how far each species in a list got
bash scripts/check_runs_status.sh [species_list.txt] [base_dir]

# build the aggregate result tables once species have finished
python3 scripts/make_summary_tables.py [summary_dir]
```

`evaluation.sh` can also be run/sbatch'd directly per species once IsoQuant has produced output.

## Output

- `<species_name>/output/<pacbio|nanopore>/` — raw per-platform IsoQuant output
- `<species_name>/output/eval/` — per-species working files (merged GTF, longest-isoform transcriptome, predicted proteome, BUSCO run dirs)
- `summary/` — central, persistent results, independent of the per-species folders:
  - `counts/` — per-species gene/transcript counts + derived metrics (density, isoforms/gene, coding fraction)
  - `busco_lineage/`, `busco_eukaryote/` — BUSCO JSON results + plots
  - `pred/` — merged predicted annotations
  - `counts_summary.tsv`, `busco_summary.tsv`, `general_summary.tsv` — final aggregate tables from `make_summary_tables.py`
