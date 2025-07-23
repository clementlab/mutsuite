# mutsuite

A framework for simulating, calling, and benchmarking indel and structural variant detection in whole genome sequencing (WGS) data.

## Features
- mutSim: Simulate indels and other mutations in BAM files
- mutRun: Orchestrate the running of many mutations with varying mutation parameters
- mutTest: Run multiple variant callers (Pindel, Mutect, Lofreq, Strelka, SomaticSniper, Vardict, Varscan, and more)
- mutAgg: Aggregate and compare results from different callers

## Directory Structure
- `docs/` — Example analysis scripts and notebooks
- `src/` — MutSim source code

## Getting Started
1. **Clone the repository**
2. **Set up the environment**
   - Use the provided `conda.env` file to create a conda environment:
     ```bash
     conda create -n mutsuite --file docs/conda.tx
     conda activate mutsuite
     ```
3. **Configure your simulation**
   - Edit a settings file (see `docs/simSettings.txt` for an example)
4. **Run the main notebook or scripts**
   - Example: Open and run `docs/demo.ipynb`

## Example Settings File
```
[Simulation]
depths: 30,100
pctMut: 0.05,0.1,0.15
addQual: 0
reps: 3
chr: chr2
loc: 72933869
simulateRange: 5000
FDRtolerance: 10,100,1000
reference: hg38/genome.fa
sourceBam: NA12878/HG001.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.300x.bam
alteredBam: SRR1046762.bam
useOnlyIndels: True
```

## Adding a New Caller
- Implement a new class inheriting from `Caller` (see `callers/pindelCaller.py` for an example)
- Add your caller to the list in the main notebook or script
