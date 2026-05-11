# PPanGGOLiN Command-Line Workflow

## Purpose

This workflow documents the command lines used to run PPanGGOLiN for the construction of the *Vibrio natriegens* pangenome from the RefSeq annotated genome records.

## Required software/environment

- Windows Subsystem for Linux
- Ubuntu/Linux terminal
- Conda
- PPanGGOLiN

## Input files

- RefSeq genome records with the `.gbff` format per strain placed within parent folders named after their repsective RefSeq strain assembly identifier.

## Output files

- `vibrio_genomes.list`: genome list file used as input for PPanGGOLiN.
- `vibrio_output/`: PPanGGOLiN output directory containing the pangenome results.

## Workflow

### 1. Navigate to the genome directory

```bash
cd /path/to/refseq/genome/records
```

### 2. Generate the genome list file

```bash
find "$(pwd)" -name "*.gbff" | awk -F/ '{print $(NF-1)"\t"$0}' > vibrio_genomes.list
```

This command searches recursively for `.gbff` files and creates a two-column genome list file. The first column contains the parent folder name and the second column contains the full path to the corresponding `.gbff` file.

### 3. Activate the PPanGGOLiN Conda environment

```bash
conda activate ppanggolin
```

### 4. Run the PPanGGOLiN pipeline

```bash
ppanggolin all --file vibrio_genomes.list --output vibrio_output
```

PPanGGOLiN automatically names the results folder in the following fomrat: ppanggolin_output_DATEYYYY-MM-DD_HOURHH.MM.SS_PIDXXXX

## Notes

The local directory path was replaced with a generic placeholder path. Users should update `/path/to/refseq/genome/records` according to their own local directory structure before running the workflow.