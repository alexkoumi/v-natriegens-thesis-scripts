# MMseqs2 Command-Line Workflow

## Purpose

This workflow documents the commands used to run MMseqs2 search on strain-specific *Vibrio natriegens* protein sequences against reference proteomes.

## Required software/environment

- Ubuntu terminal
- MMseqs2
- Standard Linux command-line tools: `bash`, `cat`, `mkdir`, `rm`

## Input files

- `query.fasta`: FASTA file containing the amino acid sequences of all strain-specific genes.
- Reference proteome FASTA files for the target database.

## Output files

- `queryDB`: MMseqs2 query database generated from `query.fasta`.
- `results/`: directory containing MMseqs2 search results for each processed target chunk.
- `result_*.m8`: tabular MMseqs2 alignment output files generated with `mmseqs convertalis`.

## Workflow

### Step 1: Create the MMseqs2 query database

```bash
mmseqs createdb /path/to/query.fasta queryDB
```

### 2. Define chunking parameters and output directories

```bash
CHUNK_SIZE=500
COUNTER=0
PART=1
TMP_DIR=tmp
OUT_DIR=results
mkdir -p "$TMP_DIR" "$OUT_DIR"
```

### 3. List all target reference proteome FASTA files

```bash
FILES=(/path/to/reference_proteomes/proteomes/*.fasta)
```

### 4. Process reference proteome FASTA files in chunks

```bash
for ((i=0; i<${#FILES[@]}; i++)); do
    cat "${FILES[i]}" >> "$TMP_DIR/target_chunk_$PART.fasta"
    let COUNTER+=1

    if [ "$COUNTER" -eq "$CHUNK_SIZE" ] || [ "$i" -eq $((${#FILES[@]} - 1)) ]; then
        echo "Processing chunk $PART..."

        mmseqs createdb "$TMP_DIR/target_chunk_$PART.fasta" "$TMP_DIR/targetDB_$PART"

        mmseqs search queryDB "$TMP_DIR/targetDB_$PART" "$OUT_DIR/result_$PART" "$TMP_DIR/tmp_$PART" --threads 4 -s 5.7 --max-seqs 300

        mmseqs convertalis queryDB "$TMP_DIR/targetDB_$PART" "$OUT_DIR/result_$PART" "$OUT_DIR/result_$PART.m8"

        rm "$TMP_DIR/target_chunk_$PART.fasta"
        mmseqs rmdb "$TMP_DIR/targetDB_$PART"
        rm -rf "$TMP_DIR/tmp_$PART"

        let PART+=1
        COUNTER=0
    fi
done
```

This loop concatenates reference proteome FASTA files into chunks of 500 files, creates a temporary MMseqs2 target database for each chunk, searches the query database against each target database, converts the results into tabular `.m8` format, and removes temporary files after each chunk is processed.

## Notes

The local directory paths were replaced with generic placeholder paths. Users should update `/path/to/query.fasta` and `/path/to/reference_proteomes/proteomes/*.fasta` according to their own local directory structure before running the workflow.
