# Convert AXELIOS BAM to FASTQ

Use this procedure when an AXELIOS alignment workflow produces a BAM file that
must be converted back to FASTQ for another alignment or analysis step.

{% stepper %}
{% step %}

## Requirements

The input BAM must contain the `YC` tag for duplex conversion.
Simplex conversion does not require this tag.

{% endstep %}
{% step %}

## Duplex reads

For duplex reads, preserve the `YC` tag in each FASTQ header.
The `tr` command changes the tab separator emitted by `samtools fastq -T` to a
space so that the tag remains part of the read header.

```bash
samtools fastq -@ "${SAMTOOLS_THREADS}" -T YC "${INPUT_BAM}" | tr '\t' ' ' | pigz -p "${PIGZ_THREADS}" > "${OUTPUT_FQ}"
```

{% endstep %}
{% step %}

## Simplex reads

For simplex reads, the `YC` tag is not needed, so omit both `-T YC` and `tr`:

```bash
samtools fastq -@ "${SAMTOOLS_THREADS}" "${INPUT_BAM}" | pigz -p "${PIGZ_THREADS}" > "${OUTPUT_FQ}"
```

{% endstep %}
{% step %}

## Shuffle FASTQ for Giraffe personalization

When the FASTQ will be used for Giraffe alignment with pangenome
personalization, shuffle the reads before building the k-mer profile.
The example below splits the input into ten chunks to limit the memory required
by each shuffle operation.

### Step 1: Split the input FASTQ

```bash
echo "Step 1/3: Splitting FASTQ file into chunks..."
seqkit split2 \
    --threads "$SEQKIT_THREADS" \
    -p 10 \
    -O "$CHUNKS_DIR" \
    "$INPUT_FASTQ"
```

### Step 2: Shuffle each chunk

Use a fixed seed when reproducibility is important.
The loop skips completed outputs and temporary outputs left by an interrupted
run, so it can be rerun safely against the same chunk directory.

```bash
set -euo pipefail

echo "Step 2/3: Shuffling chunks in memory..."
for chunk in "$CHUNKS_DIR"/*; do
    [ -f "$chunk" ] || continue
    [[ "$chunk" == *.shuf.fastq.gz || "$chunk" == *.shuf.fastq.gz.tmp ]] && continue  # Skip processed chunks

    echo "Shuffling: $(basename "$chunk")"
    tmp_output="${chunk}.shuf.fastq.gz.tmp"
    rm -f -- "$tmp_output"

    if seqkit shuffle \
        --threads "$SEQKIT_THREADS" \
        -s 42 \
        "$chunk" \
        -o "$tmp_output"; then
        mv -- "$tmp_output" "${chunk}.shuf.fastq.gz"
    else
        rm -f -- "$tmp_output"
        echo "Shuffling failed: $(basename "$chunk")" >&2
        exit 1
    fi
done
```

The loop writes one `.shuf.fastq.gz` file for each input chunk.

### Step 3: Merge the shuffled chunks

Concatenate the shuffled chunks into the FASTQ used by the personalization
workflow.

```bash
echo "Step 3/3: Merging shuffled chunks..."
cat "$CHUNKS_DIR"/*.shuf.fastq.gz > "$OUTPUT_FASTQ"
```

{% endstep %}
{% endstepper %}
