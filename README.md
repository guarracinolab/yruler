<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/yruler-logo-dark.svg">
    <img src="docs/yruler-logo.svg" alt="yruler" width="400">
  </picture>
</p>

<p align="center">Measure the human Y chromosome satellite arrays from sequencing reads.</p>

## Build

```bash
cargo build --release
```

## Run

```bash
# BAM
yruler -i sample.bam -m data/markers.ms.fa.zst -M model.tsv -o outdir

# CRAM: add the FASTA it was compressed with
yruler -i sample.cram -r GRCh38.fa -m data/markers.ms.fa.zst -M model.tsv -o outdir

# FASTQ (plain or gzip), pairs interleaved in one file
yruler -i sample.fastq.gz -m data/markers.ms.fa.zst -M model.tsv -o outdir

# FASTQ pairs in two files: list them with a comma
yruler -i R1.fastq.gz,R2.fastq.gz -m data/markers.ms.fa.zst -M model.tsv -o outdir

# SAM (plain or gzip)
yruler -i sample.sam -m data/markers.ms.fa.zst -M model.tsv -o outdir
```

With an index (`.bai` or `.crai`), yruler reads a BAM or CRAM file in parallel parts. Each read is counted on its own, so pairs give the same result interleaved or in two files. FASTQ input gives the same results as a BAM or CRAM file of the same reads.

`lengths.tsv` in `outdir` has the Y dosage and one value per class, in bp: the length of an array, or the amount of a repeat (all its units added up):

| Class | Value |
|---|---|
| `dyz3` | length of DYZ3 (Y centromere) |
| `yq12` | length of Yq12 |
| `dyz1` | amount of DYZ1 in Yq12 |
| `dyz2` | amount of DYZ2 in Yq12 |
| `dyz18s` | amount of DYZ18 and 3.1 kb repeats |
| `yq2k7` | amount of the 2.7 kb repeat |
| `dyz18sub` | amount of DYZ18, 3.1 kb and 2.7 kb repeats |
| `dyz19` | length of DYZ19 |

## Options

- `-r`: the FASTA that the CRAM was compressed with. Without `-r`, htslib looks for it through `REF_PATH`/`REF_CACHE` or downloads it by the MD5 in the CRAM header.
- `-t`: number of threads.
- `-n`: sample name.
- `-b data/regions_grch38_1kgp_dragen376.bed` (needs the index): read only the regions with marker reads, for genomes aligned to GRCh38 with DRAGEN 3.7.6. It needs the index (`.bai` or `.crai`). Without `-b`, yruler reads all reads.
- `--bins 500` (advanced): also write `bins.tsv`, with the reads and the marker hits of each class in every 500 bp of the reference genome. It needs aligned input (BAM, CRAM or SAM).
