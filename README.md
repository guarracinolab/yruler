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

With an index (`.bai` or `.crai`), yruler reads a BAM or CRAM file in parallel parts. Each read is counted on its own, so pairs give the same result interleaved or in two files.

`lengths.tsv` in `outdir` has one value per class, in bp: the length of an array, or the amount of a repeat (all its units added up):

| Class | Value |
|---|---|
| `dyz3` | length of DYZ3 (Y centromere) |
| `yq12` | length of Yq12 (DYZ1 + DYZ2 + the rest of Yq12) |
| `dyz1` | amount of DYZ1 in Yq12 |
| `dyz2` | amount of DYZ2 in Yq12 |
| `dyz18s` | amount of DYZ18 and 3.1 kb repeats, already inside Yq12: do not add it to `yq12` |
| `dyz19` | length of DYZ19 |
| `dyz17` | length of DYZ17 |

The other columns:

- `relative`: the class marker signal divided by the Y single-copy (ruler) signal, without the model scale factor; `length_bp` is this value times the class factor of the model. For `yq12`, `length_bp` is the sum of the `dyz1`, `dyz2` and rest-of-Yq12 lengths, and `relative` is the signal of the separate Yq12 marker class. It does not depend on sequencing depth or on the fraction of cells with a Y.
- `y_dose`: Y copies relative to a typical male (1 = one Y in every cell), from the Y single-copy k-mers against autosomal reference k-mers.
- `estimate_basis`: `ruler` when `y_dose` is above the no-Y cut of the model (the highest `y_dose` of the females in the model cohort); with too few Y reads for a G+C curve, these values have no G+C correction and yruler prints a warning. `haploid_depth` when `y_dose` is at or below the cut (females, samples without Y). Then the values are background signal in bp at haploid autosomal depth, not array lengths.

## Options

- `-r`: the FASTA that the CRAM was compressed with. Without `-r`, htslib looks for it through `REF_PATH`/`REF_CACHE` or downloads it by the MD5 in the CRAM header.
- `-t`: number of threads.
- `-n`: sample name.
- `-x illumina` (default): the technology of the model. The released model is calibrated on Illumina short reads only, and yruler does not detect other technologies: do not use it for other reads.
- `-b data/regions_grch38_1kgp_dragen376.bed`: read only the regions with marker reads and the autosomal reference window (chr8:91-101 Mb). Use it only for genomes aligned to GRCh38 with DRAGEN 3.7.6: other aligners place Yq12 reads elsewhere, and their Yq12 and DYZ2 values come out too low. It needs the index (`.bai` or `.crai`). Without `-b`, yruler reads all reads.
- `--reads marker_reads.fa` (advanced): also write the sequences of the reads with a marker hit. yruler does not write read sequences otherwise, so the output holds counts only. Do not use this option on controlled-access data unless your data agreement allows it.
- `--bins 500` (advanced): also write `bins.tsv`, with the reads and the marker hits of each class in every 500 bp of the reference genome. It needs aligned input (BAM, CRAM or SAM).
