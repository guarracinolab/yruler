<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/yruler-logo-dark.svg">
    <img src="docs/yruler-logo.svg" alt="yruler" width="400">
  </picture>
</p>

<p align="center">Measure the human Y chromosome satellite arrays from sequencing reads.</p>

## Install

```bash
git clone https://github.com/guarracinolab/yruler
cd yruler
cargo build --release   # binary: target/release/yruler
```

## Demo

The demo streams the public 1000 Genomes DRAGEN 3.7.6 alignment of HG00126 from Amazon S3 and reads only the marker regions (`-b`), so it needs no download beforehand:

```bash
yruler -i https://s3.us-west-2.amazonaws.com/1000genomes-dragen-3.7.6/data/individuals/hg38-graph-based/HG00126/HG00126.bam \
  -b data/regions_grch38_1kgp_dragen376.bed -m data/markers.ms.fa.zst -M model.tsv -n HG00126 -t 8 -o demo
```

Expected `demo/estimates.tsv` (columns `sample`, `class`, `estimate_bp`, `y_dose`, `y_status`):

| class | estimate_bp | y_dose | y_status |
|---|---|---|---|
| dyz3 | 309126 | 0.9698 | y_found |
| yq12 | 24159198 | 0.9698 | y_found |
| dyz1 | 13477259 | 0.9698 | y_found |
| dyz18s | 301675 | 0.9698 | y_found |
| dyz19 | 408414 | 0.9698 | y_found |
| dyz17 | 838697 | 0.9698 | y_found |
| dyz2 | 10119601 | 0.9698 | y_found |

## Run

```bash
yruler -i sample.bam -m data/markers.ms.fa.zst -M model.tsv -o outdir
```

yruler counts the k-mers of each array in the reads and compares them with the k-mers of Y regions that have one copy (the ruler). A longer array gives a larger count.

yruler writes `outdir/estimates.tsv` and prints it. The file has one row per reported class (7 rows) and 9 columns:

| Column | Value |
|---|---|
| `sample` | sample name (`-n`) |
| `class` | the array or repeat (table below) |
| `estimate_bp` | the length of the array, or the sum of all units of the repeat, in bp |
| `relative` | `marker_sum_gc` / `ysinglecopy_sum` ([Details](#how-estimates-are-computed)); `.` for `yq12`, a sum of classes |
| `marker_sum_gc` | count of the class, G+C corrected; `.` for `yq12` |
| `ysinglecopy_sum` | count of the ruler |
| `y_dose` | how much Y the sample has, same value in all rows |
| `y_status` | `y_found` or `no_y`, same value in all rows |
| `terms` | the classes whose lengths make `estimate_bp` |

| Class | `estimate_bp` |
|---|---|
| `yq12` | Yq12 array: DYZ1 + DYZ2 + DYZ18S + the rest of Yq12 |
| `dyz1` | DYZ1 in Yq12 |
| `dyz2` | DYZ2 in Yq12 |
| `dyz18s` | DYZ18 and 3.1 kb repeats, a part of Yq12 |
| `dyz3` | DYZ3 array (Y centromere) |
| `dyz17` | DYZ17 array |
| `dyz19` | DYZ19 array |

Before you use `estimate_bp`, read `y_status` and `y_dose`:

- `y_status` = `y_found`: yruler found a Y, and `estimate_bp` is the array length.
- `y_status` = `no_y`: yruler found no Y (for example, a female), and `estimate_bp` is `NA`. `relative` and `marker_sum_gc` still hold the background counts.
- `y_dose` = 1: one Y in every cell, as in a typical male. 0.4: a Y in 40% of the cells. 0: no Y. If `y_dose` is below 0.1, or the Y coverage is below 0.5x, the estimates have large errors. yruler does not look for deletions, so `y_dose` = 1 does not mean that the whole Y is present.

`model.tsv` is calibrated on PCR-free Illumina 2x150 bp reads from the 1000 Genomes high-coverage project. PCR libraries, other read lengths and other projects shift the values, so compare only samples with the same protocol.

### Nanopore reads

`data/ont/` holds a second model, for nanopore R10.4.1 reads called with dorado in super-accuracy mode. It has its own marker file, so pass both files and `-x ont`:

```bash
yruler -i sample.bam -x ont -m data/ont/markers.ms.fa.zst -M data/ont/model.tsv -o outdir
```

Use the files of one model together. The Illumina model and the nanopore model are calibrated apart, so compare nanopore samples only with other nanopore samples.

## Input

```bash
# CRAM: add the FASTA it was compressed with
yruler -i sample.cram -r GRCh38.fa -m data/markers.ms.fa.zst -M model.tsv -o outdir

# FASTQ (plain or gzip), pairs interleaved in one file
yruler -i sample.fastq.gz -m data/markers.ms.fa.zst -M model.tsv -o outdir

# FASTQ pairs in two files: list them with a comma
yruler -i R1.fastq.gz,R2.fastq.gz -m data/markers.ms.fa.zst -M model.tsv -o outdir

# SAM (plain or gzip), or SAM on stdin with -i -
samtools view -h sample.bam | yruler -i - -m data/markers.ms.fa.zst -M model.tsv -o outdir
```

yruler also reads FASTA. Each read counts on its own, so pairs give the same result interleaved or in two files. With an index (`.bai` or `.crai`), yruler reads parts of a BAM or CRAM in parallel. If the file has bases written as `=` (same as the reference), yruler stops and writes no output.

## Options

- `-r`: the FASTA the CRAM was compressed with. Without `-r`, htslib uses `REF_PATH`/`REF_CACHE` or downloads the FASTA by the MD5 in the CRAM header.
- `-t`: number of threads (default: all CPUs).
- `-n`: sample name (default: the input file name up to the first `.`).
- `-b data/regions_grch38_1kgp_dragen376.bed`: read only the marker regions, the unmapped reads and the autosomal window chr8:91-101 Mb. It needs the index. Use it only for GRCh38 alignments made with DRAGEN 3.7.6: other aligners put Yq12 reads elsewhere, and Yq12 and DYZ2 come out too low.
- `-M`: the model. Without `-M`, yruler writes only `counts.tsv` and `summary.tsv` (see `yruler estimate` below).
- `-x illumina` (default): the technology of the model. yruler does not detect the technology. `model.tsv` has K values for `illumina` and `data/ont/model.tsv` has them for `ont` ([Nanopore reads](#nanopore-reads)); `hifi` and `hic` are experimental and no released model has them, so yruler stops with an error for those.

## Details

### Reads counted

yruler counts every primary record, duplicates and reads that fail quality checks (0x200) included, and skips secondary and supplementary records. It counts the 31-mers of each read that match the marker k-mers in `data/markers.ms.fa.zst`. The model is calibrated on files that keep duplicate reads, so do not remove duplicates before you run yruler: files without them give shorter DYZ2 and Yq12 lengths.

### How estimates are computed

The marker file has one k-mer set per class: the classes above (for `yq12`, all Yq12 markers, counted but not used for `estimate_bp`), `ysinglecopy` (Y single-copy k-mers, the ruler) and `depth1` (autosomal k-mers in the chr8 window).

1. G+C correction: yruler fits a curve of ruler counts against k-mer G+C and divides the marker counts by it. This removes the G+C difference between markers and ruler; the ruler is the reference, so its count stays as it is. With too few ruler counts, there is no curve, `sum_gc` = `sum`, and yruler prints a warning.
2. `y_dose` = 2 × mean ruler count / mean `depth1` count, divided by `y_dose_scale` in the model (the median `y_dose` of the 1000 Genomes samples above half the cohort median). It is printed with 4 decimals.
3. `y_status`: `y_found` if `y_dose` is above `no_y_dose` in the model (0.025 in `model.tsv` and `data/ont/model.tsv`, above the highest `y_dose` of the females in 1000 Genomes, HGDP, SGDP and the nanopore samples), else `no_y`. For `no_y` samples `estimate_bp` is `NA`.
4. `relative` = `marker_sum_gc` / `ysinglecopy_sum` for every sample, so it does not change with sequencing depth. If `ysinglecopy_sum` is 0, `relative` and `estimate_bp` are `NA`.
5. `estimate_bp` = K of the class in the model × `relative`. For `yq12`, `estimate_bp` is the sum of the `dyz1`, `dyz2`, `dyz18s` and `yq12rest2` lengths (`yq12rest2`: the rest of Yq12, computed but not reported on its own), as the `length` column of `model.tsv` and the `terms` column say.

### Output files

- `estimates.tsv`: see [Run](#run).
- `summary.tsv`: one row per marker class with `nreads` (all reads of the sample), `n_kmers`, `sum`, `median`, `frac_zero` (fraction of k-mers with count 0) and `sum_gc`. The first line holds the `y_dose` before the model scale, `ysinglecopy_depth`, `depth_mean` and `gc_curve` (1 if the G+C curve was fitted).
- `counts.tsv`: the count of each marker k-mer.

`yruler estimate outdir/summary.tsv model.tsv [-x tech]` computes `estimates.tsv` again from `summary.tsv`, to stdout, without the reads. Use it to apply a different model.

### Optional output

- `--reads marker_reads.fa`: also write the reads with a marker hit, with the header `contig|bin|MAPQ|flag`. Without it, yruler writes no read sequences, only counts. On controlled-access data, use it only if your data agreement allows it.
- `--bins 500`: also write `outdir/bins.tsv` with one row per 500 bp bin of the reference that has reads: `contig`, `start`, `reads`, `reads_mq20` (MAPQ 20 or more), `duplicates`, and for each class `<class>_hits`, `<class>_hits_mq20` and `<class>_reads`. It needs aligned input (BAM, CRAM or SAM).

## License

MIT, see `LICENSE`.
