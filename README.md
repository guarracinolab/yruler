# yruler

Measure the human Y chromosome satellite arrays from sequencing reads.

yruler counts marker k-mers built from complete Y assemblies and reports, for each male:

- length of the Yq12 array;
- length of the DYZ3 centromeric array (active alpha-satellite HOR);
- Y dosage (2 x Y single-copy depth / chr2 depth).

Prototype: DYZ1, DYZ2, DYZ18, DYZ19 and the non-DYZ part of Yq12 are in development.

## Build

```bash
cargo build --release   # target/release/yruler; needs Rust >= 1.80 and libc
```

`scripts/yruler-run.sh` also needs `samtools` (>= 1.15) and `zstd` in `PATH`.

## Run

```bash
scripts/yruler-run.sh sample.cram GRCh38.fa data/markers_train.ms.fa.zst model.tsv SAMPLE outdir data/regions_grch38_1kgp_dragen376.bed 8
```

Inputs:

- CRAM or BAM aligned to GRCh38 (with `.crai`/`.bai`), and the reference FASTA;
- markers: `data/markers_train.ms.fa.zst` (4,730,117 bytes, 16,584,535 canonical 31-mers: dyz3 6,721, yq12 672,021, ysinglecopy 6,940,206, chr2depth 8,965,587), one masked superstring per class (KmerCamel; an uppercase base starts a k-mer). yruler also reads one 31-mer per line with headers `>class_N`;
- `model.tsv`: K per class;
- optional regions BED: only reads in these regions plus unmapped reads are read. `data/regions_grch38_1kgp_dragen376.bed` is only for genomes aligned with DRAGEN 3.7.6 (graph-based hg38). For other aligners read the whole file (no BED): on HG00126 aligned with BWA 0.7.15, the DRAGEN BED gave Yq12 18.61% short, the whole file 1.78% short (truth 24,631,978 bp).

Outputs in `outdir`:

- `lengths.tsv`: sample, class, length_bp, marker_sum_gc, ysinglecopy_sum, y_dose;
- `summary.tsv`: per class number of markers, count sum, median, fraction of zero counts, GC-corrected sum;
- `counts.tsv`: count of every marker k-mer;
- `bins.tsv`: reads and marker hits per 500 bp bin of the alignment;
- `any.fa`: reads with at least one marker k-mer.

The two steps can also run alone:

```bash
zstd -dc data/markers_train.ms.fa.zst > markers.fa
samtools view -F 0x900 -T GRCh38.fa sample.cram | yruler count markers.fa - SAMPLE 500 8 outdir
yruler estimate outdir/summary.tsv model.tsv
```

## Method

Length of an array = K x (GC-corrected sum of its marker counts) / (sum of the Y single-copy marker counts). K is fitted on males with complete Y assemblies. GC correction uses the depth curve by k-mer G+C of the Y single-copy markers.
