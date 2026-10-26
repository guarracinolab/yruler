# yruler

Measure the human Y chromosome satellite arrays from sequencing reads.

For now yruler reports the length of the Yq12 and DYZ3 arrays and the Y dosage. More arrays are coming.

## Build

```bash
cargo build --release
```

## Run

```bash
yruler -i sample.cram -r GRCh38.fa -m data/markers_train.ms.fa.zst -M model.tsv -o outdir
```

The input needs its index (`.crai` or `.bai`). `lengths.tsv` in `outdir` has one row per array.

`-x` sets the sequencing technology: `illumina` (default), `element`, `hifi`, `ont` or `hic`. Each technology needs its own K in `model.tsv`; for now `model.tsv` has K for `illumina` only. `element` and `hic` are experimental.

With `-b data/regions_grch38_1kgp_dragen376.bed`, yruler reads only the regions that carry marker reads (genomes aligned to GRCh38 with DRAGEN 3.7.6). Regions for BWA and other DRAGEN versions are upcoming; for those, run without `-b`.
