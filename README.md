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
yruler -i sample.cram -r GRCh38.fa -m data/markers.ms.fa.zst -M model.tsv -o outdir
```

The input needs its index (`.crai` or `.bai`).

`lengths.tsv` in `outdir` has the Y dosage and one length per array:

| Class | Array |
|---|---|
| `dyz3` | DYZ3 (Y centromere) |
| `yq12` | Yq12 |
| `dyz1` | DYZ1 units of Yq12 |
| `dyz2` | DYZ2 units of Yq12 |
| `dyz18s` | DYZ18 and 3.1 kb units |
| `yq2k7` | 2.7 kb units |
| `dyz18sub` | DYZ18, 3.1 kb and 2.7 kb units |
| `dyz19` | DYZ19 |

## Options

- `-t`: number of threads.
- `-n`: sample name.
- `-b data/regions_grch38_1kgp_dragen376.bed`: read only the regions with marker reads, for genomes aligned to GRCh38 with DRAGEN 3.7.6. For other aligners, run without `-b`.
