#!/bin/bash
# yruler on one CRAM/BAM: count, then estimate
# usage: yruler-run.sh <cram|bam> <ref.fa> <markers.fa> <model.tsv> <name> <outdir> [regions.bed] [threads]
# regions.bed: read only these regions plus unmapped reads
set -uo pipefail
die(){ echo "ERROR yruler-run $*" >&2; exit 1; }
[ $# -ge 6 ] || die "usage: yruler-run.sh <cram|bam> <ref.fa> <markers.fa> <model.tsv> <name> <outdir> [regions.bed] [threads]"
IN=$1; REF=$2; CLS=$3; MOD=$4; NAME=$5; OUT=$6; BED=${7:-}; T=${8:-$(nproc)}
D=$(dirname "$(readlink -f "$0")"); Y=${YRULER:-$D/../target/release/yruler}
for f in "$IN" "$REF" "$CLS" "$MOD" ${BED:+"$BED"}; do [ -s "$f" ] || die "missing $f"; done
command -v samtools >/dev/null || die "samtools not in PATH"; [ -x "$Y" ] || die "no yruler binary at $Y (cargo build --release)"
mkdir -p "$OUT" || die "mkdir $OUT"; T0=$(date +%s)
# compressed markers: decompress once
case $CLS in *.zst) command -v zstd >/dev/null || die "zstd not in PATH"; zstd -dc "$CLS" > "$OUT/markers.fa" || die "zstd $CLS"; CLS=$OUT/markers.fa;; esac
echo "yruler-run $NAME input $IN bytes $(stat -c%s "$IN") markers $CLS model $MOD regions ${BED:-whole} threads $T start $(date)" >&2
RF='--input-fmt-option required_fields=0x21E'   # FLAG RNAME POS MAPQ SEQ
rd(){ if [ -n "$BED" ]; then grep -v '^\*' "$BED" > "$OUT/regions.tmp" && samtools view -@$T $RF -F 0x900 -M -L "$OUT/regions.tmp" -T "$REF" "$IN" && samtools view -@$T $RF -F 0x900 -T "$REF" "$IN" '*'; else samtools view -@$T $RF -F 0x900 -T "$REF" "$IN"; fi; }
rd 2> "$OUT/samtools.log" | "$Y" count "$CLS" - "$NAME" 500 $T "$OUT" 2> "$OUT/count.log"; P=("${PIPESTATUS[@]}")
[ "${P[0]}" = 0 ] || die "samtools exit ${P[0]}: $(head -2 "$OUT/samtools.log")"; [ "${P[1]}" = 0 ] || die "count exit ${P[1]}: $(tail -2 "$OUT/count.log")"
grep -qiE "error|fail|truncat|\[E::" "$OUT/samtools.log" && die "samtools: $(head -2 "$OUT/samtools.log")"
grep -E "^(markers|reads|check|gc|WARNING|y_dose)" "$OUT/count.log" | cut -c1-200 >&2
"$Y" estimate "$OUT/summary.tsv" "$MOD" > "$OUT/lengths.tsv" 2>> "$OUT/count.log" || die "estimate: $(tail -1 "$OUT/count.log")"
rm -f "$OUT/regions.tmp" "$OUT/markers.fa"; cat "$OUT/lengths.tsv"; echo "yruler-run $NAME seconds $(( $(date +%s)-T0 ))" >&2
