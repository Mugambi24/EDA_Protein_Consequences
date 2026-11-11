#!/bin/sh
set -eu
umask 022
mkdir -p /app/submission /app/workbench /app/output /logs/agent
LOG=/logs/agent/oracle-install.log
exec 2> "$LOG"
printf 'Oracle reference installation started\n'
for name in cli.py coding.py evidence.py audit.py inference.py; do
    src="/solution/$name"
    dst="/app/submission/$name"
    if [ ! -s "$src" ]; then printf 'ERROR missing source: %s\n' "$src"; exit 17; fi
    cp "$src" "$dst"
    chmod 644 "$dst"
    if ! cmp -s "$src" "$dst"; then printf 'ERROR byte mismatch: %s\n' "$name"; exit 18; fi
done
if grep -q 'NotImplementedError' /app/submission/coding.py; then
    printf 'ERROR attempted to submit incomplete coding.py\n'; exit 19
fi
printf 'Reference code copied byte for byte to /app/submission\n'
python3 -I -B /app/submission/cli.py \
    --reference /data/reference.json \
    --variants /data/practice.json \
    --observations /data/practice_observations.json \
    --fasta /data/NP_001390.1.faa \
    --calibration /data/spikein_calibration.json \
    --curated /data/public_eda_variants.json \
    --batch-assay /data/batch_assay.json \
    --phased-assay /data/phased_assay.json \
    --longitudinal-assay /data/longitudinal_assay.json \
    --phase-graph /data/phase_graph_assay.json \
    --output-dir /app/output
for result in coding.json evidence.json likelihood.json audit.json batch_inference.json phased_rna.json longitudinal.json phase_graph.json; do
    if [ ! -s "/app/output/$result" ]; then printf 'ERROR missing result: %s\n' "$result"; exit 20; fi
done
printf 'Oracle completed reference calculation, output and source ready for artifact collection\n'
