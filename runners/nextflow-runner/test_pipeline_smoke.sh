#!/usr/bin/env bash
# test_pipeline_smoke.sh — Nextflow runner structural smoke (CPU-only).
#
# Builds the full DAG with -preview against dummy inputs: validates main.nf,
# subworkflows/, modules/, param validation logic and the user-config override
# layering without executing any pipeline tool. NOTE: -stub-run would still RUN
# the real scripts here — it stubs only processes that declare a stub: block,
# and these local modules have none (nextflow semantics). -preview is the pure
# DAG build. Tool-level execution is proven by the runners/bash-recipes manual
# walkthrough + test_smoke.sh lock-step.)
set -u
cd "$(dirname "$0")"
ROOT=$(cd ../.. && pwd)                       # skill root
ENV_MANIFEST="envs/nextflow/pixi.toml"
NF() { pixi run --manifest-path "$ENV_MANIFEST" nextflow "$@"; }
pass=0; fail=0
ok()  { echo "  PASS: $1"; pass=$((pass+1)); }
bad() { echo "  FAIL: $1"; fail=$((fail+1)); }

TMP=$(mktemp -d)

echo "== 1. dummy input generation =="
python3 - "$TMP" <<'PYEOF'
import sys, random
random.seed(7)
p = sys.argv[1]
seq = "".join(random.choice("ACGT") for _ in range(16596))
with open(f"{p}/ref.fasta", "w") as f:
    f.write(">dummy_ref_mitogenome\n")
    for i in range(0, len(seq), 60):
        f.write(seq[i:i+60] + "\n")
for name, n in (("ont.fastq", 20), ("r1.fastq", 100), ("r2.fastq", 100)):
    with open(f"{p}/{name}", "w") as f:
        for i in range(n):
            s = "".join(random.choice("ACGT") for _ in range(200))
            f.write(f"@d_{i}\n{s}\n+\n{'I'*len(s)}\n")
PYEOF
for f in "$TMP"/ont.fastq "$TMP"/r1.fastq "$TMP"/ref.fasta; do [ -s "$f" ] || { echo "dummy generation failed: $f"; exit 1; }; done
ok "dummy ref + fastqs"

echo "== 2. params schemas =="
python3 - "$TMP" <<'PYEOF'
import json, sys
p = sys.argv[1]
base = {"execution_mode": "nextflow", "ref_mito": f"{p}/ref.fasta", "size": "16k",
        "taxon": "vertebrate", "rounds": 15, "sample_id": "smoke",
        "ref_gff": "", "ref_gff_fasta": "", "ref_gb": "", "trnascan": "", "skip_trna": False}
json.dump({**base, "platform": "ont", "reads": f"{p}/ont.fastq"}, open(f"{p}/ont.json", "w"))
json.dump({**base, "platform": "illumina", "reads": f"{p}/r1.fastq", "reads_r2": f"{p}/r2.fastq"},
          open(f"{p}/illumina.json", "w"))
PYEOF
[ -s "$TMP/ont.json" ] && [ -s "$TMP/illumina.json" ] && ok "params files written" || bad "params generation"

echo "== 3. NF mode: ONT DAG (preview) =="
if NF run main.nf -profile local -preview -params-file "$TMP/ont.json" \
     -work-dir "$TMP/work_ont" --outdir "$TMP/out_ont" > "$TMP/ont.log" 2>&1; then
  ok "ONT stub run"
else
  tail -12 "$TMP/ont.log"; bad "ONT stub run aborted"
fi

echo "== 4. NF mode: Illumina DAG (preview) =="
if NF run main.nf -profile local -preview -params-file "$TMP/illumina.json" \
     -work-dir "$TMP/work_il" --outdir "$TMP/out_il" > "$TMP/il.log" 2>&1; then
  ok "Illumina stub run"
else
  tail -12 "$TMP/il.log"; bad "Illumina stub run aborted"
fi

echo "== 5. user-config override layering =="
python3 - "$TMP" <<'PYEOF'
import sys
open(sys.argv[1] + "/override.local", "w").write(
    "profiles { local { process { withName: 'FIRSTGENE' { cpus = 4; memory = '8 GB' } } } }\n")
PYEOF
if NF run main.nf -profile local -preview -params-file "$TMP/illumina.json" \
     -c "$TMP/override.local" -work-dir "$TMP/work_ov" --outdir "$TMP/out_ov" > "$TMP/ov.log" 2>&1; then
  ok "user override accepted via -c"
else
  tail -8 "$TMP/ov.log"; bad "user override rejected"
fi

rm -rf "$TMP"
echo
echo "PIPELINE SMOKE TEST: $pass pass, $fail fail"
exit $([ "$fail" -eq 0 ] && echo 0 || echo 1)