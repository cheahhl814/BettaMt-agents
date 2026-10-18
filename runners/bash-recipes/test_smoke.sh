#!/usr/bin/env bash
# test_smoke.sh — bash-recipes structural smoke test.
#
# Checks (CPU-only, no data needed):
#   1. recipes parse (bash -n)
#   2. lock-step: every tool invocation in a module script block appears in the
#      mirroring recipe (this is the guard that caught the .py-entry-point and
#      pypolca _corrected.fasta drifts — a one-sided edit MUST fail here)
#   3. all tools referenced by recipes exist as deps in runners/nextflow-runner/pixi.toml
#   4. recipes use the pinned manifest, not bare PATH
set -u
cd "$(dirname "$0")"
ROOT=$(cd ../.. && pwd)                       # skill root
NF="$ROOT/runners/nextflow-runner"
ENV_MANIFEST="$NF/pixi.toml"
pass=0; fail=0
ok()  { echo "  PASS: $1"; pass=$((pass+1)); }
bad() { echo "  FAIL: $1"; fail=$((fail+1)); }

echo "== 1. bash -n parse check =="
for f in run_ont.sh run_illumina.sh annotate_genes.sh; do
  bash -n "$f" && ok "$f parses" || bad "$f: syntax error"
done

echo "== 2. module <-> recipe lock-step =="
# module -> recipe pairs; for each, pull "${PIXI} <cmd...>" lines out of the .nf
# script block, strip Groovy ${...} bits, and require the first two tokens
# (command [+ subcommand]) to appear in the recipe.
check_pair() {  # $1=module  $2=recipe
  local mod="$1" rec="$2"
  python3 - "$NF/modules/local/$mod" "$rec" <<'PYEOF'
import re, sys
mod, rec = open(sys.argv[1]).read(), open(sys.argv[2]).read()
missing = []
for line in re.findall(r'\$\{PIXI\}\s+(.+)', mod):
    raw = line.replace('\\', ' ')
    toks = [t for t in raw.split()
            if '${' not in t and not t.startswith('$')]
    toks = [t.rstrip(",\\;") for t in toks if t.strip(",\\;")]
    if not toks:
        continue
    cmd = toks[0]
    # tolerate bbduk.sh vs bbduk and bare vs .py entry point names
    def has(t):
        return (t in rec or t + '.py' in rec
                or t.replace('.py', '') in rec
                or t.replace('bbduk', 'bbduk.sh') in rec)
    missing = [] if has(cmd) else [cmd]
    for t in toks[1:]:
        if not has(t):
            missing.append(t)
    missing = [f"{mod}: '{t}'" for t in missing]
print('\n'.join(missing))
sys.exit(1) if missing else sys.exit(0)
PYEOF
  if [ $? -eq 0 ]; then ok "$mod <-> $(basename "$rec") lock-step"; else bad "$mod <-> $(basename "$rec") drift"; fi
}
check_pair bait_mito.nf          run_ont.sh
check_pair flye.nf               run_ont.sh
check_pair getorganelle_filter.nf run_ont.sh
check_pair firstgene.nf          run_ont.sh
check_pair circlator.nf          run_ont.sh
check_pair racon.nf              run_ont.sh
check_pair qc_short.nf           run_illumina.sh
check_pair getorganelle_asm.nf   run_illumina.sh
check_pair firstgene.nf          run_illumina.sh
check_pair circlator.nf          run_illumina.sh
check_pair polca.nf              run_illumina.sh
check_pair annotate_genes.nf     annotate_genes.sh

echo "== 3. recipe tools present in $ENV_MANIFEST deps =="
for tool in seqkit bbmap minimap2 samtools flye getorganelle circlator trnascan-se bedtools racon fastp pypolca blast python; do
  grep -qE "^\s*${tool}\s*=" "$ENV_MANIFEST" && ok "$tool pinned" || bad "$tool NOT in pixi.toml"
done

echo "== 4. recipes use the pinned manifest =="
for f in run_ont.sh run_illumina.sh annotate_genes.sh; do
  grep -q "pixi run --manifest-path" "$f" && ok "$f routes through pixi" || bad "$f: bare PATH invocation"
done

echo
echo "SMOKE TEST: $pass pass, $fail fail"
exit $([ "$fail" -eq 0 ] && echo 0 || echo 1)