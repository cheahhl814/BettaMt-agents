# docs-corpus — offline CLI references (pinned versions)

Every tool invoked by the BettaMt runners has a verbatim `--help` capture taken
from the **pinned pixi env** (`runners/nextflow-runner/pixi.toml`), taken
2026-10-18 during the rebuild battle-test. Captures are authoritative for the
pinned versions — read these instead of guessing flags.

| Tool | Corpus entry | Captured invocation | Pin (pixi.toml) |
|---|---|---|---|
| seqkit | `seqkit/cli-reference.md` | `seqkit --help` | conda latest in lock |
| bbmap | `bbmap/cli-reference.md` | `bbduk.sh -h` | conda latest in lock |
| minimap2 | `minimap2/cli-reference.md` | `minimap2 --help` | conda latest in lock |
| samtools | `samtools/cli-reference.md` | `samtools --help` | `>=1.18` |
| flye | `flye/cli-reference.md` | `flye --help` | conda latest in lock |
| racon | `racon/cli-reference.md` | `racon --help` | conda latest in lock |
| bedtools | `bedtools/cli-reference.md` | `bedtools --help` | conda latest in lock |
| fastp | `fastp/cli-reference.md` | `fastp --help` | `>=0.24` |
| sequali | `sequali/cli-reference.md` | `sequali --help` | `>=0.14` |
| pypolca | `pypolca/cli-reference.md` | `pypolca --help` | `>=0.4.0` |
| blastn | `blastn/cli-reference.md` | `blastn --help` | `>=2.14` |
| circlator | `circlator/cli-reference.md` | `circlator --help` | conda latest in lock |
| getorganelle | `getorganelle/cli-reference.md` | `get_organelle_from_reads.py -h` | `>=1.7.7.1,<2` |
| trnascan-se | `trnascan-se/cli-reference.md` | `tRNAscan-SE -h` | conda latest in lock |

Refresh a capture after any manifest pin change:
```bash
pixi run --manifest-path runners/nextflow-runner/pixi.toml <tool> --help \
    > docs-corpus/<tool>/cli-reference.md
```

## Battle notes — manual walkthrough findings (2026-10-18)

The rebuild's manual walkthrough (CPU-only, synthetic reads off NC_002333.2)
executed all ONT stages + all Illumina stages against the pinned env. Bugs found
and fixed in the SAME commit (recipes + modules in lock-step per
`runners/bash-recipes/test_smoke.sh`):

1. `fastp` and `sequali` were absent from the runner manifest — QC_SHORT
   resolved a host-global fastp (which read 0 reads and silently emitted empty
   outputs on this host) and would fail outright on hosts without one. Both now
   pinned.
2. `getorganelle_filter.nf` called bare `get_organelle_config` /
   `get_organelle_from_assembly`; the pinned env ships only `.py`-suffixed
   entry points. Module + recipe now use `.py` names.
3. `bait_mito.nf` called bare `bbduk`; the env installs `bbduk.sh` only.
   Module + recipe now use `bbduk.sh`.
4. `polca.nf` copied `pypolca_out/<prefix>.fasta`; pypolca 0.4+ emits
   `<prefix>_corrected.fasta`. Module + recipe updated.
5. `sequali --force` does not exist in the pinned sequali; flag removed from
   module + recipe.

Smoke-data lesson: synthetic reads with zero errors, uniform coverage and
constant `~I` quality break GetOrganelle's graph disentangling ("Incomplete/
Complicated graph: EDGE_nn" → 279 bp scaffold even at `-R 15`). Realistic
per-base quality (Q20–40) + indel error profile resolved it (16,583 bp path at
`-R 1`). Don't judge the pipeline on perfect synthetic data.