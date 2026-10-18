---
name: bettamt
description: Agent skill bundle for the BettaMt mitogenome workflow (vertebrate mitogenome assembly from ONT long reads or Illumina paired short reads, with optional BLAST-based gene annotation). Preflight produces params.json + params.rationale.md with pinned-tool evidence; the default runner is manual bash recipes (single machine, per-stage headers); an opt-in Nextflow runner wraps the same commands for cluster/resume use; qc verifies contiguity, circularity, trnF/tRNA set, coverage and NUMT risk; debug is a signature-library failure interpreter. Use when the user provides raw mito reads plus a reference mitogenome and asks to assemble, annotate, or QC a mitochondrial genome.
version: 1.0.0
updated: "2026-10-18"
triggers:
  - "assemble mitogenome"
  - "BettaMt"
  - "bettamt"
  - "mitogenome from ONT reads"
  - "mitogenome from Illumina reads"
  - "getorganelle assembly"
  - "mito genome polish"
  - "circlator fixstart"
---

# bettamt — BettaMt mitogenome workflow (skill router)

> **v1.0.0** — rebuilt to the skill-forge convention: root router + phase sub-skills + dual runners with **manual bash recipes as the default source of truth** and the Nextflow pipeline as an opt-in executor (`runners/nextflow-runner/`). The recipe↔module lock-step is enforced by `runners/bash-recipes/test_smoke.sh` (31 checks, 0 fail at ship time — it is what caught the `.py`-entry-point, `bbduk.sh` and `_corrected.fasta` drifts during the CPU walkthrough).

## Audience

Two audiences, like every skill-forge bundle:

1. **AI agents** — follow the phase table below in order; every phase has Inputs/Outputs contracts and ask-user stop points; the filesystem is the handoff, never agent memory.
2. **Human operators** — each phase SKILL.md doubles as an operator runbook with copy-pastable commands.

## Layout

```
bettamt/  (this repo = skill root; $BETTAMT_HOME defaults here)
├── SKILL.md                      # this router — entry point
├── README.md                     # human-facing overview
├── CITATION.cff                  # citation metadata (upstream, kept)
├── pixi.toml                     # root workspace: run tasks (update-check, smoke)
├── bin/
│   └── skill-update-check.py     # deployed-copy vs upstream SHA diff
├── preflight/bettamt-preflight/  # phase 1: params.json + rationale + inventory
├── run/bettamt-run/              # phase 2: execution-mode fork → recipes (default) | NF runner
├── qc/bettamt-qc/                # phase 3: assembly QC + annotation QC (folded from bettamt-annotate-qc)
├── debug/bettamt-debug/          # phase 4: signature-library failure interpreter
├── runners/
│   ├── bash-recipes/             # DEFAULT runner (manual, single machine)
│   │   ├── run_ont.sh
│   │   ├── run_illumina.sh
│   │   ├── annotate_genes.sh
│   │   ├── test_smoke.sh         # structural + module-recipe lock-step checks
│   │   └── README.md
│   └── nextflow-runner/          # OPT-IN runner (same env, cluster/resume)
│       ├── main.nf + subworkflows/ + modules/ + conf/
│       ├── pixi.toml + pixi.lock # the pinned tool env BOTH runners share
│       └── bin/blast_extract_cds.py
└── docs-corpus/                  # offline tool-doc snapshots (fight flag hallucination)
```

## Execution model (ask at preflight; default `manual`)

| Mode | When | Driver |
|---|---|---|
| **manual** (default) | single sample, workstation, interactive stage-visibility | `runners/bash-recipes/*.sh` — per-stage `[stage k/N]` headers, `set -euo pipefail`, all tools via the pinned pixi env |
| **nextflow** (opt-in) | cluster/SLURM, `--resume` heavy iteration, trace reports | `nextflow run runners/nextflow-runner/main.nf -profile <env> -c <override> -params-file params.json` |

The recipes mirror the NF module scripts **command-for-command**; `bash-recipes/test_smoke.sh` fails the build if the two drift apart. One pinned env (`runners/nextflow-runner/pixi.toml`) serves both modes — never fork the manifest per mode.

## Phases (invoke the phase's SKILL.md, not this router)

| # | Phase | Sub-skill | Produces | Consumes |
|---|---|---|---|---|
| 1 | Preflight | `preflight/bettamt-preflight/SKILL.md` | `params.json`, `params.rationale.md`, optional `nextflow.config.<env>.local` (NF mode) | user's reads + reference mitogenome |
| 2 | Run | `run/bettamt-run/SKILL.md` | polished `$RUN_DIR/<sample>_mitogenome.fasta` (+ stage logs) | params.json, execution mode |
| 3 | QC | `qc/bettamt-qc/SKILL.md` | `report.md` (contiguity, circularity, tRNA set, coverage/NUMT, identity) + `annotation-report.md` if annotated | run outputs |
| 4 | Debug | `debug/bettamt-debug/SKILL.md` | `diagnosis.md` (signature match + fix) | failing stage/process stderr |

## Stop-point discipline

Every sub-skill encodes its ask-user stop points as `SP<n>` tables using **Evidence → Recommend → Options**. Do not auto-pick what a stop point owns (platform fork, taxon model, execution mode, annotation route); do not ask what preflight already recorded in `params.json`.

## External dependencies

**[read-qc-trimming](https://github.com/cheahhl814/read-qc-trimming)** (v6+; deployed at `~/.agents/skills/read-qc-trimming/`) — the agent-side read QC/trimming sandwich (sequali + multiqc raw QC → fastp / fastplong / trim_galore → post-trim QC → Go/No-Go loop).

- **When to delegate**: preflight evidence shows adapter contamination, 3'-end quality decay, poly-G tails, degenerate/short-fragment libraries, or the user wants kit-specific stringency or batch-wide multiqc reporting. `preflight/bettamt-preflight` proposes the routing (SP — Evidence → Recommend → Options); `run/bettamt-run` §3c carries it out.
- **Default remains in-run**: the pinned `QC_SHORT` stage (`sequali` + `fastp` via the shared pixi manifest) is always available and is what the smoke tests verify — the bundle degrades gracefully when the external skill is not deployed.
- **ONT caveat**: delegating pre-cleaning for ONT reads (read-qc-trimming's `fastplong`) shortens reads *before* BettaMt's 14–20 kb bait window — user-directed only, flagged by preflight.

## Self-update check

```bash
pixi run update-check        # from the repo or the deployed copy
```

`bin/skill-update-check.py` diffs the deployed copy's git SHA against `origin/main` (the source repo is `github.com/cheahhl814/BettaMt-agents`). Verdicts: `UP-TO-DATE` / `LOCAL-AHEAD` / `BEHIND-BY-N` / `OFFLINE` / `NO-ORIGIN`.

## What NOT to do

- Do not run the Nextflow runner by default on a single machine — the manual recipes are the default; NF mode is opt-in.
- Do not edit `runners/nextflow-runner/conf/*.config` in place (user-side `-c` override only), and never diverge recipes from modules without updating both in the same commit.
- Do not invent tool flags — check `docs-corpus/<tool>/` first; it exists because tools drift.
- Do not run `annotate_genes.sh` without a reference annotation (GFF3 or GenBank) — that input gate is deliberate.
- No GPU anywhere in this workflow: Flye/racon/pypolca/GetOrganelle are CPU tools; never recommend GPU.