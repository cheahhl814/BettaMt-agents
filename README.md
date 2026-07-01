# BettaMt — betta fish mitogenome assembly pipeline + agentic skills

A Nextflow pipeline for assembling and polishing mitogenomes from long-read (ONT) and short-read (Illumina) sequencing data, paired with a set of agentic skills that operate around the pipeline.

## What's in this repo

| Path | What it is | Ship it? |
|---|---|---|
| `BettaMt/` | Nextflow DSL2 pipeline (pixi-managed) | ✅ |
| `.agents/skills/` | Agent Skills (SKILL.md) for pi, OpenCode, Codex | ✅ |
| `.gitignore` | Excludes Nextflow `work/`, `.nextflow/`, test data | ✅ |

## Quick start — pipeline

```bash
cd BettaMt
pixi install
pixi run nextflow run . -profile local --input samples.csv --outdir results/
```

See [`BettaMt/README.md`](BettaMt/README.md) for the full manual, profiles (local / slurm), and input format.

## Quick start — agentic skills

The four skills in `.agents/skills/` are auto-discovered by `pi`, **OpenCode**, and **Codex** (the standard `.agents/skills/` convention). They are not auto-discovered by **Claude Code** (which only reads `.claude/skills/`); add a symlink if needed.

Recommended entry point: **`bettamt-run`** — it detects the current stage and routes to:

- `bettamt-preflight` — produce `params.json` + rationale from FASTQ + reference
- `bettamt-debug` — diagnose a failed run, write `diagnosis.md` with top-3 hypotheses
- `bettamt-qc` — check the polished mitogenome, write `report.md` with verdicts

Each skill resolves `BETA_MT_HOME` (the pipeline directory) via `$BETA_MT_HOME` env var, or by looking for the sibling `../BettaMt` next to the skill folder.

## Design principle

> The agent reasons; the pipeline executes.

The pipeline is intentionally unaware of the agentic layer — and the agentic layer is intentionally unaware of the pipeline's internals. The handoff surface is a small, bounded JSON schema (`params.json` in, `diagnosis.md` / `report.md` out).

## License

Add a LICENSE file appropriate to your needs (MIT / Apache-2.0 / GPL-3.0) before public release.
