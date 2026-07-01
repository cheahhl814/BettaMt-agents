# BettaMt agentic skills

Companion skills for the BettaMt Nextflow pipeline. These define the agentic layer — pre-flight parameter preparation, post-run failure interpretation, and post-assembly QC. The pipeline itself is untouched; the skills consume the pipeline's inputs and outputs.

## Skills

| Skill | Purpose | When |
|---|---|---|
| `bettamt-run` | End-to-end orchestrator: detect stage, route to the right skill, hand off | User says "run BettaMt" / "assemble my mitogenome" |
| `bettamt-preflight` | Produce `params.json` + `params.rationale.md` from FASTQ + reference | Before `nextflow run` |
| `bettamt-debug` | Diagnose a failed run, write `diagnosis.md` with top-3 hypotheses | After a failure |
| `bettamt-qc` | Check the polished mitogenome, write `report.md` with verdicts | After success |

**Recommended entry point for new users:** `bettamt-run`. It detects the current stage (preflight / run / debug / qc) and routes to the appropriate skill.

## Pipeline location

Each skill resolves `BETA_MT_HOME` in this order:
1. `$BETA_MT_HOME` env var (set by user)
2. Sibling directory `../BettaMt` next to this skill folder
3. Asks the user

## Relocating

This folder is meant to live one directory level above the `BettaMt/` pipeline directory, so the working directory stays clean. Recommended layout:

```
bioinformatics/
├── .agents/skills/         # this folder
│   ├── bettamt-run/        # entry point — start here
│   ├── bettamt-preflight/
│   ├── bettamt-debug/
│   └── bettamt-qc/
└── BettaMt/                # pipeline (clean, shippable)
```

**Why `.agents/` (with the `s`)?** This is the cross-tool Agent Skills standard convention. `pi`, OpenCode, and Codex all auto-scan `.agents/skills/` from `cwd` up to the repo root. The singular `.agent/` is silently ignored by all of them. Claude Code is the one exception — it scans `.claude/skills/`; see anthropics/claude-code#33733.

To install globally for one user:

```bash
ln -s "$(pwd)" ~/.pi/agent/skills/bettamt-bundle
# or symlink each skill individually:
ln -s "$(pwd)/bettamt-run"       ~/.pi/agent/skills/
ln -s "$(pwd)/bettamt-preflight" ~/.pi/agent/skills/
ln -s "$(pwd)/bettamt-debug"     ~/.pi/agent/skills/
ln -s "$(pwd)/bettamt-qc"        ~/.pi/agent/skills/
```

## Extending

The four skills share a common principle: **the agent reasons; the pipeline executes**. When you encounter a new failure mode, add a row to the signature table in `bettamt-debug/SKILL.md` rather than writing a script. When you discover a new QC concern, add a check section to `bettamt-qc/SKILL.md`. The skill files are versioned alongside the pipeline (or separately, as you prefer).
