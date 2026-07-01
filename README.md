# BettaMt — betta fish mitogenome assembly pipeline + agentic skills

A Nextflow pipeline for assembling and polishing mitogenomes from long-read (ONT) and short-read (Illumina) sequencing data, paired with a set of agentic skills that operate around the pipeline.

## What's in this repo

| Path | What it is | Ship it? |
|---|---|---|
| `BettaMt/` | Nextflow DSL2 pipeline (pixi-managed) | ✅ |
| `.agents/skills/` | Agent Skills (SKILL.md) for pi, OpenCode, Codex | ✅ |
| `CITATION.cff` | GitHub "Cite this repository" metadata | ✅ |
| `.gitignore` | Excludes Nextflow `work/`, `.nextflow/`, test data | ✅ |

## Quick start — pipeline

```bash
cd BettaMt
pixi install
pixi run nextflow run . -profile local --input samples.csv --outdir results/
```

See [`BettaMt/README.md`](BettaMt/README.md) for the full manual, profiles (local / slurm), and input format.

## The agentic interface

The `.agents/skills/` directory holds four SKILL.md files following the [Agent Skills standard](https://agentskills.io/specification). Each skill is a self-contained Markdown capability package: a YAML frontmatter (name + description) plus the body that an LLM agent reads at call time.

### How a tool decides to load a skill

The standard `description` frontmatter in each `SKILL.md` is what the agent sees at session start. A user prompt that matches a trigger phrase in the description is what causes the agent to `read` the rest of the file and behave as the skill prescribes. The skills are *model-invoked* — the agent chooses when to load them, not the user.

**Auto-discovered by** (project-local, scanned from `cwd` upward to the repo root):

- ✅ **pi** — this is the tool writing this README. Scans `.agents/skills/` after the project is marked trusted.
- ✅ **OpenCode** — scans `.agents/skills/` and walks up.
- ✅ **Codex** — scans `.agents/skills/` and walks up.
- ❌ **Claude Code** — does not scan `.agents/skills/` (see [anthropics/claude-code#33733](https://github.com/anthropics/claude-code/issues/33733)). If you also need Claude Code to see these, add a symlink in the project root: `ln -s .agents/skills .claude/skills`.

If you want them available across all your projects, copy or symlink the whole bundle to your user-scope path:

| Tool | User-scope path |
|---|---|
| pi | `~/.pi/agent/skills/` |
| OpenCode | `~/.config/opencode/skills/` or `~/.agents/skills/` |
| Codex | `~/.agents/skills/` |

### The four skills

| Skill | Frontmatter trigger | Purpose | Output |
|---|---|---|---|
| `bettamt-run` | *"run BettaMt"*, *"assemble my mitogenome"*, *"start a BettaMt run"*, *"bettamt run"* | End-to-end orchestrator: detect the user's stage and route to the right skill | Dispatches |
| `bettamt-preflight` | *"prepare params for BettaMt"*, *"preflight"*, *"set up a run"* | Produce `params.json` + `params.rationale.md` from FASTQ + reference | `params.json`, `params.rationale.md` |
| `bettamt-debug` | *"BettaMt failed"*, *"diagnose the run"*, *"why did my mitogenome assembly fail"* | Diagnose a failed run from `.nextflow.log` / `work/`, write `diagnosis.md` with top-3 hypotheses | `diagnosis.md` |
| `bettamt-qc` | *"check my polished mitogenome"*, *"QC the assembly"*, *"did it assemble correctly?"* | Check the polished mitogenome, write `report.md` with pass/warn/fail verdicts | `report.md` |

**Recommended entry point:** `bettamt-run`. It auto-detects the current stage (preflight / run / debug / qc) and routes.

### What a typical session looks like

```
You:    Run BettaMt on this betta mitogenome. I have ONT reads and a closely
        related reference.

Agent:  [loads .agents/skills/bettamt-run/SKILL.md]
        [detects STAGE=preflight, loads bettamt-preflight/SKILL.md]
        [asks 2-3 questions about your reads and reference]
        [writes params.json + params.rationale.md]

You:    Proceed.

Agent:  [loads bettamt-run again, detects STAGE=run, asks for confirmation]
        nextflow run main.nf -profile slurm -params-file params.json
        [runs it]

You:    It failed with "Bait_Mito exit code 1".

Agent:  [loads bettamt-debug/SKILL.md, inspects work/...]
        Wrote diagnosis.md — top hypothesis: bait database was indexed with the
        wrong lineage. Recommend Option A: re-run bait_mito with corrected
        taxonomy lookup.
```

### Where the agent writes outputs

By default the agent operates in `$BETA_MT_HOME` (the pipeline directory). All handoff files are written *next to* the pipeline, not inside it, so the pipeline working directory stays clean.

```
bioinformatics/
├── .agents/skills/        # the agentic layer
│   ├── bettamt-run/
│   ├── bettamt-preflight/
│   ├── bettamt-debug/
│   └── bettamt-qc/
├── BettaMt/               # the pipeline (Nextflow, pixi)
│   ├── main.nf
│   ├── conf/
│   ├── modules/
│   └── subworkflows/
├── params.json            # produced by bettamt-preflight
├── params.rationale.md    # produced by bettamt-preflight
├── diagnosis.md           # produced by bettamt-debug (only on failure)
└── report.md              # produced by bettamt-qc (only on success)
```

If you'd rather keep the run outputs in a separate directory, set `RUN_DIR` in the shell where the agent is running; the skills will write `params.json` / `diagnosis.md` / `report.md` there instead.

### The handoff contract

The skills share a small, bounded schema — that's the entire interface between agent and pipeline:

| Direction | Artifact | Schema |
|---|---|---|
| Agent → Pipeline | `params.json` | Nextflow `-params-file` JSON; `params.rationale.md` is the human-readable companion explaining every choice |
| Pipeline → Agent | `.nextflow.log`, `work/*/.command.*`, `results/**` | Nextflow-native files; the agent reads them, doesn't parse them |
| Agent → Human | `diagnosis.md`, `report.md` | Top-3-hypotheses / pass-warn-fail tables in the agentic-format style |

The pipeline never reads the Markdown files. The agent never reads `.nf` internals. The boundary is the JSON + the run log.

### Design principle

> The agent reasons; the pipeline executes.

When you encounter a new failure mode, add a row to the signature table in `bettamt-debug/SKILL.md` rather than writing a script. When you discover a new QC concern, add a check section to `bettamt-qc/SKILL.md`. The skill files are versioned alongside the pipeline and improved by the same person who runs the pipeline — that's how the agentic layer accumulates domain expertise.

### Customizing for your environment

- **`BETA_MT_HOME`** — if the pipeline is not at `../BettaMt` relative to the skills folder, set this in your shell: `export BETA_MT_HOME=/path/to/BettaMt`. All four skills honor it.
- **`RUN_DIR`** — defaults to `$BETA_MT_HOME`. Override to write `params.json` / `diagnosis.md` / `report.md` somewhere else.
- **Adding a new skill** — copy any of the four as a template, give it a unique `name` in the frontmatter, and write a clear `description` listing the trigger phrases. The agent will pick it up at next session start.

## Citing

If you use BettaMt in published research, see [`CITATION.cff`](CITATION.cff). The file is rendered by GitHub as a "Cite this repository" button. Fill in the placeholder author / ORCID fields before release.

## License

No license file is included; default copyright applies in jurisdictions that recognize it. If you intend this to be reusable, add a `LICENSE` file (MIT / Apache-2.0 / GPL-3.0) before pointing anyone at the repo.
