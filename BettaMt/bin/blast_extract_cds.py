#!/usr/bin/env python3
"""
blast_extract_cds.py — Homology-based full mitogenome annotation

Annotates all genes in a target mitogenome FASTA using a reference
annotation (GFF3 + FASTA or GenBank) and tRNAscan-SE:

  * 13 PCGs + 2 rRNAs (12S, 16S) — BLASTN against reference gene sequences
  * 22 tRNAs                      — tRNAscan-SE (-M vert, vertebrate mito)
  * D-loop (control region)       — largest intergenic gap detection

Outputs three files sorted by genomic position:
  * FASTA  (.fasta)  — nucleotide sequences
  * GFF3   (.gff3)   — standard gene annotations
  * BED    (.bed)    — genomic coordinates

Gene names are compatible with core_mito_phylo.py and MitoAnnotator format.

Usage:
  # Using GFF3 + reference FASTA
  python blast_extract_cds.py \\
      --target assembly.fasta \\
      --ref-gff annotation.gff3 --ref-fasta ref.fasta \\
      --out annotated_genes.fasta

  # Using GFF3 with embedded FASTA (##FASTA section)
  python blast_extract_cds.py \\
      --target assembly.fasta \\
      --ref-gff annotation_with_seqs.gff3 \\
      --out annotated_genes.fasta

  # Using GenBank
  python blast_extract_cds.py \\
      --target assembly.fasta \\
      --ref-gb reference.gb \\
      --out annotated_genes.fasta

  # Optional tRNAscan-SE
  python blast_extract_cds.py \\
      --target assembly.fasta \\
      --ref-gb reference.gb \\
      --out annotated_genes.fasta \\
      --trnascan /path/to/tRNAscan-SE
"""

import argparse
import os
import re
import subprocess
import sys
import shutil
import tempfile
from collections import defaultdict

# ─── Gene sets ────────────────────────────────────────────────────────────────

PCG_NAMES = {
    "ND1", "ND2", "COXI", "COXII", "ATPase8", "ATPase6",
    "COXIII", "ND3", "ND4L", "ND4", "ND5", "ND6", "Cytb",
}
RRNA_NAMES = {"12S", "16S"}
BLAST_GENE_NAMES = PCG_NAMES | RRNA_NAMES

# D-loop (control region) detection
DLOOP_MIN_SIZE = 100  # Minimum gap (bp) to consider as a candidate D-loop

# Name normalization: lowercase alias → canonical name
# Handles the wide variety of naming conventions across GFF/GenBank files.
GENE_ALIASES = {
    # PCGs
    "nd1": "ND1", "nad1": "ND1",
    "nd2": "ND2", "nad2": "ND2",
    "cox1": "COXI", "coxi": "COXI", "co1": "COXI",
    "cox2": "COXII", "coxii": "COXII", "co2": "COXII",
    "atp8": "ATPase8", "atpase8": "ATPase8",
    "atp6": "ATPase6", "atpase6": "ATPase6",
    "cox3": "COXIII", "coxiii": "COXIII", "co3": "COXIII",
    "nd3": "ND3", "nad3": "ND3",
    "nd4l": "ND4L", "nad4l": "ND4L",
    "nd4": "ND4", "nad4": "ND4",
    "nd5": "ND5", "nad5": "ND5",
    "nd6": "ND6", "nad6": "ND6",
    "cytb": "Cytb", "cob": "Cytb",
    # rRNAs — NCBI uses RNR1/RNR2, MitoAnnotator uses 12S/16S
    "12s": "12S", "12srna": "12S", "rns": "12S",
    "rnr1": "12S", "mtrnr1": "12S", "srrna": "12S",
    "16s": "16S", "16srna": "16S", "rnl": "16S",
    "rnr2": "16S", "mtrnr2": "16S", "lrrna": "16S",
}

# BLASTN e-value cutoffs by gene type
EVALUE = {"PCG": "1e-20", "rRNA": "1e-20"}

TRNASCAN_SEARCH_PATHS = [
    "/home/cheahhl814/.pixi/envs/trnascan-se/bin/tRNAscan-SE",
]


def normalize_gene_name(raw):
    """Map a raw gene name (any case/format) to canonical form, or None."""
    key = raw.strip().lower().replace(" ", "").replace("-", "").replace("_", "")
    return GENE_ALIASES.get(key)


def gene_type(name):
    if name == "D-loop":
        return "CR"
    base = name.rsplit("_", 1)[0] if name[-1].isdigit() and "_" in name else name
    if base in PCG_NAMES:
        return "PCG"
    if base in RRNA_NAMES:
        return "rRNA"
    return "tRNA"


# ─── Sequence utilities ──────────────────────────────────────────────────────

def read_first_seq(path):
    """Read the first (only) sequence from a single-record FASTA.
    Returns (seqid, sequence).
    """
    seqid = None
    parts = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                seqid = line[1:].split()[0]
            else:
                parts.append(line.strip())
    return seqid or "seq", "".join(parts)


def read_fasta(path):
    """Read a multi-sequence FASTA into {name: sequence} dict."""
    seqs = {}
    name = None
    parts = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(parts)
                name = line[1:].split()[0]
                parts = []
            elif name is not None:
                parts.append(line)
    if name is not None:
        seqs[name] = "".join(parts)
    return seqs


def reverse_complement(seq):
    comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(comp)[::-1]


# ─── Reference gene extraction: GFF3 ────────────────────────────────────────

def _parse_gff_attributes(attr_str):
    """Parse GFF3 attribute column into a dict."""
    attrs = {}
    for item in attr_str.split(";"):
        if "=" not in item:
            continue
        key, val = item.split("=", 1)
        attrs[key.strip()] = val.strip()
    return attrs


def _extract_fasta_from_gff(gff_path):
    """Read sequences from an embedded ##FASTA section in a GFF3 file."""
    seqs = {}
    name = None
    parts = []
    in_fasta = False
    with open(gff_path) as f:
        for line in f:
            if line.startswith("##FASTA"):
                in_fasta = True
                continue
            if in_fasta:
                if line.startswith(">"):
                    if name is not None:
                        seqs[name] = "".join(parts)
                    name = line[1:].split()[0]
                    parts = []
                else:
                    parts.append(line.rstrip("\n"))
    if name is not None:
        seqs[name] = "".join(parts)
    return seqs


def read_ref_genes_from_gff(gff_path, fasta_path=None):
    """
    Extract reference gene sequences from GFF3 annotations.

    Requires either a separate --ref-fasta or a ##FASTA section in the GFF3.
    Prefers CDS/rRNA features over generic 'gene' features when both exist.
    """
    # Load reference sequences
    if fasta_path:
        ref_seqs = read_fasta(fasta_path)
    else:
        ref_seqs = _extract_fasta_from_gff(gff_path)
        if not ref_seqs:
            print("  Error: no --ref-fasta provided and no ##FASTA section "
                  "found in GFF3 file", file=sys.stderr)
            return {}

    # Determine genome size (for wrap-around detection)
    genome_size = 0
    for seqid, seq in ref_seqs.items():
        genome_size = max(genome_size, len(seq))

    # Parse GFF3 annotations
    gene_regions = {}  # canonical_name → (seqid, start, end, strand, feat_type)
    with open(gff_path) as f:
        for line in f:
            if line.startswith("##sequence-region"):
                parts = line.strip().split()
                if len(parts) >= 4:
                    genome_size = max(genome_size, int(parts[3]))
                continue
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue

            seqid  = parts[0]
            ftype  = parts[2]
            start  = int(parts[3])
            end    = int(parts[4])
            strand = parts[6] if parts[6] != "." else "+"
            attrs  = _parse_gff_attributes(parts[8])

            # Accept gene, CDS, rRNA, tRNA, and D_loop features
            if ftype not in ("gene", "CDS", "rRNA", "tRNA", "D_loop"):
                continue

            # Handle D_loop feature type (NCBI encodes wrap-around as end > genome_size)
            if ftype == "D_loop":
                if genome_size > 0 and end > genome_size:
                    end = end - genome_size  # now end < start → wrap-around
                gene_name = "D-loop"
            else:
                # Try to extract gene name from attributes
                gene_name = None
                for key in ("gene", "Name", "label", "product"):
                    if key in attrs:
                        gene_name = normalize_gene_name(attrs[key])
                        if gene_name:
                            break

            if gene_name is None or gene_name not in BLAST_GENE_NAMES:
                continue

            # Prefer CDS over gene, rRNA over gene for the same gene name
            existing = gene_regions.get(gene_name)
            if existing is not None:
                _, _, _, _, existing_type = existing
                if ftype == "gene" and existing_type in ("CDS", "rRNA"):
                    continue  # keep the more specific feature

            gene_regions[gene_name] = (seqid, start, end, strand, ftype)

    # Extract sequences from reference
    result = {}
    for gene_name, (seqid, start, end, strand, _) in gene_regions.items():
        if seqid not in ref_seqs:
            print(f"  Warning: seqid '{seqid}' not found in reference FASTA",
                  file=sys.stderr)
            continue
        seq = ref_seqs[seqid][start - 1 : end]  # GFF3 is 1-based inclusive
        if strand == "-":
            seq = reverse_complement(seq)
        result[gene_name] = seq

    missing = BLAST_GENE_NAMES - set(result)
    if missing:
        print(f"  Warning: reference missing: {', '.join(sorted(missing))}",
              file=sys.stderr)
    return result


# ─── Reference gene extraction: GenBank ──────────────────────────────────────

def read_ref_genes_from_gb(gb_path):
    """
    Extract reference gene sequences from a GenBank (.gb/.gbk) file.

    Parses both the FEATURES section (for coordinates) and the ORIGIN section
    (for sequence data).  Prefers CDS/rRNA features over generic 'gene' features.
    """
    with open(gb_path) as f:
        content = f.read()

    # ── Parse ORIGIN section for the full sequence ──────────────────────────
    origin_idx = content.find("ORIGIN")
    if origin_idx == -1:
        print("  Error: no ORIGIN section found in GenBank file",
              file=sys.stderr)
        return {}

    seq_text = content[origin_idx + 6:]
    full_seq = re.sub(r'[\d\s]', '', seq_text)
    full_seq = re.sub(r'[^ACGTNacgtn]', '', full_seq).upper()

    # ── Parse FEATURES section ──────────────────────────────────────────────
    feat_idx = content.find("FEATURES")
    if feat_idx == -1:
        print("  Error: no FEATURES section found in GenBank file",
              file=sys.stderr)
        return {}

    features_text = content[feat_idx:origin_idx]

    # Split into feature blocks (each starts with 5 spaces + a keyword)
    blocks = re.split(r'\n(?=     \S)', features_text)

    gene_regions = {}  # canonical_name → (start, end, strand, feat_type)

    for block in blocks:
        block = block.strip()
        if not block or block.startswith("Location"):
            continue

        lines = block.split("\n")
        # First line: "keyword   location"
        first = re.match(r'(\S+)\s+(.+)', lines[0])
        if not first:
            continue

        feat_type = first.group(1)
        location_str = first.group(2).strip()

        # Handle multiline location (rare but possible)
        i = 1
        while i < len(lines) and not lines[i].strip().startswith("/"):
            location_str += lines[i].strip()
            i += 1

        if feat_type not in ("gene", "CDS", "rRNA", "tRNA"):
            continue

        # Find gene name from qualifiers
        block_text = "\n".join(lines)
        gene_name = None
        for pattern in [r'/gene="([^"]+)"', r'/Name="([^"]+)"',
                        r'/label="([^"]+)"']:
            match = re.search(pattern, block_text)
            if match:
                canonical = normalize_gene_name(match.group(1))
                if canonical:
                    gene_name = canonical
                    break

        if gene_name is None or gene_name not in BLAST_GENE_NAMES:
            continue

        # Parse location
        complement = "complement" in location_str
        strand = "-" if complement else "+"

        loc_clean = location_str
        if complement:
            loc_clean = re.sub(r'^complement\(', '', loc_clean)
            loc_clean = re.sub(r'\)$', '', loc_clean)

        # Handle join() — may be multi-exon or wrap-around (D-loop)
        if "join" in loc_clean:
            loc_clean = re.sub(r'^join\(', '', loc_clean)
            loc_clean = re.sub(r'\)$', '', loc_clean)
            segments = re.findall(r'[<?]?(\d+)\.\.[>?]?(\d+)', loc_clean)
            if not segments:
                continue
            if len(segments) > 1 and int(segments[-1][1]) < int(segments[0][0]):
                # Wrap-around feature (e.g. D-loop spanning the origin)
                start = int(segments[0][0])
                end = int(segments[-1][1])   # end < start signals wrap-around
            else:
                start = min(int(s[0]) for s in segments)
                end = max(int(s[1]) for s in segments)
        else:
            # Strip < and > (partial boundaries, e.g. <447..551)
            range_match = re.search(r'[<]?(\d+)\.\.[>]?(\d+)', loc_clean)
            if not range_match:
                continue
            start = int(range_match.group(1))
            end = int(range_match.group(2))

        # Prefer CDS over gene, rRNA over gene
        existing = gene_regions.get(gene_name)
        if existing is not None:
            _, _, _, existing_type = existing
            if feat_type == "gene" and existing_type in ("CDS", "rRNA"):
                continue

        gene_regions[gene_name] = (start, end, strand, feat_type)

    # ── Extract sequences ───────────────────────────────────────────────────
    result = {}
    for gene_name, (start, end, strand, _) in gene_regions.items():
        seq = full_seq[start - 1 : end]  # GB is 1-based inclusive
        if strand == "-":
            seq = reverse_complement(seq)
        result[gene_name] = seq

    missing = BLAST_GENE_NAMES - set(result)
    if missing:
        print(f"  Warning: reference missing: {', '.join(sorted(missing))}",
              file=sys.stderr)
    return result


# ─── GeneRecord ──────────────────────────────────────────────────────────────

class GeneRecord:
    __slots__ = ("name", "start", "end", "strand", "sequence")

    def __init__(self, name, start, end, strand, sequence):
        self.name     = name
        self.start    = start
        self.end      = end
        self.strand   = strand
        self.sequence = sequence


# ─── BLAST extraction ─────────────────────────────────────────────────────────

def blast_extract_genes(target_seq, ref_seqs, blast_db, tmpdir):
    """Extract PCGs and rRNAs from the target genome via BLASTN."""
    records = []
    for gene, ref_seq in sorted(ref_seqs.items()):
        gtype   = "rRNA" if gene in RRNA_NAMES else "PCG"
        evalue  = EVALUE[gtype]
        qfile   = os.path.join(tmpdir, f"{gene}.fa")
        with open(qfile, "w") as f:
            f.write(f">{gene}\n{ref_seq}\n")

        result = subprocess.run(
            ["blastn", "-query", qfile, "-db", blast_db,
             "-outfmt", "6 sseqid sstart send pident length qcovs bitscore",
             "-evalue", evalue, "-max_target_seqs", "1"],
            capture_output=True, text=True,
        )
        hits = [h for h in result.stdout.strip().split("\n") if h]

        if not hits:
            print(f"  {gene:10s}  WARNING: no BLAST hit")
            continue

        fields   = hits[0].split("\t")
        sstart   = int(fields[1])
        send     = int(fields[2])
        rev      = sstart > send
        lo, hi   = (send, sstart) if rev else (sstart, send)
        strand   = "-" if rev else "+"

        seq = target_seq[lo - 1 : hi]
        if rev:
            seq = reverse_complement(seq)

        # Trim PCG to first in-frame ATG (within 9 bp)
        if gtype == "PCG":
            atg = seq.lower().find("atg")
            if 0 <= atg < 9:
                seq = seq[atg:]

        diff = len(seq) - len(ref_seq)
        flag = "OK" if abs(diff) <= 9 else f"{diff:+d} bp"
        print(f"  {gene:10s}  {gtype:4s}  {lo:7d}–{hi:7d} {strand}  "
              f"{len(seq):5d} bp  [{flag}]")
        records.append(GeneRecord(gene, lo, hi, strand, seq))

    return records


# ─── tRNAscan-SE ─────────────────────────────────────────────────────────────

def find_trnascan(user_path=None):
    if user_path:
        return user_path
    for p in TRNASCAN_SEARCH_PATHS:
        if os.path.exists(p):
            return p
    return shutil.which("tRNAscan-SE")


def run_trnascan(target_fasta, target_seq, tmpdir, trnascan_bin):
    """
    Run tRNAscan-SE in vertebrate mitochondrial mode (-M vert).
    Parses the tab output for coordinates and amino acid types.
    Extracts sequences directly from the genome (robust against header format
    variations between tRNAscan-SE versions).
    """
    tab_out = os.path.join(tmpdir, "trnascan.txt")

    result = subprocess.run(
        [trnascan_bin, "-M", "vert", "--brief",
         "-o", tab_out,
         target_fasta],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        print(f"  tRNAscan-SE failed:\n{result.stderr[-600:]}", file=sys.stderr)
        return []

    if not os.path.exists(tab_out):
        print("  tRNAscan-SE produced no output file", file=sys.stderr)
        return []

    records = []
    with open(tab_out) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                continue
            try:
                start   = int(parts[2])
                end     = int(parts[3])
            except ValueError:
                continue

            aa_type = parts[4].strip()  # e.g. "Phe", "Val", "Leu"
            if aa_type in ("Undet", "Sup", "Pseudo", ""):
                continue

            strand = "-" if start > end else "+"
            lo, hi = (end, start) if start > end else (start, end)

            seq = target_seq[lo - 1 : hi]
            if strand == "-":
                seq = reverse_complement(seq)

            name = f"tRNA-{aa_type}"
            records.append(GeneRecord(name, lo, hi, strand, seq))

    # Sort by position so duplicates are numbered in genomic order
    records.sort(key=lambda r: r.start)
    print(f"  {len(records)} tRNAs annotated")
    return records


# ─── D-loop (control region) detection ─────────────────────────────────────

def detect_dloop(records, target_seq):
    """
    Detect the D-loop (control region) as the largest intergenic gap.

    In vertebrate mitogenomes the D-loop is almost always the largest
    non-coding region, typically located between tRNA-Pro and tRNA-Phe.

    For circular genomes the D-loop may span the origin; this is handled
    by also checking the wrap-around gap (end of last gene → start of
    first gene via the sequence boundary).

    Returns a GeneRecord or None if no sufficiently large gap is found.
    """
    if not records:
        return None

    seq_len = len(target_seq)
    sorted_recs = sorted(records, key=lambda r: r.start)

    gaps = []

    # Internal gaps between consecutive genes
    for i in range(len(sorted_recs) - 1):
        gap_start = sorted_recs[i].end + 1
        gap_end   = sorted_recs[i + 1].start - 1
        gap_len   = gap_end - gap_start + 1
        if gap_len >= DLOOP_MIN_SIZE:
            gaps.append({
                'len': gap_len, 'start': gap_start, 'end': gap_end,
                'wrap': False,
            })

    # Wrap-around gap (circular genome): last gene → first gene
    # via the end of the sequence
    tail_len = seq_len - sorted_recs[-1].end
    head_len = sorted_recs[0].start - 1
    wrap_len = tail_len + head_len
    if wrap_len >= DLOOP_MIN_SIZE:
        gaps.append({
            'len': wrap_len,
            'start': sorted_recs[-1].end + 1,
            'end': sorted_recs[0].start - 1,
            'wrap': True,
        })

    if not gaps:
        return None

    best = max(gaps, key=lambda g: g['len'])

    if best['wrap']:
        # Wrap-around: D-loop spans end → start across the origin
        seq = target_seq[best['start'] - 1:] + target_seq[:best['end']]
    else:
        seq = target_seq[best['start'] - 1 : best['end']]

    return GeneRecord("D-loop", best['start'], best['end'], "+", seq)


# ─── Output ──────────────────────────────────────────────────────────────────

def disambiguate(records):
    """
    Append _1, _2 … to names that appear more than once (e.g. tRNA-Leu x2).
    Applied after sorting by position so numbering reflects genomic order.
    """
    counts = defaultdict(int)
    for r in records:
        counts[r.name] += 1
    seen = defaultdict(int)
    for r in records:
        if counts[r.name] > 1:
            seen[r.name] += 1
            r.name = f"{r.name}_{seen[r.name]}"
    return records


def write_output(records, out_path):
    """Sort by genomic position, disambiguate names, write FASTA."""
    records.sort(key=lambda r: r.start)
    disambiguate(records)
    with open(out_path, "w") as f:
        for r in records:
            if r.sequence:
                f.write(f">{r.name}\n{r.sequence}\n")
    return records


def write_gff3(records, out_path, seqid, genome_size):
    """Write GFF3 annotation file."""
    type_map = {"PCG": "CDS", "rRNA": "rRNA", "tRNA": "tRNA", "CR": "D_loop"}
    with open(out_path, "w") as f:
        f.write("##gff-version 3\n")
        f.write(f"##sequence-region {seqid} 1 {genome_size}\n")
        for r in records:
            if not r.sequence:
                continue
            gtype = gene_type(r.name)
            ftype = type_map.get(gtype, "gene")
            # GFF3 is 1-based, inclusive — same as our coordinates
            attrs = f"ID={r.name};Name={r.name};gene_biotype={gtype}"
            if gtype == "PCG":
                attrs += ";gbkey=CDS"
            elif gtype == "rRNA":
                attrs += ";gbkey=rRNA"
            elif gtype == "tRNA":
                attrs += ";gbkey=tRNA"
            elif gtype == "CR":
                attrs += ";gbkey=D-loop"
            f.write(f"{seqid}\tblast_extract_cds\t{ftype}\t{r.start}\t{r.end}\t."
                    f"\t{r.strand}\t.\t{attrs}\n")


def write_bed(records, out_path, seqid):
    """Write BED (6-column) coordinate file."""
    with open(out_path, "w") as f:
        f.write(f"track name=mitogenome_genes description=\"Annotated mitogenome genes\"\n")
        for r in records:
            if not r.sequence:
                continue
            # BED is 0-based, half-open: start-1, end (same as end in 1-based inclusive)
            bed_start = r.start - 1
            bed_end = r.end if r.end >= r.start else r.end  # wrap-around: keep as-is
            gtype = gene_type(r.name)
            score = len(r.sequence)
            f.write(f"{seqid}\t{bed_start}\t{bed_end}\t{r.name}\t{score}\t{r.strand}\n")


# ─── Main ────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ref_group = ap.add_mutually_exclusive_group(required=True)
    ref_group.add_argument("--ref-gff",
                          help="Reference GFF3 annotation file")
    ref_group.add_argument("--ref-gb",
                          help="Reference GenBank (.gb/.gbk) file")

    ap.add_argument("--ref-fasta", default=None,
                    help="Reference genome FASTA (required with --ref-gff unless "
                         "GFF3 has a ##FASTA section)")
    ap.add_argument("--target",    required=True,
                    help="Target mitogenome FASTA (single sequence)")
    ap.add_argument("--out",       required=True,
                    help="Output base path (extensions .fasta/.gff3/.bed auto-appended)")
    ap.add_argument("--trnascan",  default=None,
                    help="Path to tRNAscan-SE binary (default: auto-detect)")
    ap.add_argument("--no-trna",   action="store_true",
                    help="Skip tRNA annotation (output PCGs + rRNAs only)")
    args = ap.parse_args()

    # Validate reference arguments
    if args.ref_gff and not args.ref_fasta:
        # Check for embedded FASTA
        has_embedded = False
        with open(args.ref_gff) as f:
            for line in f:
                if line.startswith("##FASTA"):
                    has_embedded = True
                    break
                if not line.startswith("#") and line.strip():
                    break  # past header, no FASTA ahead
        if not has_embedded:
            ap.error("--ref-fasta is required when using --ref-gff "
                     "(unless GFF3 has a ##FASTA section)")

    trnascan_bin = None
    if not args.no_trna:
        trnascan_bin = find_trnascan(args.trnascan)
        if not trnascan_bin:
            print("  Warning: tRNAscan-SE not found — tRNA annotation skipped.\n"
                  "  Install via pixi or pass --trnascan PATH.", file=sys.stderr)

    # Derive output paths
    base = os.path.splitext(args.out)[0]  # strip .fasta/.fa extension
    fasta_path = base + ".fasta" if not args.out.endswith(".fasta") else args.out
    # Use args.out directly for FASTA; derive GFF3/BED from same base
    if args.out.endswith(".fasta") or args.out.endswith(".fa") or args.out.endswith(".fna"):
        base = args.out.rsplit(".", 1)[0]
    else:
        base = args.out
    fasta_path = args.out

    print(f"Target   : {args.target}")
    target_id, target_seq = read_first_seq(args.target)
    print(f"           {len(target_seq):,} bp")

    # Load reference gene sequences
    if args.ref_gff:
        print(f"Reference: {args.ref_gff}"
              + (f" + {args.ref_fasta}" if args.ref_fasta else ""))
        ref_seqs = read_ref_genes_from_gff(args.ref_gff, args.ref_fasta)
    else:
        print(f"Reference: {args.ref_gb}")
        ref_seqs = read_ref_genes_from_gb(args.ref_gb)

    print(f"           {len(ref_seqs)} reference genes loaded")

    tmpdir = tempfile.mkdtemp(prefix="mito_annotate_")
    try:
        db = os.path.join(tmpdir, "blastdb")
        subprocess.run(
            ["makeblastdb", "-in", args.target, "-dbtype", "nucl", "-out", db],
            check=True, capture_output=True,
        )

        print("\nBLAST extraction (PCGs + rRNAs):")
        records = blast_extract_genes(target_seq, ref_seqs, db, tmpdir)

        if trnascan_bin:
            print("\ntRNAscan-SE annotation:")
            records += run_trnascan(args.target, target_seq, tmpdir, trnascan_bin)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    # D-loop (control region) detection — largest intergenic gap
    print("\nD-loop detection:")
    dloop = detect_dloop(records, target_seq)
    if dloop:
        print(f"  D-loop       CR   {dloop.start:>7,d}–{dloop.end:>7,d} {dloop.strand:>3}  "
              f"{len(dloop.sequence):>6} bp  [gap-based]")
        records.append(dloop)
    else:
        print("  No D-loop detected (no gap >= "
              f"{DLOOP_MIN_SIZE} bp)")

    records = write_output(records, fasta_path)

    # Write GFF3 and BED alongside FASTA
    gff_path = base + ".gff3"
    bed_path = base + ".bed"
    write_gff3(records, gff_path, target_id, len(target_seq))
    write_bed(records, bed_path, target_id)

    # Summary table
    print(f"\n{'Gene':<14} {'Type':<5} {'Start':>8} {'End':>8} {'Str':>3} {'Length':>7}")
    print("─" * 52)
    for r in records:
        gtype = gene_type(r.name)
        end_str = f"{r.end:>8,}"
        if r.name == "D-loop" and r.start > r.end:
            # Wrap-around D-loop (spans sequence origin)
            end_str = f"{r.end:>7,}→"
        print(f"  {r.name:<12} {gtype:<5} {r.start:>8,} {end_str} {r.strand:>3}  "
              f"{len(r.sequence):>6} bp")

    n  = len([r for r in records if r.sequence])
    bp = sum(len(r.sequence) for r in records if r.sequence)
    print(f"\n{n} genes  ({bp:,} bp total)")
    print(f"  FASTA → {fasta_path}")
    print(f"  GFF3  → {gff_path}")
    print(f"  BED   → {bed_path}")


if __name__ == "__main__":
    main()