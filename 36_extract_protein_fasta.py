#!/usr/bin/env python3
"""36 - Protein FASTA for all 7 genomes, extracted from the submitted records.

Translations are read from the .gbk/.gbf rendering of the .sqn, so the sequences and
their product names are exactly what is deposited - after the EC resolution, the
NCBI-safe pass and any curated override.  Taking them from an upstream *.proteins.fa
would give pre-correction names instead.

Output: 23_protein_fasta/<Organism>_<strain>.faa
"""
import os, re, sys

P    = "/home/highfieldc/FUNANNOTATE_YEAST_PROJECT"
SUB  = f"{P}/22_ncbi_submission"
OUT  = f"{P}/23_protein_fasta"
REG  = f"{SUB}/SAMPLE_REGISTRY.tsv"

def log(m): sys.stderr.write(m + "\n"); sys.stderr.flush()

def registry():
    rows, hdr = [], []
    for l in open(REG, errors="replace"):
        if not l.strip(): continue
        c = l.rstrip("\n").split("\t")
        if l.startswith("#"): hdr = [x.lstrip("#") for x in c]; continue
        rows.append({hdr[i]: (c[i] if i < len(c) else "") for i in range(len(hdr))})
    return rows

def flatfile(d):
    """The .gbk/.gbf for a sample directory."""
    for f in sorted(os.listdir(d)):
        if f.endswith(".gbk") or f.endswith(".gbf"):
            return os.path.join(d, f)
    return None

def cds_records(path):
    """Yield {locus_tag, protein_id, product, translation} for every CDS.

    Uses Biopython's GenBank parser rather than reading the flatfile by hand.  A
    hand-rolled parser got this wrong: where a CDS is the last feature before ORIGIN
    it kept appending and welded 'ORIGIN' plus the nucleotide sequence onto the end of
    the translation, corrupting 2,004 of 27,480 sequences.  The format has enough
    edge cases - wrapped qualifiers, quotes inside values, features at the block
    boundary - that a maintained parser is the safer choice.
    """
    from Bio import SeqIO
    for rec in SeqIO.parse(path, "genbank"):
        for f in rec.features:
            if f.type != "CDS":
                continue
            tr = f.qualifiers.get("translation", [""])[0]
            if not tr:
                continue                      # pseudogenes carry no translation
            yield {"locus_tag":  f.qualifiers.get("locus_tag",  [""])[0],
                   "protein_id": f.qualifiers.get("protein_id", [""])[0],
                   "product":    f.qualifiers.get("product",    [""])[0],
                   "translation": tr}

def main():
    os.makedirs(OUT, exist_ok=True)
    total = 0
    index = []
    for row in registry():
        d = f"{SUB}/{row['id']}"
        gb = flatfile(d)
        if not gb:
            log(f"  [{row['id']}] no .gbk/.gbf - skipped"); continue
        base = re.sub(r'_+', '_', re.sub(r'[^A-Za-z0-9]+', '_',
                      f"{row['organism']} {row['strain']}")).strip("_")
        dst = f"{OUT}/{base}.faa"
        n = 0
        with open(dst, "w") as out:
            for r in cds_records(gb):
                seq = re.sub(r'[^A-Za-z]', '', r["translation"])
                if not seq: continue
                name = r["locus_tag"] or r["protein_id"] or f"unnamed_{n}"
                prod = r["product"] or "hypothetical protein"
                out.write(f">{name} {prod} [{row['organism']} {row['strain']}]\n")
                for i in range(0, len(seq), 60):
                    out.write(seq[i:i+60] + "\n")
                n += 1
        total += n
        index.append((row["id"], row["organism"], row["strain"], base, n,
                      os.path.getsize(dst)))
        log(f"  [{row['id']}] {n:,} proteins -> {base}.faa")

    with open(f"{OUT}/README.txt", "w") as fh:
        fh.write("Protein sequences for the 7 submitted genomes\n")
        fh.write("=" * 46 + "\n\n")
        fh.write("Extracted from the /translation qualifiers of each genome's .gbk/.gbf,\n"
                 "which is the GenBank rendering of the .sqn actually submitted.  The\n"
                 "sequences and the product names on the header lines are therefore what\n"
                 "is deposited, after EC resolution, the NCBI-safe pass and any curated\n"
                 "override.  funannotate's own *.proteins.fa carries pre-correction names\n"
                 "and is NOT the source used here.\n\n")
        fh.write("Header format:  >locus_tag product [Organism strain]\n\n")
        fh.write(f"{'sample':<7}{'organism':<32}{'strain':<10}{'proteins':>9}  file\n")
        for sid, org, strain, base, n, sz in index:
            fh.write(f"{sid:<7}{org:<32}{strain:<10}{n:>9,}  {base}.faa\n")
        fh.write(f"\n{'total':<49}{total:>9,}\n")
    log(f"\n  {total:,} proteins across {len(index)} genomes -> {OUT}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
