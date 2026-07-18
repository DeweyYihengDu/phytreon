"""Download a large public 16S rRNA set from NCBI (Entrez) and align it.

A broad sampler of the prokaryotic tree: ~90 type strains and model organisms
spanning most named bacterial phyla plus a range of archaea. Bigger and more
taxonomically even than the 18-taxon set in ``fetch_example_data.py``, so it
makes a substantial tree for the tanglegram demo and the gallery.

All sequences are public RefSeq/GenBank 16S rRNA. Writes:

    examples/data/big16S.fasta          unaligned 16S sequences
    examples/data/big16S_aligned.fasta  aligned + gap-trimmed (built-in MSA)
    examples/data/big16S_metadata.csv   per-tip metadata (domain, phylum, ...)

The alignment step is pure Python and quadratic in the number of sequences --
budget ~10-20 minutes for the full set. It runs once; the outputs are
committed so the demo itself starts from the aligned file.

Run:  python examples/data/fetch_large_16S.py
"""
import csv
import os
import sys
import time

from Bio import Entrez, SeqIO

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))))
import phytreon as pt                                          # noqa: E402

Entrez.email = "yiheng.du@slu.edu"          # required by NCBI
HERE = os.path.dirname(os.path.abspath(__file__))

# organism -> (tip label, domain, phylum). Phylum names follow the current
# (post-2021) nomenclature, matching tol_metadata.csv.
ORGANISMS = [
    # -- Pseudomonadota -------------------------------------------------
    ("Escherichia coli", "Escherichia_coli", "Bacteria", "Pseudomonadota"),
    ("Salmonella enterica", "Salmonella_enterica", "Bacteria", "Pseudomonadota"),
    ("Klebsiella pneumoniae", "Klebsiella_pneumoniae", "Bacteria", "Pseudomonadota"),
    ("Yersinia pestis", "Yersinia_pestis", "Bacteria", "Pseudomonadota"),
    ("Serratia marcescens", "Serratia_marcescens", "Bacteria", "Pseudomonadota"),
    ("Proteus mirabilis", "Proteus_mirabilis", "Bacteria", "Pseudomonadota"),
    ("Vibrio cholerae", "Vibrio_cholerae", "Bacteria", "Pseudomonadota"),
    ("Aeromonas hydrophila", "Aeromonas_hydrophila", "Bacteria", "Pseudomonadota"),
    ("Pseudomonas aeruginosa", "Pseudomonas_aeruginosa", "Bacteria", "Pseudomonadota"),
    ("Pseudomonas putida", "Pseudomonas_putida", "Bacteria", "Pseudomonadota"),
    ("Acinetobacter baumannii", "Acinetobacter_baumannii", "Bacteria", "Pseudomonadota"),
    ("Haemophilus influenzae", "Haemophilus_influenzae", "Bacteria", "Pseudomonadota"),
    ("Neisseria meningitidis", "Neisseria_meningitidis", "Bacteria", "Pseudomonadota"),
    ("Bordetella pertussis", "Bordetella_pertussis", "Bacteria", "Pseudomonadota"),
    ("Burkholderia cepacia", "Burkholderia_cepacia", "Bacteria", "Pseudomonadota"),
    ("Ralstonia solanacearum", "Ralstonia_solanacearum", "Bacteria", "Pseudomonadota"),
    ("Legionella pneumophila", "Legionella_pneumophila", "Bacteria", "Pseudomonadota"),
    ("Coxiella burnetii", "Coxiella_burnetii", "Bacteria", "Pseudomonadota"),
    ("Francisella tularensis", "Francisella_tularensis", "Bacteria", "Pseudomonadota"),
    ("Agrobacterium tumefaciens", "Agrobacterium_tumefaciens", "Bacteria", "Pseudomonadota"),
    ("Bradyrhizobium japonicum", "Bradyrhizobium_japonicum", "Bacteria", "Pseudomonadota"),
    ("Caulobacter vibrioides", "Caulobacter_vibrioides", "Bacteria", "Pseudomonadota"),
    ("Rhodobacter sphaeroides", "Rhodobacter_sphaeroides", "Bacteria", "Pseudomonadota"),
    ("Brucella melitensis", "Brucella_melitensis", "Bacteria", "Pseudomonadota"),
    ("Rickettsia prowazekii", "Rickettsia_prowazekii", "Bacteria", "Pseudomonadota"),
    ("Desulfovibrio vulgaris", "Desulfovibrio_vulgaris", "Bacteria", "Thermodesulfobacteriota"),
    ("Geobacter sulfurreducens", "Geobacter_sulfurreducens", "Bacteria", "Thermodesulfobacteriota"),
    ("Myxococcus xanthus", "Myxococcus_xanthus", "Bacteria", "Myxococcota"),
    ("Bdellovibrio bacteriovorus", "Bdellovibrio_bacteriovorus", "Bacteria", "Bdellovibrionota"),
    # -- Campylobacterota -----------------------------------------------
    ("Helicobacter pylori", "Helicobacter_pylori", "Bacteria", "Campylobacterota"),
    ("Campylobacter jejuni", "Campylobacter_jejuni", "Bacteria", "Campylobacterota"),
    # -- Bacillota ------------------------------------------------------
    ("Bacillus subtilis", "Bacillus_subtilis", "Bacteria", "Bacillota"),
    ("Bacillus cereus", "Bacillus_cereus", "Bacteria", "Bacillota"),
    ("Staphylococcus aureus", "Staphylococcus_aureus", "Bacteria", "Bacillota"),
    ("Staphylococcus epidermidis", "Staphylococcus_epidermidis", "Bacteria", "Bacillota"),
    ("Streptococcus pneumoniae", "Streptococcus_pneumoniae", "Bacteria", "Bacillota"),
    ("Streptococcus pyogenes", "Streptococcus_pyogenes", "Bacteria", "Bacillota"),
    ("Enterococcus faecalis", "Enterococcus_faecalis", "Bacteria", "Bacillota"),
    ("Lactococcus lactis", "Lactococcus_lactis", "Bacteria", "Bacillota"),
    ("Listeria monocytogenes", "Listeria_monocytogenes", "Bacteria", "Bacillota"),
    ("Clostridium perfringens", "Clostridium_perfringens", "Bacteria", "Bacillota"),
    ("Clostridioides difficile", "Clostridioides_difficile", "Bacteria", "Bacillota"),
    ("Veillonella parvula", "Veillonella_parvula", "Bacteria", "Bacillota"),
    # -- Mycoplasmatota -------------------------------------------------
    ("Mycoplasma genitalium", "Mycoplasma_genitalium", "Bacteria", "Mycoplasmatota"),
    ("Mycoplasma pneumoniae", "Mycoplasma_pneumoniae", "Bacteria", "Mycoplasmatota"),
    # -- Actinomycetota -------------------------------------------------
    ("Mycobacterium tuberculosis", "Mycobacterium_tuberculosis", "Bacteria", "Actinomycetota"),
    ("Mycobacterium leprae", "Mycobacterium_leprae", "Bacteria", "Actinomycetota"),
    ("Bifidobacterium longum", "Bifidobacterium_longum", "Bacteria", "Actinomycetota"),
    ("Corynebacterium glutamicum", "Corynebacterium_glutamicum", "Bacteria", "Actinomycetota"),
    ("Corynebacterium diphtheriae", "Corynebacterium_diphtheriae", "Bacteria", "Actinomycetota"),
    ("Streptomyces coelicolor", "Streptomyces_coelicolor", "Bacteria", "Actinomycetota"),
    ("Streptomyces griseus", "Streptomyces_griseus", "Bacteria", "Actinomycetota"),
    ("Micrococcus luteus", "Micrococcus_luteus", "Bacteria", "Actinomycetota"),
    ("Nocardia asteroides", "Nocardia_asteroides", "Bacteria", "Actinomycetota"),
    ("Frankia alni", "Frankia_alni", "Bacteria", "Actinomycetota"),
    # -- Bacteroidota ---------------------------------------------------
    ("Bacteroides fragilis", "Bacteroides_fragilis", "Bacteria", "Bacteroidota"),
    ("Bacteroides thetaiotaomicron", "Bacteroides_thetaiotaomicron", "Bacteria", "Bacteroidota"),
    ("Porphyromonas gingivalis", "Porphyromonas_gingivalis", "Bacteria", "Bacteroidota"),
    ("Prevotella melaninogenica", "Prevotella_melaninogenica", "Bacteria", "Bacteroidota"),
    ("Flavobacterium johnsoniae", "Flavobacterium_johnsoniae", "Bacteria", "Bacteroidota"),
    ("Cytophaga hutchinsonii", "Cytophaga_hutchinsonii", "Bacteria", "Bacteroidota"),
    # -- Chlorobiota / Chloroflexota ------------------------------------
    ("Chlorobaculum tepidum", "Chlorobaculum_tepidum", "Bacteria", "Chlorobiota"),
    ("Chloroflexus aurantiacus", "Chloroflexus_aurantiacus", "Bacteria", "Chloroflexota"),
    ("Dehalococcoides mccartyi", "Dehalococcoides_mccartyi", "Bacteria", "Chloroflexota"),
    # -- Cyanobacteriota ------------------------------------------------
    ("Synechocystis sp. PCC 6803", "Synechocystis_PCC6803", "Bacteria", "Cyanobacteriota"),
    ("Synechococcus elongatus", "Synechococcus_elongatus", "Bacteria", "Cyanobacteriota"),
    ("Prochlorococcus marinus", "Prochlorococcus_marinus", "Bacteria", "Cyanobacteriota"),
    ("Nostoc punctiforme", "Nostoc_punctiforme", "Bacteria", "Cyanobacteriota"),
    ("Microcystis aeruginosa", "Microcystis_aeruginosa", "Bacteria", "Cyanobacteriota"),
    ("Gloeobacter violaceus", "Gloeobacter_violaceus", "Bacteria", "Cyanobacteriota"),
    ("Trichodesmium erythraeum", "Trichodesmium_erythraeum", "Bacteria", "Cyanobacteriota"),
    # -- Spirochaetota --------------------------------------------------
    ("Borreliella burgdorferi", "Borreliella_burgdorferi", "Bacteria", "Spirochaetota"),
    ("Treponema pallidum", "Treponema_pallidum", "Bacteria", "Spirochaetota"),
    ("Leptospira interrogans", "Leptospira_interrogans", "Bacteria", "Spirochaetota"),
    # -- PVC superphylum ------------------------------------------------
    ("Chlamydia trachomatis", "Chlamydia_trachomatis", "Bacteria", "Chlamydiota"),
    ("Chlamydia pneumoniae", "Chlamydia_pneumoniae", "Bacteria", "Chlamydiota"),
    ("Parachlamydia acanthamoebae", "Parachlamydia_acanthamoebae", "Bacteria", "Chlamydiota"),
    ("Rhodopirellula baltica", "Rhodopirellula_baltica", "Bacteria", "Planctomycetota"),
    ("Planctopirus limnophila", "Planctopirus_limnophila", "Bacteria", "Planctomycetota"),
    ("Gemmata obscuriglobus", "Gemmata_obscuriglobus", "Bacteria", "Planctomycetota"),
    ("Akkermansia muciniphila", "Akkermansia_muciniphila", "Bacteria", "Verrucomicrobiota"),
    ("Verrucomicrobium spinosum", "Verrucomicrobium_spinosum", "Bacteria", "Verrucomicrobiota"),
    ("Opitutus terrae", "Opitutus_terrae", "Bacteria", "Verrucomicrobiota"),
    # -- deep-branching bacteria ----------------------------------------
    ("Deinococcus radiodurans", "Deinococcus_radiodurans", "Bacteria", "Deinococcota"),
    ("Thermus thermophilus", "Thermus_thermophilus", "Bacteria", "Deinococcota"),
    ("Aquifex aeolicus", "Aquifex_aeolicus", "Bacteria", "Aquificota"),
    ("Hydrogenobacter thermophilus", "Hydrogenobacter_thermophilus", "Bacteria", "Aquificota"),
    ("Thermotoga maritima", "Thermotoga_maritima", "Bacteria", "Thermotogota"),
    ("Fervidobacterium nodosum", "Fervidobacterium_nodosum", "Bacteria", "Thermotogota"),
    ("Fusobacterium nucleatum", "Fusobacterium_nucleatum", "Bacteria", "Fusobacteriota"),
    ("Acidobacterium capsulatum", "Acidobacterium_capsulatum", "Bacteria", "Acidobacteriota"),
    ("Leptospirillum ferrooxidans", "Leptospirillum_ferrooxidans", "Bacteria", "Nitrospirota"),
    # -- Archaea --------------------------------------------------------
    ("Methanocaldococcus jannaschii", "Methanocaldococcus_jannaschii", "Archaea", "Euryarchaeota"),
    ("Methanococcus maripaludis", "Methanococcus_maripaludis", "Archaea", "Euryarchaeota"),
    ("Methanosarcina barkeri", "Methanosarcina_barkeri", "Archaea", "Euryarchaeota"),
    ("Methanobrevibacter smithii", "Methanobrevibacter_smithii", "Archaea", "Euryarchaeota"),
    ("Halobacterium salinarum", "Halobacterium_salinarum", "Archaea", "Euryarchaeota"),
    ("Haloferax volcanii", "Haloferax_volcanii", "Archaea", "Euryarchaeota"),
    ("Archaeoglobus fulgidus", "Archaeoglobus_fulgidus", "Archaea", "Euryarchaeota"),
    ("Thermococcus kodakarensis", "Thermococcus_kodakarensis", "Archaea", "Euryarchaeota"),
    ("Pyrococcus furiosus", "Pyrococcus_furiosus", "Archaea", "Euryarchaeota"),
    ("Thermoplasma acidophilum", "Thermoplasma_acidophilum", "Archaea", "Euryarchaeota"),
    ("Saccharolobus solfataricus", "Saccharolobus_solfataricus", "Archaea", "Thermoproteota"),
    ("Sulfolobus acidocaldarius", "Sulfolobus_acidocaldarius", "Archaea", "Thermoproteota"),
    ("Pyrobaculum aerophilum", "Pyrobaculum_aerophilum", "Archaea", "Thermoproteota"),
    ("Thermoproteus tenax", "Thermoproteus_tenax", "Archaea", "Thermoproteota"),
    ("Nitrosopumilus maritimus", "Nitrosopumilus_maritimus", "Archaea", "Nitrososphaerota"),
]


def fetch_16s(organism: str):
    """esearch a near-full-length 16S, return (accession, sequence) or None."""
    term = (f'{organism}[ORGN] AND 16S ribosomal RNA[Title] '
            f'AND 1200:1600[SLEN] AND biomol_rrna[PROP]')
    h = Entrez.esearch(db="nucleotide", term=term, retmax=3)
    ids = Entrez.read(h)["IdList"]
    h.close()
    if not ids:
        # relax: drop the rRNA biomol filter
        term2 = f'{organism}[ORGN] AND 16S ribosomal RNA[Title] AND 1200:1600[SLEN]'
        h = Entrez.esearch(db="nucleotide", term=term2, retmax=3)
        ids = Entrez.read(h)["IdList"]
        h.close()
    if not ids:
        return None
    h = Entrez.efetch(db="nucleotide", id=ids[0], rettype="fasta", retmode="text")
    rec = SeqIO.read(h, "fasta")
    h.close()
    return rec.id, str(rec.seq)


def write_fasta(path, pairs):
    with open(path, "w") as f:
        for name, seq in pairs:
            f.write(f">{name}\n")
            for i in range(0, len(seq), 70):
                f.write(seq[i:i + 70] + "\n")


def main():
    rows, fasta = [], []
    print(f"fetching {len(ORGANISMS)} 16S sequences from NCBI ...")
    for organism, label, domain, phylum in ORGANISMS:
        try:
            res = fetch_16s(organism)
        except Exception as e:                      # network / parse hiccup
            print(f"  ! {organism}: {e}")
            res = None
        if res is None:
            print(f"  - {organism}: no hit, skipped")
            continue
        acc, seq = res
        fasta.append((label, seq))
        rows.append({"name": label, "organism": organism, "domain": domain,
                     "phylum": phylum, "accession": acc, "length": len(seq)})
        print(f"  + {label:34s} {acc:14s} {len(seq)} bp  ({phylum})")
        time.sleep(0.4)                             # be nice to NCBI

    write_fasta(os.path.join(HERE, "big16S.fasta"), fasta)
    with open(os.path.join(HERE, "big16S_metadata.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["name", "organism", "domain", "phylum",
                                          "accession", "length"])
        w.writeheader()
        w.writerows(rows)
    print(f"\nfetched {len(fasta)} sequences")

    print("aligning (pure Python, quadratic -- this is the slow part) ...")
    t0 = time.time()
    aln = pt.align(fasta)
    aln = pt.trim(aln, max_gap=0.5)
    print(f"  aligned in {time.time() - t0:.0f}s -> "
          f"{aln.nseq} seqs x {aln.ncol} columns")
    write_fasta(os.path.join(HERE, "big16S_aligned.fasta"),
                list(zip(aln.names, aln.seqs)))
    print("wrote big16S.fasta / big16S_aligned.fasta / big16S_metadata.csv")


if __name__ == "__main__":
    main()
