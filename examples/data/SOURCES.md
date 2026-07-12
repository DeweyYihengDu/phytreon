# Example data sources

Small public **16S rRNA** set �� a "microbial tree of life" sampler used
by phytreon's examples and tests. Downloaded from the NCBI `nucleotide`
database (16S ribosomal RNA records, preferring RefSeq `NR_*` type strains)
with `examples/data/fetch_example_data.py` (Biopython Entrez).

**License:** public domain. NCBI sequence records are not subject to
copyright; redistribution as example data is unrestricted.

**Note:** accession versions (`NR_*.x`) may advance over time, so a fresh
download need not byte-match the snapshot shipped here; the cached
`tol_16S_aligned.fasta` keeps the examples reproducible offline.

| tip label | organism | domain | phylum | accession | length |
|---|---|---|---|---|---|
| `Escherichia_coli` | *Escherichia coli* | Bacteria | Pseudomonadota | NR_114042.1 | 1467 |
| `Salmonella_enterica` | *Salmonella enterica* | Bacteria | Pseudomonadota | NR_074910.1 | 1544 |
| `Pseudomonas_aeruginosa` | *Pseudomonas aeruginosa* | Bacteria | Pseudomonadota | NR_117678.1 | 1527 |
| `Vibrio_cholerae` | *Vibrio cholerae* | Bacteria | Pseudomonadota | NR_117894.1 | 1461 |
| `Helicobacter_pylori` | *Helicobacter pylori* | Bacteria | Campylobacterota | NR_044761.1 | 1503 |
| `Bacillus_subtilis` | *Bacillus subtilis* | Bacteria | Bacillota | NR_102783.2 | 1550 |
| `Staphylococcus_aureus` | *Staphylococcus aureus* | Bacteria | Bacillota | NR_037007.2 | 1552 |
| `Streptococcus_pneumoniae` | *Streptococcus pneumoniae* | Bacteria | Bacillota | NR_117496.1 | 1426 |
| `Mycobacterium_tuberculosis` | *Mycobacterium tuberculosis* | Bacteria | Actinomycetota | NR_044826.2 | 1532 |
| `Bifidobacterium_longum` | *Bifidobacterium longum* | Bacteria | Actinomycetota | NR_145535.1 | 1494 |
| `Bacteroides_fragilis` | *Bacteroides fragilis* | Bacteria | Bacteroidota | NR_074784.2 | 1529 |
| `Synechocystis_PCC6803` | *Synechocystis sp. PCC 6803* | Bacteria | Cyanobacteriota | AY224195.1 | 1238 |
| `Thermus_thermophilus` | *Thermus thermophilus* | Bacteria | Deinococcota | NR_113293.1 | 1457 |
| `Aquifex_aeolicus` | *Aquifex aeolicus* | Bacteria | Aquificota | NR_075056.2 | 1584 |
| `Methanocaldococcus_jannaschii` | *Methanocaldococcus jannaschii* | Archaea | Euryarchaeota | NR_113292.1 | 1458 |
| `Halobacterium_salinarum` | *Halobacterium salinarum* | Archaea | Euryarchaeota | NR_025555.1 | 1434 |
| `Saccharolobus_solfataricus` | *Saccharolobus solfataricus* | Archaea | Thermoproteota | NR_119198.1 | 1436 |
| `Pyrococcus_furiosus` | *Pyrococcus furiosus* | Archaea | Euryarchaeota | NR_113294.1 | 1321 |

## Single-cell CRISPR lineage-tracing sample

`lineage_alleletable.txt` (raw allele table) and `lineage_reference_tree.nwk`
(published reconstruction, topology only) are one real sample -- **3432_NT_T1**,
226 cells, 10 CRISPR integration barcodes -- from the **KP-Tracer** mouse lung
adenocarcinoma lineage-tracing study (Yang, Jones, Chan, et al. 2022,
*Lineage tracing reveals the phylodynamics, plasticity, and paths of tumor
evolution*, Cell 185:1905-1923, doi:10.1016/j.cell.2022.04.015), redistributed
by the **Cassiopeia** package (Jones, Khodaverdian, Quinn, et al. 2020,
*Inference of single-cell phylogenies from lineage tracing data using
Cassiopeia*, Genome Biology 21:92; `github.com/YosefLab/Cassiopeia`, MIT
license).

Downloaded with `examples/data/fetch_lineage_data.py`, which also strips the
reference tree's custom ete3 per-node annotation format down to a plain
Newick topology (see that script and `examples/lineage_demo.py` for details).
