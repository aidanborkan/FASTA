#UniProt provides extremely powerful query capabilities, but the output is typically a **raw FASTA file**. FASTA files are ideal for sequence-based tools, but they are difficult to work with in downstream data analysis, annotation, and integration workflows.

This application fills the gap between **UniProt FASTA datasets** and **proteomics analysis tables**.

Specifically, it is designed for situations where you want to answer questions like:

- *Which proteins from a curated UniProt set are actually present in my experiment?*
- *How many of my significant MS hits belong to a specific functional class (e.g., kinases, metal-binding proteins, reviewed human proteins)?*
- *Which UniProt proteins should I carry forward into structural, enrichment, or targeted analyses?*

---

### What problem does this solve?

In many proteomics workflows, you work with **two disconnected data sources**:

1. **A UniProt-derived FASTA set**  
   Examples:
   - Reviewed human proteome
   - Proteins annotated with a specific keyword (e.g., copper-binding)
   - Proteins from a specific organism or proteome
   - Proteins matching a functional or annotation-based query

2. **An experimental or analytical protein list**  
   Examples:
   - MS-identified proteins
   - Differential abundance hits
   - LiP-MS or TPP significant proteins
   - Curated protein sets from prior analyses

This app:
- Fetches FASTA data *directly* from UniProt using a reproducible query
- Converts semi-structured FASTA headers into a **tidy, joinable table**
- Allows you to **intersect UniProt-derived protein sets with your own data**
- Produces outputs that are immediately usable in downstream analysis

---

### Why not just download FASTA manually?

Manually downloading FASTA files and inspecting headers is:
- Error-prone
- Difficult to reproduce
- Hard to integrate with data frames and pipelines
- Unscalable for iterative analysis

This app:
- Makes UniProt queries **explicit and documented**
- Ensures identifier parsing is **consistent and deterministic**
- Eliminates copy-paste and manual filtering
- Encourages reproducible, query-driven workflows

---

### When is this especially useful?

This tool is particularly helpful when you want to:

- Filter large UniProt protein sets down to **experimentally observed proteins**
- Validate whether expected functional classes appear in your data
- Prepare protein subsets for:
  - Structural modeling
  - Enrichment analysis
  - Targeted follow-up experiments
  - Custom FASTA databases
- Rapidly iterate on different UniProt queries without re-downloading files manually

---

### UniProt query examples

Below are example queries that can be entered directly into the app:

- **Human reviewed proteome**  
  `proteome: UP000005640 AND reviewed: true`

- **Human kinases (by protein name)**  
  `organism_id:9606 AND (protein_name:kinase)`

- **Copper-binding proteins**  
  `keyword: "Copper-binding" AND reviewed: true.`

The application uses UniProt’s REST **stream endpoint** (`uniprotkb/stream`) to retrieve FASTA records directly from UniProt, ensuring results are current and reproducible.
