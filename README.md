# NGS Variant Filtering & Prioritization Pipeline

An automated, rule-based clinical informatics script designed to filter, annotate, and prioritize high-confidence pathogenic variants from Next-Generation Sequencing (NGS) gene panels. This tool simulates a core bioinformatics step in translational oncology and precision medicine workflows.

## Key Features
- **Frequency Filtering**: Dynamically filters out common polymorphisms using population data (`gnomAD_AF` cutoff set to $< 1\%$).
- **Multi-Omics Assessment**: Integrates functional consequence predictions, in silico impact metrics (`CADD Phred` $\ge 20$), and ClinVar clinical significance classifications.
- **ACMG-aligned Prioritization**: Automatically flags variants classified as Pathogenic or Likely Pathogenic, and prioritizes novel or high-impact protein-truncating variants (e.g., stop-gained, frameshift mutations).
- **Clinical Readiness**: Outputs a tab-separated list sorted by clinical urgency for easy import into downstream reporting software or review by clinical geneticists.

---

## Repository Structure
```text
├── data/
│   ├── example_variants.tsv      # Input mock NGS variant dataset
│   └── prioritized_variants.tsv  # Output generated clinical-ready table
├── scripts/
│   └── filter_variants.py        # Core Python processing logic
├── .gitignore
├── LICENSE
└── requirements.txt              # Standard package requirements
