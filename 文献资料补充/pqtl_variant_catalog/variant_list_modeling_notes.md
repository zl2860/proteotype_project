# Nature pQTL supplement: variant list and model stratification notes

## Generated files

- `pqtl_variant_catalog.csv`: merged candidate variants and rsIDs extracted from
  pQTL, fine-mapping, PAV, regulatory, colocalization, and trans-network sheets.
- `pqtl_variant_catalog_summary.csv`: counts by source sheet and evidence type.
- `supplement_sheet_catalog.csv`: sheet-level row/column catalog.

The extraction script is `scripts/extract_pqtl_variant_catalog.py`.

## Sheets most useful for variant matrix extraction

| Priority | Sheets | Use |
|---|---|---|
| Core lead pQTL list | ST10, ST9, ST11, ST15 | Direct SNP/protein instruments. ST10 is the best default because it uses the combined cohort. ST9 is discovery-only. ST11 is ancestry-specific. ST15 labels novel vs replicated. |
| Fine-mapped variants | ST16, ST17 | ST16 gives top PIP variants and credible sets. Use top variants for compact models; use credible sets for expanded sensitivity models. ST17 marks credible sets containing PAVs. |
| Functional variants | ST12, ST13, ST14 | ST12 high-impact coding/splicing pQTLs; ST13 regulatory noncoding pQTLs; ST14 lead variants in LD with protein-altering variants. |
| Network/trans biology | ST18, ST20, ST21, ST23, ST29 | pQTL-pQTL colocalization, protein interaction at trans loci, reciprocal trans relationships, receptor-ligand axes, inflammasome examples. |
| Multi-omics colocalization | ST26, ST28, ST30 | ST26 cis pQTL-eQTL colocalization from GTEx. ST28/ST30 are disease-specific examples rather than broad training lists. |
| Covariate robustness | ST24, ST25, ST27 | ST24 tests whether pQTL effects remain after blood-cell/BMI/season/fasting adjustment. ST25/27 are useful for downstream interpretation rather than primary variant extraction. |
| Protein annotation | ST3, ST19 | ST3 maps protein IDs to gene symbols, UniProt, panel, genomic coordinates, QC. ST19 gives SNP heritability components for prioritizing genetically predictable proteins. |

## Recommended variant tiers

1. **Tier A: compact lead pQTL backbone**
   - Source: ST10.
   - Size in extracted catalog: 12,366 unique variant keys.
   - Purpose: first-pass genetically influenced proteotype (GIP) model.
   - Suggested split: cis lead pQTLs vs trans lead pQTLs.

2. **Tier B: discovery/replication labels**
   - Source: ST15 joined to ST9/ST10.
   - Purpose: mark variants as novel or replicated; useful as confidence labels
     or analysis strata.

3. **Tier C: functional high-priority pQTLs**
   - Sources: ST12, ST13, ST14, ST17.
   - Purpose: biologically interpretable models and mechanistic discussion.
   - Subclasses: coding/splicing, regulatory noncoding, PAV-linked, credible-set PAV.

4. **Tier D: fine-mapped expanded set**
   - Source: ST16 credible sets.
   - Size in extracted catalog: 255,116 unique variants.
   - Purpose: sensitivity analysis or attention/PIP-informed sequence models.
   - Not recommended as the first server extraction unless storage and model size are already planned.

5. **Tier E: trans-network module set**
   - Sources: ST18, ST20, ST21, ST23, ST29.
   - Purpose: downstream multi-omics modules: ligand-receptor, PPI, inflammasome,
     reciprocal trans effects, and protein-protein colocalized axes.

6. **Tier F: eQTL-colocalized cis set**
   - Source: ST26.
   - Purpose: classify GIP signals likely mediated by tissue gene expression and
     support transcriptome-proteome interpretation.

## Suggested model families for cancer trajectory work

| Model family | Variant input | Biological interpretation |
|---|---|---|
| GIP-cis proteotype | ST10/ST16 cis variants | Direct genetically regulated protein abundance; easiest to interpret causally. |
| GIP-trans regulatory proteotype | ST10/ST20/ST21/ST23/ST29 trans variants | Network-level regulation; useful for immune, inflammation, receptor-ligand, and systemic cancer biology. |
| Functional-PAV proteotype | ST12/ST14/ST17 variants | Protein-altering or splice-impact mechanisms; prioritize for mechanistic validation. |
| Regulatory-noncoding proteotype | ST13 plus cis pQTLs | Enhancer/promoter/TF motif effects; bridge to epigenomics and chromatin annotations. |
| eQTL-colocalized GIP | ST26 variants/proteins | Expression-mediated protein regulation; bridge to GTEx/tissue-specific transcriptomics. |
| Fine-mapped probabilistic GIP | ST16 top variants and credible sets | Model uncertainty around causal variants; use PIP/credible set size as priors or weights. |
| Cancer-panel GIP | ST3 Oncology/Oncology_II proteins joined to pQTL variants | Cancer-focused protein space; useful for pre-cancer to cancer to metastasis trajectory discussion. |

## Coordinate caution for UKBB extraction

The main pQTL `Variant ID` format is `CHROM:GENPOS(hg37):A0:A1:imp:v1`.
The same tables often also include a separate `GENPOS (hg38)` column. For UKBB
imputed genotype extraction, prefer the hg37 position embedded in `Variant ID`
unless your server-side genotype resource is explicitly indexed on GRCh38.
Do not mix hg37 and hg38 positions without liftover.

## Practical first server extraction

Start with:

1. ST10 lead variants: all, plus cis/trans labels.
2. ST16 top-PIP variants: add if not already present in ST10.
3. ST12/ST13/ST17/ST14 functional subsets: force-include even if absent from
   the compact lead list.
4. ST26 eQTL-colocalized cis variants: add for multi-omics interpretation.

Keep ST16 full credible sets as a second extraction batch.
