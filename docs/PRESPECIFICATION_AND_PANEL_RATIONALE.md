# Pre-Specification and Literature-Guided Panel Documentation

## Practical answer

If the panel was not formally pre-registered before all modeling, do not claim formal pre-registration. Instead, document that the panel was literature-guided and fixed before the final 100-repeat model-comparison run.

The strongest defensible wording is:

> The literature-guided panels were defined from prior ALS CSF/progression biology and assay availability before the final repeated cross-validation analysis. Because the study is small and internally validated only, all findings are interpreted as exploratory and require external validation.

## How to support this in the manuscript/repository

1. Include a dated literature-review document or PDF used to define the candidate proteins.
2. Include a panel-rationale table with one row per protein: gene/protein label, mechanism axis, evidence tier, assay-matched feature, and reason for inclusion.
3. Include the exact final code commit or timestamp before the final 100-repeat run.
4. Include `run_metadata.json`, which records the fixed model set, panel terms, and panel sizes used for the run.
5. Keep all old/stale files with outdated ad-hoc panel names out of the repository and supplement.
6. Avoid wording that implies post-result optimization, such as "best panel," "optimized panel," or "updated after performance review."

## Recommended manuscript wording

> The primary protein model used a predefined seven-protein mechanism panel comprising CHIT1, RBP4, AMBP, SERPINA4, CST3, TNR, and NCAN. This panel was defined before final model evaluation and was designed to represent ALS-relevant biological mechanisms rather than to maximize model performance. A broader 14-protein literature-supported comparator panel and one exploratory training-fold data-driven top-10 protein panel were also evaluated.

## Wording to avoid

- ad hoc panel chosen after inspecting results
- manually revised after performance review
- optimized final panel
- validated classifier
- clinical-grade prediction model
- definitive biomarker discovery

## NCAN caveat

For NCAN, state explicitly that the assay likely measures total/intact NCAN abundance and does not specifically isolate the disease-linked 150 kDa NCAN cleavage fragment reported in ALS CSF literature. This makes NCAN biologically motivated but not a direct replication of that fragment-level finding.
