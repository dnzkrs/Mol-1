# Mol-1 session: activity logs

4 activities, written after the work, in the order it happened.

| # | Activity | What it covers |
| --- | --- | --- |
| 01 | [Inspect the Mol-1 repository](01-inspect-mol1-repository.md) | A read-only survey of the Mol-1 structure-prediction codebase that mapped its layout, its Protenix lineage and its install paths, and found that CI never ran because it watched a branch that does not exist |
| 02 | [Fix the CI trigger to main](02-fix-ci-trigger.md) | A two-line change to `.github/workflows/ci.yml` pointed the lint-and-test workflow at `main` instead of the non-existent `master`; it was committed and pushed to the session branch but not yet merged |
| 03 | [Predict the PhyA structure on Colab](03-predict-phya-on-colab.md) | Mol-1 predicted the 319-residue mature *Selenomonas ruminantium* phytase on a Colab T4 GPU (pTM 0.96, mean pLDDT about 95), after four failures were worked through: a hidden error, a missing `ninja`, a 502 from the MSA server and a LayerNorm build without T4 code |
| 04 | [Compare the prediction with the PhyA crystal structures](04-compare-with-crystal-structures.md) | Against 11 chains from five crystal structures, the Mol-1 model scored a median 0.98 Å Cα RMSD and 0.987 lDDT, as close as the crystals are to one another, with the catalytic P-loop within 0.54 Å all-atom; the structures were probably in its training data |

Figures live in `figures/` and are referenced relatively, so keep the
folder next to the Markdown files when moving them.
