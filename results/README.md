# Local outputs (not distributed)

Only this README is tracked. Scripts create `results/private/`, `derived_data/`, and `figures/` locally. These directories and common data formats are ignored by Git.

Outputs include locked TCR and multimodal model objects, fitted coefficients, feature rankings, PET candidate AIC comparisons, participant-level dynamic features and prediction tables. Keep all such outputs local. Script 05 calculates ROC and survival objects; it does not implement a complete file export of all plots or statistics.

Never force-add generated outputs. An ignore rule is a safeguard, not authorization to share study data.
