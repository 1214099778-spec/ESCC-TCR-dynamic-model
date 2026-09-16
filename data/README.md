# Input data (not distributed)

Only this README belongs in version control. Obtain authorized inputs independently and keep them local. No real or synthetic participant examples are supplied.

| Local filename | Required fields |
| --- | --- |
| `discovery_clonotypes.tsv`, `validation_clonotypes.tsv`, `external_clonotypes.tsv` | `Patient_ID`, `Timepoint`, `CDR3b`, `TRBV`, `TRBJ`, `count` |
| `discovery_cluster_dictionary.tsv` | `Cluster_Tag`, `CDR3b`, `TRBV`, `TRBJ` |
| `discovery_clinical.csv`, `validation_clinical.csv` | `Patient_ID`, `Target`; PET columns below for multimodal analysis |
| `external_clinical.csv` | `Patient_ID`, `TRS`; optional survival columns |

`Patient_ID` is a schema field, not an actual participant identifier distributed in this repository. Use only approved pseudonymous linkage values in local inputs. Clinical tables should have one row per participant. `Target` is 1 for pCR and 0 for non-pCR. Paired `Pre` and `Post` timepoints are required. Input counts and missing values require study-specific quality control.

PET candidate fields are `SUVmax_post`, `SUVmean_post`, `SUVTBR_post`, `TLG_post`, and `MTV_post`. The final multimodal model uses `MTV_post`. Optional survival fields are `PFS_time`, `PFS_event`, `CSS_time`, `CSS_event`, and external `OS_time`, `OS_event`. Time units and event coding must be confirmed against the study definitions.

External `TRS` labels must be TRS0/TRS1/TRS2/TRS3: TRS0/1 maps to response 1 and TRS2/3 to 0. Other labels become missing.

Processed bulk TCR data: OMIX020653. Single-cell RNA/TCR data: PRJCA028740. Deposition is not a statement of public release or permission to redistribute. Raw sequencing, clinical records, spreadsheets, model objects, and participant-level supplementary files must remain outside version control.

See the workflow review notes concerning overlapping `Target` columns before running script 04.
