# External labels

Labels downloaded from outside the GDC export. Raw downloads are kept unchanged next to the cleaned tables.

## PAM50 intrinsic subtypes: `pam50_tcga_pancan_atlas_2018.tsv`

- **Source:** cBioPortal study `brca_tcga_pan_can_atlas_2018` ("Breast Invasive Carcinoma (TCGA, PanCancer Atlas)", TCGA, Cell 2018), clinical attribute `SUBTYPE`. These are the TCGA PanCancer Atlas subtype calls.
- **Downloaded:** 2026-11-22 from `https://www.cbioportal.org/api/studies/brca_tcga_pan_can_atlas_2018/clinical-data?clinicalDataType=PATIENT&attributeId=SUBTYPE`. Raw response: `pam50_cbioportal_raw_2026-11-22.json`.
- **Level:** one call per **patient** (`case_id` = TCGA case barcode, same as `Case ID` in `Data/gdc_sample_sheet.2026-09-28.tsv` and `submitter_id` in the clinical JSON). Patients with two or three tumour files get the same call for each.
- **Columns:** `case_id`, `pam50` (LumA, LumB, Her2, Basal, Normal), `cbioportal_value` (original string, e.g. `BRCA_LumA`), `in_local_tumour_sheet`.
- **Coverage:** 981 patients, all in the local tumour sample sheet. 114 of the 1,095 local tumour patients have no call. 105 of the 113 normal-tissue patients have a call (it describes their tumour, not the normal sample).
- **Class counts:** LumA 499, LumB 197, Basal 171, Her2 78, Normal 36 (Normal = PAM50 "Normal-like", a tumour class, not normal tissue).
