# پروژه کارشناسی · Breast Cancer Gene-Expression Project

Last updated 2026-11-22 (full folder review: every file listed, every expression TSV header-checked, sheets, metadata and clinical JSON parsed). Project owner Mahdyar Sadeghi. The advisor referred to in project conversations is دکتر قربانعلی.

Facts marked **(Drive)** come from the Google Drive folder "G4 Uni Project" and were not re-checked locally. Facts marked **(verified locally)** were recomputed from the files in this folder.

## Start here

- This undergraduate computer science and bioinformatics project uses processed TCGA-BRCA bulk RNA-seq and existing labels to study breast-cancer classification. The intended outcome is a reproducible research analysis with a documented cohort, interpretable baseline models, and a written report and presentation.
- Start technical work in [`eda.ipynb`](eda.ipynb). The data are downloaded and initial exploratory data analysis (EDA) runs. The joined expression matrix, final target labels, evaluation splits, feature-selection runs, and trained classifiers are still pending locally.
- The written plan specifies **ductal versus lobular histology**. Recent paper work and classifier discussions focus on **five PAM50 intrinsic molecular subtypes**. **Update 2026-11-22: Mahdyar is moving to PAM50 subtypes, but tentatively ("don't make it certain for now"). Work on subtypes, keep the histology route documented and reversible, and do not rewrite the proposal or delete histology material until he confirms.** His next task is feature selection for the subtypes (see "Where to resume"). File inventory, sample mapping, and expression QC can proceed for either target.

## Research question and goals

| Direction | Question and label source | Status |
| --- | --- | --- |
| Histology | Can expression distinguish invasive ductal carcinoma from invasive lobular carcinoma? The plan uses recorded pathology diagnoses from TCGA-BRCA and restricts the cohort to unambiguously labelled primary tumours. | Original written scope. The reference document explicitly leaves PAM50 prediction outside this version of the project. [Google Doc](https://docs.google.com/document/d/1_PFOLz06WVv_oO43swtqO2S7VYTTyh255NwNlLqwylM/edit?usp=drivesdk) |
| Intrinsic molecular subtyping | Can a model recover published PAM50-derived calls for Luminal A, Luminal B, HER2-enriched, Basal-like, and Normal-like? This needs a documented subtype-label source, an identifier match to local samples, and coverage counts for each class. | **Current tentative direction (2026-11-22), not final.** The paper presentation and modeling chats explore it. Formal replacement of the written histology scope is unconfirmed. [Notion presentation notes](https://app.notion.com/p/404e0a24205d80c0b43ae98cb7f532a9?pvs=204) |

**Written plan (Drive).** The General Reference doc frames this as a reproduction study: invasive ductal (IDC) vs invasive lobular (ILC) carcinoma from TCGA-BRCA expression, comparing logistic regression, SVM and random forest against a majority-class baseline, with leakage-safe cross-validation (all preprocessing and gene selection inside training folds) and imbalance-aware metrics. METABRIC external validation is optional. PAM50 subtype prediction is explicitly out of scope in that plan, which is why the recent Cascianelli 2020 slide decks do not by themselves change the target.

Histology, receptor status such as ER or HER2, and PAM50 molecular subtype are separate targets. Normal-like is a PAM50 class, while the files under `DataNormal/` contain normal tissue samples. For either direction, document how the cohort and labels were built, compare simple models on the same held-out patients, and interpret influential genes as associations. The planned deliverables include the proposal, a data dictionary with sample eligibility and label rules, reproducible analysis code, evaluation figures and tables, and a final report and presentation. The completed paper presentation is a literature-review milestone. Project results and the final project presentation are separate deliverables.

## Local workspace

Root on Mahdyar’s Mac is `/Users/mahdyar/uniproject` (`~/uniproject`). Paths below are relative to that root. The directory was checked as a local workspace with no Git repository.

| Path | Contents and use |
| --- | --- |
| [`eda.ipynb`](eda.ipynb) | Section 1 (six code cells, from October): Loads and combines the two sample sheets, counts samples by tissue type, opens the tumour clinical JSON, counts distinct tumour expression files per case, and previews one expression file. All six cells were re-executed in order on 2026-11-22. The notebook currently has a hard-coded root path. Section 2, "Labels available in the data" (added and run 2026-11-22), explains every label the local files have and builds a per-patient `labels` table (histology, ER, PR, HER2, receptor group, stage, T/N/M, vital status) with count tables. |
| `Data/<File ID>/<File Name>` | GDC TCGA-BRCA tumour expression files named `*.rna_seq.augmented_star_gene_counts.tsv`. The primary-tumour cohort is defined by `Data/gdc_sample_sheet.2026-09-28.tsv`. The directory also contains GDC logs, extra expression files, and clinical/biospecimen supplement directories. |
| `Data/clinical.cart.2026-09-28.json` | Case-level clinical records with nested diagnoses, demographics, and follow-up. `Data/biospecimen.cart.2026-09-28.json` contains the sample hierarchy. `Data/metadata.cart.2026-09-28.json` and `Data/gdc_manifest.2026-09-28.123514.txt` describe the primary-tumour export. |
| `DataNormal/<File ID>/<File Name>` | Normal-tissue expression files with the same STAR TSV layout. The sheet, manifest, and file metadata are `gdc_sample_sheet.normal.2026-10-13.tsv`, `gdc_manifest.normal.2026-10-13.txt`, and `metadata.cart.normal.2026-10-13.json`. The accompanying case exports are `clinical.cart.normal.2026-09-28.json` and `biospecimen.cart.normal.2026-09-28.json`. |
| `Data/<File ID>/annotations.txt` | GDC case annotations present for 33 tumour files, all in the sample sheet. Mostly prior-malignancy or prior-treatment notifications, plus barcode corrections. Two flag case `TCGA-A7-A0DC` (samples `-01A` and `-01B`) as "Tumor is DCIS and not invasive", so exclude that case from an invasive IDC/ILC cohort. Review the rest when writing eligibility rules. |
| `Data/<File ID>/logs/*.parcel` | One per tumour file (1,118). Binary pickled checksum logs from the GDC download client. Not analysis data; ignore. |
| `Data/clinical.cart.2026-09-28/` | Flattened clinical TSVs: `clinical.tsv`, `follow_up.tsv`, `pathology_detail.tsv`, `exposure.tsv`, `family_history.tsv`. |
| `Data/biospecimen.cart.2026-09-28/` | Flattened biospecimen TSVs: `sample.tsv`, `portion.tsv`, `analyte.tsv` (includes RNA integrity), `aliquot.tsv`, `slide.tsv`. |
| `Labels/` | External labels, documented in `Labels/README.md`. `pam50_tcga_pancan_atlas_2018.tsv` holds TCGA PanCancer Atlas PAM50 calls from cBioPortal (downloaded 2026-11-22): one row per patient, 981 patients, all in the tumour sheet (LumA 499, LumB 197, Basal 171, Her2 78, Normal-like 36). 114 local tumour patients have no call. The raw API response is kept beside it. |
| `.DS_Store` (root, `Data/`, three UUID dirs) | macOS Finder metadata. Ignore. |
| `.venv/` | Local Python 3.14.6 environment (macOS). The notebook is configured for its IPython kernel. Installed: pandas 3.0, numpy 2.5, scipy, AnnData 0.13, h5py, zarr, marimo 0.25, ipykernel. **Not installed: scikit-learn and matplotlib/seaborn**, which the modelling and plotting steps need. No project dependency lockfile was found. |

### Verified data inventory

| Selected by | Expression files in sheet | Distinct sample IDs | Distinct case IDs |
| --- | ---: | ---: | ---: |
| Primary-tumour sample sheet | 1,111 | 1,106 | 1,095 |
| Normal sample sheet | 113 | 113 | 113 |
| Both sheets combined | 1,224 | 1,219 | 1,095 |

Every normal case in the local sheet also appears in the primary-tumour sheet. The notebook finds 1,084 primary-tumour cases with one file, six with two files, and five with three files. All 1,224 sheet-listed expression files exist. There are 1,118 count TSVs under `Data/`, including seven outside the primary-tumour sheet that the Notion inventory describes as metastatic. Use the sample sheets to select the cohort. [Notion data notes](https://app.notion.com/p/3f3e0a24205d80b6b432f8fd2640d662?pvs=204)

More checks (verified locally, 2026-11-22):

- The seven extra tumour TSVs under `Data/` are not in the sample sheet, the manifest, or `metadata.cart.2026-09-28.json` (all three list exactly the 1,111 sheet files), so their sample type cannot be confirmed from local files. Their File IDs start `216378aa`, `e0ec172b`, `a10463d8`, `64f43094`, `d2659fd9`, `11396f70`, `ecbec59f`. Exclude them.
- Tumour sheet: every row is `Tissue Type = Tumor`, `Tumor Descriptor = Primary`. Sample barcode suffixes are `01A` (1,086), `01B` (24), `01C` (1). Preservation is OCT 679, FFPE 13, Unknown 419; consider flagging FFPE samples in QC. Normal sheet suffixes are `11A` (99) and `11B` (14).
- All 1,231 count TSVs (1,118 tumour, 113 normal) have the same layout: the GENCODE v36 comment line, a 9-column header, 4 `N_` summary rows and **60,660 genes**. None failed to open.
- The clinical JSON has exactly 1,095 cases, one per tumour-sheet case, matched on `submitter_id`. 913 cases have one diagnosis record, 181 have several (metastasis, recurrence, prior or synchronous primaries), and one has none, so pick the record with `classification_of_tumor = primary`.
- **Preliminary histology counts** (one per case, from the primary diagnosis record): `Infiltrating duct carcinoma, NOS` 777, `Lobular carcinoma, NOS` 201, `Infiltrating duct and lobular carcinoma` 28, ductal or lobular "mixed with other types" 25, other special types (mucinous, metaplastic, medullary and so on) 63, no diagnosis 1. A pure IDC vs ILC cohort is therefore about 777 vs 201 cases (roughly 79% / 21%) before the DCIS and annotation exclusions. This is a working count, not the final label table.
- **PAM50 labels are not in the GDC files.** They were downloaded separately into `Labels/pam50_tcga_pancan_atlas_2018.tsv` (see the Local workspace table). Join on TCGA case barcode.
- **Receptor labels** live in the clinical JSON under `follow_ups[].molecular_tests[]` (`gene_symbol` ESR1 / PGR / ERBB2, `molecular_analysis_method` IHC or FISH, `test_result`). Per patient (any Positive wins; HER2 uses FISH when present): ER 808+ / 236− / 51 missing or equivocal; PR 699+ / 342− / 54; HER2 182+ / 768− / 145. Receptor groups: HR+/HER2− 605, HER2+ 181, triple-negative 159, unknown 150. These only roughly track PAM50 and are not a substitute.
- `DataNormal/clinical.cart.normal.2026-09-28.json` holds 113 cases. The `2026-09-28` in the normal clinical/biospecimen file names does not match their 2026-10-13 download date; treat the normal sample sheet as authoritative.

- Use `File ID` and `File Name` to resolve each expression path. Keep `Sample ID` and `Case ID` separately. The sheet’s TCGA case barcode corresponds to the clinical JSON’s `submitter_id`. The JSON’s `case_id` is a UUID. Validate joins and report duplicates or ambiguous diagnoses before producing class counts.
- The GENCODE v36 STAR TSV starts with a gene-model comment. Skip that line when reading the headers and drop the four STAR summary rows whose `gene_id` starts with `N_`. The last cell in the current notebook previews the file without skipping the comment, so its displayed header is misparsed. The verified example parsed correctly with `pd.read_csv(path, sep="\t", skiprows=1)`.
- Expression fields include raw `unstranded`, `stranded_first`, and `stranded_second` counts, plus `tpm_unstranded`, `fpkm_unstranded`, and `fpkm_uq_unstranded`. Preserve the original downloads, match genes by Ensembl `gene_id`, and record the expression measure and gene filters used in each analysis.

## Documents and references

| Where | What is stored |
| --- | --- |
| [G4 Uni Project · Google Drive](https://drive.google.com/drive/folders/1c3ttUwuIJz1uyHD46zVg56hofgMhKNC_) | Main cloud folder for the reference plan, comparison spreadsheet, papers, and slide files. |
| [Breast Cancer Gene-Expression Project — General Reference · Google Doc](https://docs.google.com/document/d/1_PFOLz06WVv_oO43swtqO2S7VYTTyh255NwNlLqwylM/edit?usp=drivesdk) | Original ductal-versus-lobular plan, cohort rules, workflow, deliverables, and evaluation criteria. |
| [Gene Expression — ML Paper Comparison · Google Sheet](https://docs.google.com/spreadsheets/d/1VJ1VncEPPgDcQn8frjgRUvAD7eU_JQhXysVFg0kfiWA/edit?usp=drivesdk) | Literature comparison spreadsheet covering papers 1 to 12 **(Drive)**. Paper 11, on molecular classification of lobular carcinoma, is the closest benchmark for the histology task and was still unread. Papers 08 (Cascianelli 2020) and 13 (PAM50) are in the Papers folder but not in the sheet. |
| [Breast Project · Notion](https://app.notion.com/p/3f3e0a24205d80b6b432f8fd2640d662?pvs=204) | Literature and data notes, including the STAR columns and cohort inventory. |
| [Paper Features · Notion database](https://app.notion.com/p/3f3e0a24205d80f48bf3d6f9b87e65b1?pvs=204) | Paper-by-paper notes on prediction target, feature selection, models, and reading status. Includes [Cascianelli et al.](https://app.notion.com/p/3f3e0a24205d80da962def2e15ba09de?pvs=204), [FSFSP](https://app.notion.com/p/3f3e0a24205d80b2b46ed8e34f801e87?pvs=204), [TNBC and CNNs](https://app.notion.com/p/3f3e0a24205d80bd8f8ae2e32f80c47f?pvs=204), and [ER-status classification](https://app.notion.com/p/403e0a24205d802d99d3c274b935f5cd?pvs=204). |
| [Subtype paper presentation · Notion](https://app.notion.com/p/404e0a24205d80c0b43ae98cb7f532a9?pvs=204) | Detailed notes and slide outline on PAM50 reference centering, AWCA, single-sample prediction, mLR, gene selection, results, and limitations. |
| [Paper presentation PDF](https://drive.google.com/file/d/1lJwH4L87Ciwbif_fCxk5iJb9Zw8o7PwU/view?usp=drivesdk), [English slides](https://docs.google.com/presentation/d/1OIW49leL8ZefzAlDTuFMJzsySvGKHgGf/edit?usp=drivesdk), [Persian slides](https://docs.google.com/presentation/d/1bjbOlkz00POaAZX9ZuyrC7Kr9lSWrqIF/edit?usp=drivesdk) | Existing slide deliverables for the literature presentation. |
| [Papers · Google Drive](https://drive.google.com/drive/folders/1-NIBGyRObZ9fr59slbI-rlgzobgxsJNS) | Source PDFs, including [Cascianelli et al., RNA-seq intrinsic subtyping](https://drive.google.com/file/d/1GJdAaKc8NMZoJ7v4-G_crIJSHxX8Lqc7/view?usp=drivesdk), the [original PAM50 paper](https://drive.google.com/file/d/1569yTe_1Z0-W4yjDCVm0Elp4Ph-hwiCL/view?usp=drivesdk), and [molecular classification of lobular carcinoma](https://drive.google.com/file/d/1TFuyr3iqEA0pgJeUxO5PFFpZxRjxSeIM/view?usp=drivesdk). |
| [🎓Uni · TickTick](https://ticktick.com/webapp/#p/65d3b9710232d10263d0b3d1) | Task tracker for the project and other coursework. Use the breast-project entries below. |

## Task history

Dates for completed tasks reflect TickTick. Completion refers to the task itself.

| Date | Task | Outcome |
| --- | --- | --- |
| 2026-08-09 | مشخص کردن پروژه | Topic definition was marked complete after an in-person discussion was planned. |
| 2026-08-09 | هماهنگی با دکتر زارع برای صحبت در رابطه با پروژه کارشناسی | Early coordination about the undergraduate project was marked complete. |
| 2026-08-15 | آماده کردن پروپوزال برای پروژه | Proposal preparation was marked complete. The Google Doc records the histology plan. Proposal submission and approval are unverified. |
| 2026-10-01 | مطالعه مقاله های مربوط به برست | Initial reading on breast cancer and gene expression was marked complete. |
| 2026-10-15 | Read another paper | Reading milestone marked complete. The task title gives no specific citation. |
| 2026-10-19 | Read the 3rd paper | Reading milestone marked complete. The task title gives no specific citation. |
| 2026-10-19 | Create a table comparing methods of papers | Compared algorithms and data-reduction methods. The Google Sheet and Notion database hold the literature comparisons. |
| 2026-10-19 | شناسایی داده | Identified the GDC file layout and expression columns and performed initial EDA in the notebook. |
| 2026-10-19 | روش های مختلف data reduction و کاهش ابعاد رو ببینم. | Reviewed ways to reduce genes, including the FSFSP reading notes on Boruta, RNA-type grouping, and tree-based ranking. This was literature work. |
| 2026-11-21 | تکمیل ارائه مقاله بر اساس outline | Marked complete. The Notion outline and Drive copies of the PDF and slide decks are available. [TickTick task](https://ticktick.com/webapp/#p/65d3b9710232d10263d0b3d1/tasks/6ae9ea9b3fad110ec742ce6b) |
| Open | Read Ciriello et al. 2015 (Cell), the TCGA lobular paper | Requested 2026-11-22. It has PAM50 calls and expert-reviewed IDC/ILC histology for about 800 TCGA tumours, useful as a cross-check for either target. Not yet added to TickTick (no TickTick access from Claude). |
| Open | Feature engineering for gene expressions | Define the expression representation and build the sample-by-gene matrix and metadata needed for experiments. |
| Open | Apply a DEG based feature selection on TCGA data | Try a differential-expression-based feature-selection workflow once the cohort and labels are defined. The task notes mention an R package and a volcano plot. |

## Project chats and discussion notes

This log covers the project discussions available in the conversation excerpts and project notes. Original chat permalinks were unavailable. The matching project artifacts are linked where available.

| Date | Discussion | Takeaway |
| --- | --- | --- |
| Sep–Oct 2026 | Daily Schedule Forecast and advisor preparation | Project sessions were used for reading papers, comparing methods, understanding gene expression files, and EDA. The user also requested dedicated time for article and proposal work. The meeting with دکتر قربانعلی is relevant to the thesis. |
| 2026-10-17 | Continuing the university breast-cancer project | Mahdyar wanted to continue the existing advisor-supported project. The discussion recommended agreeing on the next concrete milestone with the advisor. |
| Oct–Nov 2026 | Subtype paper presentation | Prepared an outline explaining why PAM50 calls depend on the reference cohort, how AWCA averages subtype means to construct a more stable reference, and how regularized multiclass logistic regression compares with PAM50-based approaches. The paper’s reported metrics measure agreement with its published labels. [Presentation notes](https://app.notion.com/p/404e0a24205d80c0b43ae98cb7f532a9?pvs=204) |
| 2026-11-22 | Data Warehouse Comparison, university project follow-up | Mahdyar wanted spreadsheet-like exploration plus Python. The suggested workflow was pandas for validated ingestion, AnnData `.h5ad` for the matrix and metadata, marimo for interactive exploration, and optionally JupyterLab with Mito for spreadsheet-style interaction. An analysis-ready dataset was proposed. |
| 2026-11-22 | Train Cancer Subtype Model | Mahdyar meant Jev-style decision models. The discussion favored a specialized classifier for the numeric RNA-seq matrix. Elastic-net logistic regression was recommended first, then a linear SVM as a comparison. This is a modeling recommendation. |
| 2026-11-22 | Differential Expression Explained | Worked through why five subtypes produce ten pairwise contrasts. Take up to the top `N` genes from each contrast and deduplicate their union to form `limmaN`. In the published paper, `limma50` yielded 277 genes. The authors tested 11 values of `N` from 10 to 1,000. These are literature results. [Paper PDF](https://drive.google.com/file/d/1GJdAaKc8NMZoJ7v4-G_crIJSHxX8Lqc7/view?usp=drivesdk) |

## Where to resume

0. Housekeeping: install scikit-learn and matplotlib (or seaborn) in `.venv`, replace the notebook's hard-coded `ROOT` with a relative path, and fix the last cell to use `skiprows=1`. For the histology task, read paper 11 in the Papers folder, the closest published benchmark.
**Next task (Mahdyar, 2026-11-22): feature selection on the data for PAM50 subtypes.** It needs, in order: a subtype label source matched to the local samples (step 1), the expression matrix (step 2), and patient-disjoint splits (step 3), because supervised selection such as limma/DEG must run on training data only. Candidate feature sets to compare are the PAM50 50-gene panel, the `limmaN` union of top genes from the ten pairwise contrasts, and unsupervised filters such as variance or expression thresholds. Not started.

1. Confirm the primary research question and label source (tentatively PAM50 subtypes). For histology, document how eligible diagnosis records map to ductal or lobular cases and resolve mixed or ambiguous records explicitly. For intrinsic subtype, obtain published PAM50-derived labels and verify the match at the correct TCGA identifier level. Record unmatched samples and class counts before defining the cohort.
2. Create a reproducible loader from the sample sheets and expression TSVs. Resolve case, sample, and file relationships, handle repeated tumours with a documented rule, and write an analysis-ready expression matrix and metadata. The current preference is a notebook-friendly workflow. AnnData and marimo are available locally.
3. Run cohort and expression QC, including missingness, class balance, and PCA where useful. Create patient-disjoint train and test sets and use patient-aware cross-validation when repeated samples are retained. Fit preprocessing, reference construction where applicable, supervised gene selection, and tuning on training data only, inside each training fold when cross-validating.
4. Establish a majority-class baseline, elastic-net logistic regression, and a linear SVM. Evaluate with balanced accuracy, macro-F1, per-class precision and recall, and a confusion matrix. For the molecular task, a PAM50 gene-panel comparison can follow label verification. Keep external validation and more complex models optional.
5. Compare feature sets, including the open DEG/limma experiment if it fits the chosen target. Treat the paper’s ten pairwise contrasts as specific to a five-class setup. Report label coverage, split rules, feature-selection results, model errors, and limitations alongside the final figures and write-up.

## Maintaining this file

Keep this handoff current when the target is confirmed, a new dataset or output is created, or an experiment is run. Record its path, inputs, and status so another agent can reproduce it. When you change it, update the date at the top, mark where a fact came from (verified locally, Drive, Notion, or a chat), and add a line to the log below.

| Date | Change |
| --- | --- |
| 2026-11-22 | Full folder review. Added annotation, log and supplement-file descriptions, the extra-file check, the STAR layout check across all 1,231 files, sample-type and preservation counts, preliminary IDC/ILC counts, the DCIS exclusion, missing scikit-learn, and the Drive plan details (RF, METABRIC, paper sheet, paper 11). |
| 2026-11-22 | Recorded the tentative switch to PAM50 subtypes and subtype feature selection as the next task. |
| 2026-11-22 | Added section 2 to `eda.ipynb` (labels available, with explanations and counts); recorded that no PAM50 labels exist locally and the receptor-status counts. |
| 2026-11-22 | Added `Labels/` with the TCGA PanCancer Atlas PAM50 calls; added reading Ciriello 2015 to open tasks. |
