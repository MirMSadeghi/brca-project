# /// script
# requires-python = ">=3.14"
# dependencies = [
#     "anndata>=0.13.4",
#     "marimo>=0.23.3",
#     "matplotlib>=3.11.2",
#     "numpy>=2.5.5",
#     "pandas>=3.1.1",
#     "plotly>=7.2.0",
#     "pylimma>=0.3.0",
#     "scikit-learn>=1.9.1",
#     "seaborn>=0.13.2",
# ]
# ///

import marimo

__generated_with = "0.25.2"
app = marimo.App()


@app.cell
def _():
    import marimo as mo
    import anndata as ad
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    import sklearn as sk
    import pylimma as pl    # python port of limma
    from itertools import combinations
    import plotly.express as px

    PROJECT_ROOT = "/Users/mahdyar/uniproject"
    DATA_FILE = PROJECT_ROOT + "/Processed/tcga_brca.h5ad"
    return DATA_FILE, ad, combinations, mo, np, pd, pl, plt, px, sk


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    We create our train/test split. Each patient has one file, so we can split by patient. We only include patients that have both pam50 and tissue_type information in our initial cohort.
    """)
    return


@app.cell
def _(DATA_FILE, ad):

    adata = ad.read_h5ad(DATA_FILE)

    # get cohort (patients with pam50 and tumor)
    cohort_adata = adata[(adata.obs["pam50"].notnull()) & (adata.obs["tissue_type"] == "Tumor")]

    # table showing files per patient
    files_per_patient = cohort_adata.obs.groupby("case_id").size()
    files_per_patient.value_counts()
    return (cohort_adata,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Train/Test Split
    """)
    return


@app.cell
def _(cohort_adata, sk):

    # split
    train_ids, test_ids = sk.model_selection.train_test_split(
        cohort_adata.obs_names,
        test_size=0.2,
        stratify=cohort_adata.obs["pam50"],
        random_state=42,
    )

    train_adata = cohort_adata[train_ids].copy()
    test_adata = cohort_adata[test_ids].copy()
    return (train_adata,)


@app.cell
def _(train_adata):
    # use only protein coding genes for feature selection
    protein_coding_genes = train_adata.var.index[train_adata.var["gene_type"] == "protein_coding"]
    len(protein_coding_genes)
    return (protein_coding_genes,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Get log2 mean of all paitents for each gene
    Also filter out genes that have a low expression. mean < 1.
    """)
    return


@app.cell
def _(np, pd, protein_coding_genes, train_adata):
    # get log2(tpm + 1)
    train_log = pd.DataFrame(
        np.log2(train_adata[:, protein_coding_genes].layers["tpm"] + 1),
        index=train_adata.obs.index,
        columns=protein_coding_genes,
    )
    train_log.head()
    return (train_log,)


@app.cell
def _(pd, train_adata, train_log):
    # use genes with log2(tpm + 1) mean of at least1
    train_expression_mean = train_log.mean(axis=0)
    meaningful_genes_index = train_expression_mean[train_expression_mean >= 1].index

    # show mean expression of meaningful genes
    pd.DataFrame({
        "gene_name": train_adata.var.loc[meaningful_genes_index, "gene_name"],
        "mean_log2_tpm_expression": train_expression_mean[meaningful_genes_index],
    }).sort_values("mean_log2_tpm_expression", ascending=False)
    return (meaningful_genes_index,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    **`train_expression`**: The training set expression data, with samples as rows and genes as columns. `sample * gene`
    """)
    return


@app.cell
def _(meaningful_genes_index, train_log):
    train_expressions = train_log[meaningful_genes_index]
    train_expressions.head()
    return (train_expressions,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Define subtype contrasts
    **`limma_design`**: The design matrix for the limma model, which includes the subtype information for each sample. `sample * subtype`
    """)
    return


@app.cell
def _(pd, train_adata):
    SHORT_SUBTYPE_NAMES = {
        "Luminal A": "LumA",
        "Luminal B": "LumB",
        "HER2-enriched": "Her2",
        "Basal-like": "Basal",
        "Normal-like": "Normal",
    }

    train_groups = train_adata.obs["pam50"].map(SHORT_SUBTYPE_NAMES)

    # one hot encoding
    limma_design = pd.get_dummies(train_groups).astype(float)
    limma_design.head()
    return (limma_design,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Here we create the combinations limma should evaluate. Each column is a combination and that column shows the formula for that combination. For example, the first column is the LumA-LumB comparison, and the formula is "LumA - LumB".

    Limma accepts contrasts as "TypeA-TypeB".
    """)
    return


@app.cell
def _(combinations, limma_design, pl):
    # create acceptable contrast format
    contrast_names = [f"{a}-{b}" for a, b in combinations(list(limma_design.columns), 2)]

    # feed contrasts to limma
    limma_contrasts = pl.make_contrasts(*contrast_names, levels=list(limma_design.columns))
    limma_contrasts
    return contrast_names, limma_contrasts


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Fit limma on the `limma_design` matrix.
    The genes are in `train_expressions`, which is a matrix of `samples * genes`. The design matrix is `limma_design`, which is a matrix of `samples * subtypes`.
    `lm_fit` wants genes as rows and samples as columns, hence the .T
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Each contrast will give a list of all genes and related statistics. Here we show the LumA-LumB contrast as an example. The columns are:
    - **gene_id**: the gene identifier
    - **log_fc**: log2 fold change between the two subtypes in the contrast. log2 of LumA mean minus the LumB mean. +1 means a twice as high expression in LumA and −1 means twice as high in LumB.
    - **ave_expr**: the gene's average log2(TPM + 1) over all training patients. regardless of subtype.
    - **t**: the moderated t-statistic for the contrast. log_fc divided by its uncertainty after limma's variance sharing. The further from 0, the stronger the evidence. Its sign matches log_fc.
    - **p_value**: the p-value for the contrast.
    - **adj_p_value**: the p-value adjusted for multiple testing using the Benjamini-Hochberg method. This is the value we use to determine significance.
    - **b**: the B-statistic, the log-odds that the gene really differs. Above 0 means more likely than not. It ranks genes the same way as p_value.
    """)
    return


@app.cell
def _(limma_contrasts, limma_design, pl, train_adata, train_expressions):
    _fit = pl.lm_fit(train_expressions.T, limma_design)    # five subtype means per gene
    _fit = pl.contrasts_fit(_fit, limma_contrasts)  # the 10 pairwise differences
    limma_fit = pl.e_bayes(_fit, trend=True)        # moderated t-statistics and p-values
    pl.top_table(limma_fit, coef="LumA-LumB", number=10).join(train_adata.var["gene_name"])
    return (limma_fit,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Get one gene list of **significant genes** for each contrast.

    Significant genes are defined as:
    - p-value below 0.05
    - At least a two-fold difference, which is an absolute log2 fold change above 1 (log_fc=1).

    With about 776 patients nearly every real difference is significant, so the two-fold rule is what keeps the lists to a sensible length
    """)
    return


@app.cell
def _(contrast_names, limma_fit, np, pd, pl):
    # a dictionary of dataframes, one for each contrast
    de_tables = {
        name: pl.top_table(
            limma_fit, coef=name, number=np.inf, sort_by="p", p_value=0.05, lfc=1
        )
        for name in contrast_names
    }

    # get number of DE genes for each contrast from the dictionary of dataframes
    pd.Series({name: len(table) for name, table in de_tables.items()}, name="No. of DE genes")
    return (de_tables,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Build limmaN: Top N genes of each list.
    """)
    return


@app.cell
def _(de_tables, pd):
    # return gene_id of top N genes of each contrast
    def limma_n(n):
        genes = set()
        for table in de_tables.values():
            genes.update(table.index[:n])  # tables are already sorted, most significant first
        return sorted(genes)

    N_values = [10, 20, 50, 100, 200]

    pd.Series({n: len(limma_n(n)) for n in N_values}, name="genes in limmaN")
    return (limma_n,)


@app.cell
def _(limma_n, train_adata):
    limma50_table = train_adata.var.loc[limma_n(50), ["gene_name"]]
    limma50_table
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Display Volcano plots
    The x-axis shows the log2 fold change (log_fc), and the y-axis shows the -log10 p-value. Genes that are significantly differentially expressed will appear in the upper left and upper right corners of the plot.
    Genes on the right are more highly expressed in the first subtype of the contrast, and vice versa.
    The further the gene is from the x=0 line, the larger the difference in expression between the two subtypes. The further up the gene is, the more statistically significant the difference is.
    """)
    return


@app.cell
def _(limma_fit, np, pl, train_adata):
    # builds the volcano table for one contrast; used by the interactive plot and the grid below
    def volcano_table(contrast):
        volcano = pl.top_table(limma_fit, coef=contrast, number=limma_fit.nrow)
        volcano = volcano.join(train_adata.var["gene_name"])   # add gene symbols
        volcano["neg_log10_p"] = -np.log10(volcano["p_value"].clip(lower=1e-300))  # very small p-values round to 0
        volcano["DE"] = (volcano["adj_p_value"] < 0.05) & (volcano["log_fc"].abs() > 1)
        return volcano

    return (volcano_table,)


@app.cell
def _(px, volcano_table):
    # interactive volcano plot for one contrast (hover shows the gene name)
    def volcano_figure(contrast):
        figure = px.scatter(
            volcano_table(contrast).reset_index(),   # Ensembl ID becomes a normal column
            x="log_fc",
            y="neg_log10_p",
            color="DE",
            color_discrete_map={True: "crimson", False: "lightgrey"},
            hover_name="gene_name",             # bold title of the hover box
            hover_data=["gene", "adj_p_value"], # extra lines in the hover box
            render_mode="webgl",                # keeps ~13k points smooth
            title=contrast,
        )
        return figure

    return (volcano_figure,)


@app.cell
def _(contrast_names, mo):
    contrast_picker = mo.ui.dropdown(contrast_names, value=contrast_names[0], label="Contrast")
    contrast_picker
    return (contrast_picker,)


@app.cell
def _(contrast_picker, mo, volcano_figure):
    volcano_plot = mo.ui.plotly(volcano_figure(contrast_picker.value))
    volcano_plot
    return


@app.cell
def _(contrast_names, np, plt, volcano_table):
    # all 10 contrasts at once as one static image: 2 rows of 5 panels.
    # matplotlib sends a picture, not every point, so the output stays small.
    volcano_grid, _axes = plt.subplots(2, 5, figsize=(20, 8), sharex=True, sharey=True)

    # walk through the panels and the contrasts together
    for _ax, _contrast in zip(_axes.flat, contrast_names):
        _volcano = volcano_table(_contrast)
        _ax.scatter(
            _volcano["log_fc"],
            _volcano["neg_log10_p"],
            s=3,
            c=np.where(_volcano["DE"], "crimson", "lightgrey"),
        )
        _ax.axvline(-1, linestyle="--", color="grey")
        _ax.axvline(1, linestyle="--", color="grey")
        _ax.set_title(f"{_contrast}  ({_volcano['DE'].sum()} DE genes)")

    # one shared label per axis instead of ten
    volcano_grid.supxlabel("log2 fold change (first subtype minus second)")
    volcano_grid.supylabel("-log10 p-value")
    volcano_grid.tight_layout()
    volcano_grid
    return (volcano_grid,)


if __name__ == "__main__":
    app.run()
