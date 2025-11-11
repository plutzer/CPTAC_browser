"""
CPTAC Browser - Python Shiny GUI for CPTAC Data Analysis
"""

from shiny import App, ui, render, reactive
import cptac
import pandas as pd
import numpy as np
from scipy import stats
import plotly.graph_objects as go
import plotly.express as px
from typing import Tuple, List, Optional

# Available cancer types
CANCER_TYPES = {
    "brca": "Breast Cancer",
    "coad": "Colon Adenocarcinoma",
    "hnscc": "Head and Neck Squamous Cell Carcinoma",
    "luad": "Lung Adenocarcinoma",
    "ovarian": "Ovarian Cancer",
    "ccrcc": "Clear Cell Renal Cell Carcinoma",
    "gbm": "Glioblastoma",
    "lscc": "Lung Squamous Cell Carcinoma",
    "pdac": "Pancreatic Ductal Adenocarcinoma"
}

# UI Definition
app_ui = ui.page_navbar(
    ui.nav_panel(
        "Protein Analysis",
        ui.layout_sidebar(
            ui.sidebar(
                ui.input_select(
                    "protein_cancer",
                    "Cancer Type:",
                    choices=CANCER_TYPES
                ),
                ui.input_text_area(
                    "protein_genes",
                    "Gene Names (comma-separated):",
                    placeholder="e.g., TP53, EGFR, KRAS",
                    rows=3
                ),
                ui.input_action_button("protein_run", "Run Analysis", class_="btn-primary"),
                width=300
            ),
            ui.output_ui("protein_results_ui"),
            ui.output_plot("protein_plot")
        )
    ),
    ui.nav_panel(
        "Phospho Analysis",
        ui.layout_sidebar(
            ui.sidebar(
                ui.input_select(
                    "phospho_cancer",
                    "Cancer Type:",
                    choices=CANCER_TYPES
                ),
                ui.input_text_area(
                    "phospho_query",
                    "Phosphosites/Genes (comma-separated):",
                    placeholder="e.g., EGFR_Y1068, TP53, AKT1_S473",
                    rows=3
                ),
                ui.input_switch(
                    "phospho_normalized",
                    "Normalized (protein-adjusted)",
                    value=True
                ),
                ui.input_action_button("phospho_run", "Run Analysis", class_="btn-primary"),
                width=300
            ),
            ui.output_ui("phospho_results_ui"),
            ui.output_plot("phospho_plot")
        )
    ),
    ui.nav_panel(
        "Correlation Analysis",
        ui.layout_sidebar(
            ui.sidebar(
                ui.input_select(
                    "corr_cancer",
                    "Cancer Type:",
                    choices=CANCER_TYPES
                ),
                ui.input_text_area(
                    "corr_query",
                    "Genes/Phosphosites (comma-separated):",
                    placeholder="e.g., EGFR_Y1068, TP53, AKT1_S473",
                    rows=3
                ),
                ui.input_select(
                    "corr_data_type",
                    "Data Type:",
                    choices={"phospho": "Phosphoproteomics", "proteomics": "Proteomics", "both": "Both"}
                ),
                ui.input_switch(
                    "corr_normalized",
                    "Normalized Phospho (protein-adjusted)",
                    value=True
                ),
                ui.input_action_button("corr_run", "Run Analysis", class_="btn-primary"),
                width=300
            ),
            ui.output_ui("corr_results_ui"),
            ui.output_plot("corr_plot")
        )
    ),
    title="CPTAC Browser",
    id="navbar"
)


def server(input, output, session):

    # Helper function to identify tumor vs normal samples
    def is_normal_sample(sample_id: str) -> bool:
        """Check if sample ID indicates normal tissue (ends with .N)"""
        return str(sample_id).endswith('.N')

    def get_tumor_normal_pairs(df: pd.DataFrame) -> Tuple[List[str], List[str]]:
        """Extract paired tumor and normal sample IDs"""
        normal_samples = [col for col in df.columns if is_normal_sample(col)]
        tumor_samples = []

        for normal in normal_samples:
            tumor = normal[:-2]  # Remove '.N' suffix
            if tumor in df.columns:
                tumor_samples.append(tumor)

        return tumor_samples, normal_samples

    def calculate_fold_change_stats(tumor_data: pd.Series, normal_data: pd.Series) -> dict:
        """Calculate log2 fold change and statistics"""
        tumor_clean = tumor_data.replace([np.inf, -np.inf], np.nan).dropna()
        normal_clean = normal_data.replace([np.inf, -np.inf], np.nan).dropna()

        if len(tumor_clean) < 2 or len(normal_clean) < 2:
            return None

        mean_tumor = tumor_clean.mean()
        mean_normal = normal_clean.mean()
        log2_fc = mean_tumor - mean_normal

        # Paired t-test
        _, p_value = stats.ttest_rel(tumor_clean, normal_clean)

        return {
            'log2_fold_change': log2_fc,
            'p_value': p_value,
            'mean_tumor': mean_tumor,
            'mean_normal': mean_normal,
            'n_pairs': len(tumor_clean)
        }

    def apply_fdr_correction(results_df: pd.DataFrame) -> pd.DataFrame:
        """Apply Benjamini-Hochberg FDR correction"""
        from statsmodels.stats.multitest import multipletests

        if 'p_value' in results_df.columns and len(results_df) > 0:
            valid_pvals = results_df['p_value'].dropna()
            if len(valid_pvals) > 0:
                _, pvals_corrected, _, _ = multipletests(valid_pvals, method='fdr_bh')
                results_df.loc[valid_pvals.index, 'p_value_adjusted'] = pvals_corrected

        return results_df

    # Protein Analysis
    protein_results = reactive.Value(None)

    @reactive.Effect
    @reactive.event(input.protein_run)
    def run_protein_analysis():
        cancer = input.protein_cancer()
        genes = [g.strip() for g in input.protein_genes().split(',') if g.strip()]

        if not genes:
            protein_results.set(pd.DataFrame({'error': ['Please enter gene names']}))
            return

        try:
            # Load cancer data
            if cancer == "brca":
                data = cptac.Brca()
            elif cancer == "coad":
                data = cptac.Coad()
            elif cancer == "hnscc":
                data = cptac.Hnscc()
            elif cancer == "luad":
                data = cptac.Luad()
            elif cancer == "ovarian":
                data = cptac.Ovarian()
            elif cancer == "ccrcc":
                data = cptac.Ccrcc()
            elif cancer == "gbm":
                data = cptac.Gbm()
            elif cancer == "lscc":
                data = cptac.Lscc()
            elif cancer == "pdac":
                data = cptac.Pdac()

            # Get proteomics data
            proteomics = data.get_proteomics()

            # Handle multi-index
            if isinstance(proteomics.columns, pd.MultiIndex):
                proteomics = proteomics.groupby(level=list(range(proteomics.columns.nlevels)), axis=1).mean()

            # Transpose so genes are rows and samples are columns
            proteomics = proteomics.T

            results = []
            for gene in genes:
                if gene in proteomics.index:
                    gene_data = proteomics.loc[gene]
                    tumor_samples, normal_samples = get_tumor_normal_pairs(pd.DataFrame(gene_data).T)

                    if len(tumor_samples) > 0:
                        tumor_data = gene_data[tumor_samples]
                        normal_data = gene_data[[s for s in normal_samples if s[:-2] in tumor_samples]]

                        stats_result = calculate_fold_change_stats(tumor_data, normal_data)
                        if stats_result:
                            stats_result['gene'] = gene
                            results.append(stats_result)

            if results:
                results_df = pd.DataFrame(results)
                results_df = apply_fdr_correction(results_df)
                protein_results.set(results_df)
            else:
                protein_results.set(pd.DataFrame({'error': ['No results found for the specified genes']}))

        except Exception as e:
            protein_results.set(pd.DataFrame({'error': [f'Error: {str(e)}']}))

    @output
    @render.ui
    def protein_results_ui():
        results = protein_results.get()
        if results is None:
            return ui.p("Run analysis to see results")

        if 'error' in results.columns:
            return ui.p(results['error'].iloc[0], style="color: red;")

        return ui.HTML(results.to_html(index=False, float_format=lambda x: f'{x:.4f}'))

    @output
    @render.plot
    def protein_plot():
        results = protein_results.get()
        if results is None or 'error' in results.columns:
            return None

        # Create volcano plot
        results['neg_log10_pval'] = -np.log10(results['p_value'])

        fig = go.Figure()

        # Add significant points (p < 0.05 and |log2FC| > 1)
        significant = (results['p_value'] < 0.05) & (abs(results['log2_fold_change']) > 1)

        # Non-significant points
        fig.add_trace(go.Scatter(
            x=results[~significant]['log2_fold_change'],
            y=results[~significant]['neg_log10_pval'],
            mode='markers',
            name='Not Significant',
            marker=dict(color='gray', size=8),
            text=results[~significant]['gene'],
            hovertemplate='<b>%{text}</b><br>Log2 FC: %{x:.2f}<br>-log10(p): %{y:.2f}<extra></extra>'
        ))

        # Significant points
        fig.add_trace(go.Scatter(
            x=results[significant]['log2_fold_change'],
            y=results[significant]['neg_log10_pval'],
            mode='markers',
            name='Significant',
            marker=dict(color='red', size=10),
            text=results[significant]['gene'],
            hovertemplate='<b>%{text}</b><br>Log2 FC: %{x:.2f}<br>-log10(p): %{y:.2f}<extra></extra>'
        ))

        # Add threshold lines
        fig.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="gray", opacity=0.5)
        fig.add_vline(x=-1, line_dash="dash", line_color="gray", opacity=0.5)
        fig.add_vline(x=1, line_dash="dash", line_color="gray", opacity=0.5)

        fig.update_layout(
            title="Protein Tumor vs Normal - Volcano Plot",
            xaxis_title="Log2 Fold Change (Tumor - Normal)",
            yaxis_title="-Log10(p-value)",
            hovermode='closest',
            template='plotly_white',
            height=600
        )

        return fig

    # Phospho Analysis
    phospho_results = reactive.Value(None)

    @reactive.Effect
    @reactive.event(input.phospho_run)
    def run_phospho_analysis():
        cancer = input.phospho_cancer()
        query = input.phospho_query().strip()
        normalized = input.phospho_normalized()

        if not query:
            phospho_results.set(pd.DataFrame({'error': ['Please enter phosphosites or genes']}))
            return

        try:
            # Load cancer data
            if cancer == "brca":
                data = cptac.Brca()
            elif cancer == "coad":
                data = cptac.Coad()
            elif cancer == "hnscc":
                data = cptac.Hnscc()
            elif cancer == "luad":
                data = cptac.Luad()
            elif cancer == "ovarian":
                data = cptac.Ovarian()
            elif cancer == "ccrcc":
                data = cptac.Ccrcc()
            elif cancer == "gbm":
                data = cptac.Gbm()
            elif cancer == "lscc":
                data = cptac.Lscc()
            elif cancer == "pdac":
                data = cptac.Pdac()

            # Get phosphoproteomics data
            if normalized:
                phospho = data.get_phosphoproteomics(source='harmonized')
                proteomics = data.get_proteomics()

                # Normalize phospho by protein levels
                # This is a simplified normalization - subtract protein levels
                if isinstance(phospho.columns, pd.MultiIndex):
                    phospho = phospho.groupby(level=list(range(phospho.columns.nlevels)), axis=1).mean()
            else:
                phospho = data.get_phosphoproteomics(source='harmonized')
                if isinstance(phospho.columns, pd.MultiIndex):
                    phospho = phospho.groupby(level=list(range(phospho.columns.nlevels)), axis=1).mean()

            # Transpose
            phospho = phospho.T

            # Parse query
            queries = [q.strip() for q in query.split(',') if q.strip()]

            results = []
            for q in queries:
                # Check if it's a specific phosphosite (contains underscore)
                if '_' in q:
                    # Specific phosphosite query
                    matching = [idx for idx in phospho.index if q in str(idx)]
                else:
                    # Gene query - find all phosphosites for that gene
                    matching = [idx for idx in phospho.index if q in str(idx)]

                for site in matching[:25]:  # Limit to 25 sites per gene
                    site_data = phospho.loc[site]
                    tumor_samples, normal_samples = get_tumor_normal_pairs(pd.DataFrame(site_data).T)

                    if len(tumor_samples) > 0:
                        tumor_data = site_data[tumor_samples]
                        normal_data = site_data[[s for s in normal_samples if s[:-2] in tumor_samples]]

                        stats_result = calculate_fold_change_stats(tumor_data, normal_data)
                        if stats_result:
                            stats_result['site'] = str(site)
                            results.append(stats_result)

            if results:
                results_df = pd.DataFrame(results)
                results_df = apply_fdr_correction(results_df)
                phospho_results.set(results_df)
            else:
                phospho_results.set(pd.DataFrame({'error': ['No results found for the specified query']}))

        except Exception as e:
            phospho_results.set(pd.DataFrame({'error': [f'Error: {str(e)}']}))

    @output
    @render.ui
    def phospho_results_ui():
        results = phospho_results.get()
        if results is None:
            return ui.p("Run analysis to see results")

        if 'error' in results.columns:
            return ui.p(results['error'].iloc[0], style="color: red;")

        return ui.HTML(results.to_html(index=False, float_format=lambda x: f'{x:.4f}'))

    @output
    @render.plot
    def phospho_plot():
        results = phospho_results.get()
        if results is None or 'error' in results.columns:
            return None

        # Create volcano plot
        results['neg_log10_pval'] = -np.log10(results['p_value'])

        fig = go.Figure()

        # Add significant points
        significant = (results['p_value'] < 0.05) & (abs(results['log2_fold_change']) > 1)

        fig.add_trace(go.Scatter(
            x=results[~significant]['log2_fold_change'],
            y=results[~significant]['neg_log10_pval'],
            mode='markers',
            name='Not Significant',
            marker=dict(color='gray', size=8),
            text=results[~significant]['site'],
            hovertemplate='<b>%{text}</b><br>Log2 FC: %{x:.2f}<br>-log10(p): %{y:.2f}<extra></extra>'
        ))

        fig.add_trace(go.Scatter(
            x=results[significant]['log2_fold_change'],
            y=results[significant]['neg_log10_pval'],
            mode='markers',
            name='Significant',
            marker=dict(color='red', size=10),
            text=results[significant]['site'],
            hovertemplate='<b>%{text}</b><br>Log2 FC: %{x:.2f}<br>-log10(p): %{y:.2f}<extra></extra>'
        ))

        fig.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="gray", opacity=0.5)
        fig.add_vline(x=-1, line_dash="dash", line_color="gray", opacity=0.5)
        fig.add_vline(x=1, line_dash="dash", line_color="gray", opacity=0.5)

        fig.update_layout(
            title="Phospho Tumor vs Normal - Volcano Plot",
            xaxis_title="Log2 Fold Change (Tumor - Normal)",
            yaxis_title="-Log10(p-value)",
            hovermode='closest',
            template='plotly_white',
            height=600
        )

        return fig

    # Correlation Analysis
    corr_results = reactive.Value(None)
    corr_data = reactive.Value(None)  # Store raw data for plotting

    @reactive.Effect
    @reactive.event(input.corr_run)
    def run_correlation_analysis():
        cancer = input.corr_cancer()
        query = input.corr_query().strip()
        data_type = input.corr_data_type()
        normalized = input.corr_normalized()

        if not query:
            corr_results.set(pd.DataFrame({'error': ['Please enter genes or phosphosites']}))
            return

        try:
            # Load cancer data
            if cancer == "brca":
                data = cptac.Brca()
            elif cancer == "coad":
                data = cptac.Coad()
            elif cancer == "hnscc":
                data = cptac.Hnscc()
            elif cancer == "luad":
                data = cptac.Luad()
            elif cancer == "ovarian":
                data = cptac.Ovarian()
            elif cancer == "ccrcc":
                data = cptac.Ccrcc()
            elif cancer == "gbm":
                data = cptac.Gbm()
            elif cancer == "lscc":
                data = cptac.Lscc()
            elif cancer == "pdac":
                data = cptac.Pdac()

            # Parse query
            queries = [q.strip() for q in query.split(',') if q.strip()]

            # Collect data based on data type
            all_data = []
            item_labels = []

            if data_type in ['phospho', 'both']:
                phospho = data.get_phosphoproteomics(source='harmonized')
                if isinstance(phospho.columns, pd.MultiIndex):
                    phospho = phospho.groupby(level=list(range(phospho.columns.nlevels)), axis=1).mean()
                phospho = phospho.T

                for q in queries:
                    if '_' in q or data_type == 'phospho':
                        matching = [idx for idx in phospho.index if q in str(idx)]
                        for site in matching[:25]:
                            all_data.append(phospho.loc[site])
                            item_labels.append(str(site))

            if data_type in ['proteomics', 'both']:
                proteomics = data.get_proteomics()
                if isinstance(proteomics.columns, pd.MultiIndex):
                    proteomics = proteomics.groupby(level=list(range(proteomics.columns.nlevels)), axis=1).mean()
                proteomics = proteomics.T

                for q in queries:
                    if '_protein' in q or data_type == 'proteomics':
                        gene = q.replace('_protein', '')
                        if gene in proteomics.index:
                            all_data.append(proteomics.loc[gene])
                            item_labels.append(gene)

            if len(all_data) < 2:
                corr_results.set(pd.DataFrame({'error': ['Need at least 2 items for correlation analysis']}))
                return

            # Create data matrix
            data_matrix = pd.DataFrame({label: data for label, data in zip(item_labels, all_data)})

            # Filter to tumor samples only (exclude .N samples)
            tumor_cols = [col for col in data_matrix.index if not is_normal_sample(col)]
            data_matrix_tumor = data_matrix.loc[tumor_cols]

            # Calculate correlation
            corr_matrix = data_matrix_tumor.corr(method='pearson')

            # Store for plotting
            corr_data.set({
                'matrix': corr_matrix,
                'data': data_matrix,
                'labels': item_labels
            })

            corr_results.set(corr_matrix)

        except Exception as e:
            corr_results.set(pd.DataFrame({'error': [f'Error: {str(e)}']}))

    @output
    @render.ui
    def corr_results_ui():
        results = corr_results.get()
        if results is None:
            return ui.p("Run analysis to see results")

        if 'error' in results.columns:
            return ui.p(results['error'].iloc[0], style="color: red;")

        return ui.HTML(results.to_html(float_format=lambda x: f'{x:.3f}'))

    @output
    @render.plot
    def corr_plot():
        results = corr_results.get()
        data_dict = corr_data.get()

        if results is None or data_dict is None or 'error' in results.columns:
            return None

        # If we have 2 items, create scatter plot
        if len(data_dict['labels']) == 2:
            data_matrix = data_dict['data']
            label1, label2 = data_dict['labels']

            # Create scatter plot with tumor/normal coloring
            sample_types = ['Normal' if is_normal_sample(idx) else 'Tumor' for idx in data_matrix.index]

            fig = go.Figure()

            # Add tumor samples
            tumor_mask = [st == 'Tumor' for st in sample_types]
            fig.add_trace(go.Scatter(
                x=data_matrix.loc[tumor_mask, label1],
                y=data_matrix.loc[tumor_mask, label2],
                mode='markers',
                name='Tumor',
                marker=dict(color='red', size=10),
                text=data_matrix.index[tumor_mask],
                hovertemplate='<b>%{text}</b><br>%{xaxis.title.text}: %{x:.2f}<br>%{yaxis.title.text}: %{y:.2f}<extra></extra>'
            ))

            # Add normal samples
            normal_mask = [st == 'Normal' for st in sample_types]
            if any(normal_mask):
                fig.add_trace(go.Scatter(
                    x=data_matrix.loc[normal_mask, label1],
                    y=data_matrix.loc[normal_mask, label2],
                    mode='markers',
                    name='Normal',
                    marker=dict(color='blue', size=10),
                    text=data_matrix.index[normal_mask],
                    hovertemplate='<b>%{text}</b><br>%{xaxis.title.text}: %{x:.2f}<br>%{yaxis.title.text}: %{y:.2f}<extra></extra>'
                ))

            # Add correlation line
            valid_data = data_matrix[[label1, label2]].dropna()
            if len(valid_data) > 1:
                z = np.polyfit(valid_data[label1], valid_data[label2], 1)
                p = np.poly1d(z)
                x_line = np.linspace(valid_data[label1].min(), valid_data[label1].max(), 100)
                fig.add_trace(go.Scatter(
                    x=x_line,
                    y=p(x_line),
                    mode='lines',
                    name='Trend',
                    line=dict(color='gray', dash='dash'),
                    showlegend=True
                ))

            corr_val = results.loc[label1, label2]

            fig.update_layout(
                title=f"Correlation Plot (r = {corr_val:.3f})",
                xaxis_title=label1,
                yaxis_title=label2,
                hovermode='closest',
                template='plotly_white',
                height=600
            )

            return fig

        else:
            # Create heatmap for multiple items
            fig = go.Figure(data=go.Heatmap(
                z=results.values,
                x=results.columns,
                y=results.index,
                colorscale='RdBu_r',
                zmid=0,
                zmin=-1,
                zmax=1,
                text=results.values,
                texttemplate='%{text:.2f}',
                textfont={"size": 10},
                colorbar=dict(title="Correlation")
            ))

            fig.update_layout(
                title="Correlation Matrix",
                xaxis_title="",
                yaxis_title="",
                template='plotly_white',
                height=600,
                xaxis={'side': 'bottom'}
            )

            return fig


app = App(app_ui, server)
