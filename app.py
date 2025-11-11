"""
CPTAC Browser - Python Shiny GUI for CPTAC Data Analysis (Optimized)
"""

from shiny import App, ui, render, reactive
import cptac
import pandas as pd
import numpy as np
from scipy import stats
import plotly.graph_objects as go
import plotly.express as px
from typing import Tuple, List, Optional, Dict
from dataclasses import dataclass
import sys

# Debug mode - set to True to see detailed error messages
DEBUG = True

# Available cancer types
CANCER_TYPES = {
    "brca": "Breast Cancer",
    "coad": "Colon Adenocarcinoma",
    "hnscc": "Head and Neck Squamous Cell Carcinoma",
    "luad": "Lung Adenocarcinoma",
    "ov": "Ovarian Cancer",
    "ccrcc": "Clear Cell Renal Cell Carcinoma",
    "gbm": "Glioblastoma",
    "lscc": "Lung Squamous Cell Carcinoma",
    "pdac": "Pancreatic Ductal Adenocarcinoma"
}

@dataclass
class ProcessedData:
    """Container for pre-processed cancer data"""
    proteomics: pd.DataFrame
    phospho_normalized: pd.DataFrame
    phospho_unnormalized: pd.DataFrame
    tumor_samples: List[str]
    normal_samples: List[str]
    paired_tumor_samples: List[str]
    paired_normal_samples: List[str]

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
                ui.input_selectize(
                    "protein_genes",
                    "Select Genes:",
                    choices=[],
                    multiple=True,
                    options={
                        "placeholder": "Start typing to search genes (e.g., TP53, EGFR)...",
                        "maxItems": 50
                    }
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
                ui.input_switch(
                    "phospho_normalized",
                    "Normalized (protein-adjusted)",
                    value=True
                ),
                ui.input_selectize(
                    "phospho_query",
                    "Select Phosphosites:",
                    choices=[],
                    multiple=True,
                    options={
                        "placeholder": "Start typing gene or phosphosite (e.g., EGFR_Y1068, TP53)...",
                        "maxItems": 50
                    }
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
                ui.input_select(
                    "corr_data_type",
                    "Data Type:",
                    choices={"proteomics": "Proteomics", "phospho": "Phosphoproteomics", "both": "Both"}
                ),
                ui.input_switch(
                    "corr_normalized",
                    "Normalized Phospho (protein-adjusted)",
                    value=True
                ),
                ui.input_selectize(
                    "corr_query",
                    "Select Genes/Phosphosites:",
                    choices=[],
                    multiple=True,
                    options={
                        "placeholder": "Start typing to search...",
                        "maxItems": 50
                    }
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

    # ==================== Helper Functions ====================

    def is_normal_sample(sample_id: str) -> bool:
        """Check if sample ID indicates normal tissue (ends with .N)"""
        return str(sample_id).endswith('.N')

    def load_cancer_data(cancer: str):
        """Load cancer data object from cptac library"""
        cancer_map = {
            "brca": cptac.Brca,
            "coad": cptac.Coad,
            "hnscc": cptac.Hnscc,
            "luad": cptac.Luad,
            "ov": cptac.Ov,
            "ccrcc": cptac.Ccrcc,
            "gbm": cptac.Gbm,
            "lscc": cptac.Lscc,
            "pdac": cptac.Pdac
        }
        return cancer_map[cancer]()

    def process_dataframe(df: pd.DataFrame) -> pd.DataFrame:
        """Handle multi-index and transpose dataframe"""
        if DEBUG:
            print(f"[DEBUG] process_dataframe input type: {type(df)}, shape: {df.shape if hasattr(df, 'shape') else 'unknown'}", file=sys.stderr)

        if isinstance(df.columns, pd.MultiIndex):
            if DEBUG:
                print(f"[DEBUG] MultiIndex detected with {df.columns.nlevels} levels", file=sys.stderr)
            # Deduplicate multi-index columns by taking mean
            # Use column names to group, then aggregate
            level_names = list(range(df.columns.nlevels))
            df = df.groupby(level=level_names, axis=1).mean()
            if DEBUG:
                print(f"[DEBUG] After groupby type: {type(df)}", file=sys.stderr)

        # Ensure we have a proper DataFrame after groupby
        if not isinstance(df, pd.DataFrame):
            if DEBUG:
                print(f"[DEBUG] Converting {type(df)} to DataFrame", file=sys.stderr)
            df = pd.DataFrame(df)

        # Transpose so features are rows and samples are columns
        result = df.T
        if DEBUG:
            print(f"[DEBUG] process_dataframe output type: {type(result)}, shape: {result.shape}", file=sys.stderr)
        return result

    def load_phosphoproteomics(data) -> Optional[pd.DataFrame]:
        """Load phosphoproteomics data with fallback to multiple sources"""
        sources = ['harmonized', 'broad', 'mssm', 'washu', None]

        for source in sources:
            try:
                if source:
                    phospho_raw = data.get_phosphoproteomics(source=source)
                else:
                    phospho_raw = data.get_phosphoproteomics()

                if phospho_raw is not None and not phospho_raw.empty:
                    return phospho_raw
            except (AttributeError, ValueError, KeyError, Exception):
                continue

        return None

    def get_tumor_normal_pairs(columns) -> Tuple[List[str], List[str], List[str], List[str]]:
        """
        Extract all tumor/normal samples and paired samples.
        Returns: (all_tumor, all_normal, paired_tumor, paired_normal)
        """
        # Ensure columns is a proper list
        if not isinstance(columns, list):
            columns = list(columns)

        all_normal = [col for col in columns if is_normal_sample(col)]
        all_tumor = [col for col in columns if not is_normal_sample(col)]

        paired_tumor = []
        paired_normal = []

        for normal in all_normal:
            tumor = normal[:-2]  # Remove '.N' suffix
            if tumor in all_tumor:
                paired_tumor.append(tumor)
                paired_normal.append(normal)

        return all_tumor, all_normal, paired_tumor, paired_normal

    def calculate_fold_change_stats(tumor_data: pd.Series, normal_data: pd.Series) -> Optional[dict]:
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

    # ==================== Cached Data Loaders ====================

    @reactive.calc
    def protein_data() -> Optional[ProcessedData]:
        """Load and cache protein analysis data for selected cancer type"""
        cancer = input.protein_cancer()

        try:
            data = load_cancer_data(cancer)

            # Process proteomics data
            proteomics = process_dataframe(data.get_proteomics())

            # Get tumor/normal pairs
            all_tumor, all_normal, paired_tumor, paired_normal = get_tumor_normal_pairs(
                proteomics.columns
            )

            # For protein analysis, we don't need phospho data
            return ProcessedData(
                proteomics=proteomics,
                phospho_normalized=pd.DataFrame(),
                phospho_unnormalized=pd.DataFrame(),
                tumor_samples=all_tumor,
                normal_samples=all_normal,
                paired_tumor_samples=paired_tumor,
                paired_normal_samples=paired_normal
            )
        except Exception as e:
            print(f"Error loading protein data: {e}", file=sys.stderr)
            if DEBUG:
                import traceback
                traceback.print_exc()
            return None

    @reactive.calc
    def phospho_data() -> Optional[ProcessedData]:
        """Load and cache phospho analysis data for selected cancer type"""
        cancer = input.phospho_cancer()
        normalized = input.phospho_normalized()

        try:
            data = load_cancer_data(cancer)

            # Load phospho data with fallback sources
            phospho_raw = load_phosphoproteomics(data)

            if phospho_raw is None:
                print(f"No phosphoproteomics data available for {cancer}")
                return None

            phospho_processed = process_dataframe(phospho_raw)

            # Store based on normalization preference
            if normalized:
                phospho_norm = phospho_processed
                phospho_unnorm = pd.DataFrame()
            else:
                phospho_norm = pd.DataFrame()
                phospho_unnorm = phospho_processed

            # Get tumor/normal pairs
            all_tumor, all_normal, paired_tumor, paired_normal = get_tumor_normal_pairs(
                phospho_processed.columns
            )

            return ProcessedData(
                proteomics=pd.DataFrame(),
                phospho_normalized=phospho_norm,
                phospho_unnormalized=phospho_unnorm,
                tumor_samples=all_tumor,
                normal_samples=all_normal,
                paired_tumor_samples=paired_tumor,
                paired_normal_samples=paired_normal
            )
        except Exception as e:
            print(f"Error loading phospho data: {e}", file=sys.stderr)
            if DEBUG:
                import traceback
                traceback.print_exc()
            return None

    @reactive.calc
    def correlation_data() -> Optional[Dict]:
        """Load and cache correlation analysis data for selected cancer type"""
        cancer = input.corr_cancer()
        data_type = input.corr_data_type()
        normalized = input.corr_normalized()

        try:
            data = load_cancer_data(cancer)

            result = {}

            # Load proteomics if needed
            if data_type in ['proteomics', 'both']:
                result['proteomics'] = process_dataframe(data.get_proteomics())

            # Load phospho if needed
            if data_type in ['phospho', 'both']:
                phospho_raw = load_phosphoproteomics(data)

                if phospho_raw is not None:
                    result['phospho'] = process_dataframe(phospho_raw)

            return result
        except Exception as e:
            print(f"Error loading correlation data: {e}", file=sys.stderr)
            if DEBUG:
                import traceback
                traceback.print_exc()
            return None

    # ==================== Update Selectize Choices ====================

    @reactive.Effect
    def update_protein_gene_choices():
        """Update gene choices when cancer type changes"""
        data = protein_data()
        if data is not None and not data.proteomics.empty:
            genes = sorted([str(idx) for idx in data.proteomics.index])
            ui.update_selectize("protein_genes", choices=genes)

    @reactive.Effect
    def update_phospho_choices():
        """Update phosphosite choices when cancer type or normalization changes"""
        data = phospho_data()
        if data is not None:
            # Get the appropriate phospho data
            normalized = input.phospho_normalized()
            phospho = data.phospho_normalized if normalized else data.phospho_unnormalized
            if phospho.empty:
                phospho = data.phospho_normalized if not normalized else data.phospho_unnormalized

            if not phospho.empty:
                # Format phosphosites as readable strings
                phosphosites = sorted([str(idx) for idx in phospho.index])
                ui.update_selectize("phospho_query", choices=phosphosites)
            else:
                ui.update_selectize("phospho_query", choices=[])
        else:
            # No phospho data available for this cancer type
            ui.update_selectize("phospho_query", choices=[])

    @reactive.Effect
    def update_correlation_choices():
        """Update choices when cancer type or data type changes"""
        data = correlation_data()
        if data is not None:
            choices = []
            data_type = input.corr_data_type()

            if data_type in ['proteomics', 'both'] and 'proteomics' in data:
                genes = [str(idx) for idx in data['proteomics'].index]
                choices.extend(genes)

            if data_type in ['phospho', 'both'] and 'phospho' in data:
                phosphosites = [str(idx) for idx in data['phospho'].index]
                choices.extend(phosphosites)

            ui.update_selectize("corr_query", choices=sorted(set(choices)))

    # ==================== Protein Analysis ====================

    protein_results = reactive.Value(None)

    @reactive.Effect
    @reactive.event(input.protein_run)
    def run_protein_analysis():
        genes = input.protein_genes()

        # Convert to list if needed
        if isinstance(genes, str):
            genes = [g.strip() for g in genes.split(',') if g.strip()]
        elif not genes:
            genes = []

        if not genes:
            protein_results.set(pd.DataFrame({'error': ['Please enter gene names']}))
            return

        try:
            data = protein_data()
            if data is None:
                protein_results.set(pd.DataFrame({'error': ['Failed to load cancer data']}))
                return

            proteomics = data.proteomics
            paired_tumor = data.paired_tumor_samples
            paired_normal = data.paired_normal_samples

            results = []
            for gene in genes:
                if gene in proteomics.index:
                    gene_data = proteomics.loc[gene]

                    if len(paired_tumor) > 0:
                        tumor_data = gene_data[paired_tumor]
                        normal_data = gene_data[paired_normal]

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

    # ==================== Phospho Analysis ====================

    phospho_results = reactive.Value(None)

    @reactive.Effect
    @reactive.event(input.phospho_run)
    def run_phospho_analysis():
        query = input.phospho_query()
        normalized = input.phospho_normalized()

        # Convert to list if needed
        if isinstance(query, str):
            queries = [q.strip() for q in query.split(',') if q.strip()]
        else:
            queries = list(query) if query else []

        if not queries:
            phospho_results.set(pd.DataFrame({'error': ['Please select phosphosites or genes']}))
            return

        try:
            data = phospho_data()
            if data is None:
                cancer_name = CANCER_TYPES.get(input.phospho_cancer(), input.phospho_cancer())
                phospho_results.set(pd.DataFrame({'error': [f'No phosphoproteomics data available for {cancer_name}. Try a different cancer type.']}))
                return

            # Get the appropriate phospho data
            phospho = data.phospho_normalized if normalized else data.phospho_unnormalized
            if phospho.empty:
                phospho = data.phospho_normalized if not normalized else data.phospho_unnormalized

            paired_tumor = data.paired_tumor_samples
            paired_normal = data.paired_normal_samples

            results = []
            for q in queries:
                # Try direct lookup first (from selectize selection)
                if q in phospho.index:
                    sites_to_process = [q]
                else:
                    # Fallback to substring matching for manual entry
                    sites_to_process = [idx for idx in phospho.index if q in str(idx)][:25]

                for site in sites_to_process:
                    site_data = phospho.loc[site]

                    if len(paired_tumor) > 0:
                        tumor_data = site_data[paired_tumor]
                        normal_data = site_data[paired_normal]

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

    # ==================== Correlation Analysis ====================

    corr_results = reactive.Value(None)
    corr_data_store = reactive.Value(None)

    @reactive.Effect
    @reactive.event(input.corr_run)
    def run_correlation_analysis():
        query = input.corr_query()
        data_type = input.corr_data_type()

        # Convert to list if needed
        if isinstance(query, str):
            queries = [q.strip() for q in query.split(',') if q.strip()]
        else:
            queries = list(query) if query else []

        if not queries:
            corr_results.set(pd.DataFrame({'error': ['Please select genes or phosphosites']}))
            return

        try:
            data = correlation_data()
            if data is None:
                corr_results.set(pd.DataFrame({'error': ['Failed to load cancer data']}))
                return

            # Collect data
            all_data = []
            item_labels = []

            if data_type in ['phospho', 'both'] and 'phospho' in data:
                phospho = data['phospho']
                for q in queries:
                    # Try direct lookup first (from selectize selection)
                    if q in phospho.index:
                        all_data.append(phospho.loc[q])
                        item_labels.append(str(q))
                    elif '_' in q or data_type == 'phospho':
                        # Fallback to substring matching for manual entry
                        matching = [idx for idx in phospho.index if q in str(idx)]
                        for site in matching[:25]:
                            all_data.append(phospho.loc[site])
                            item_labels.append(str(site))

            if data_type in ['proteomics', 'both'] and 'proteomics' in data:
                proteomics = data['proteomics']
                for q in queries:
                    # Try direct lookup first
                    if q in proteomics.index:
                        all_data.append(proteomics.loc[q])
                        item_labels.append(str(q))
                    elif '_protein' in q or data_type == 'proteomics':
                        # Fallback to substring matching for manual entry
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
            corr_data_store.set({
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
        data_dict = corr_data_store.get()

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
