# CPTAC Browser

Python Shiny GUI for browsing and analyzing CPTAC proteomics data. This application provides interactive tools for tumor vs. normal comparisons and correlation analysis across multiple cancer types.

## Features

### 1. Protein Tumor vs. Normal Analysis
- Compare protein expression levels between tumor and normal samples
- Performs paired t-tests with FDR correction
- Interactive volcano plots showing log2 fold change vs. significance
- Supports multiple genes simultaneously

### 2. Phosphoproteomics Tumor vs. Normal Analysis
- Analyze phosphorylation site changes between tumor and normal samples
- Option for normalized (protein-adjusted) or unnormalized phospho data
- Query by specific phosphosites (e.g., EGFR_Y1068) or gene names
- Interactive volcano plots with significance highlighting

### 3. Correlation Analysis
- Calculate correlations between proteins and/or phosphosites
- Scatter plots colored by tumor (red) vs. normal (blue) samples
- Correlation heatmaps for multiple features
- Supports phosphoproteomics, proteomics, or both data types

## Supported Cancer Types

- Breast Cancer (brca)
- Colon Adenocarcinoma (coad)
- Head and Neck Squamous Cell Carcinoma (hnscc)
- Lung Adenocarcinoma (luad)
- Ovarian Cancer (ovarian)
- Clear Cell Renal Cell Carcinoma (ccrcc)
- Glioblastoma (gbm)
- Lung Squamous Cell Carcinoma (lscc)
- Pancreatic Ductal Adenocarcinoma (pdac)

## Installation

```bash
pip install -r requirements.txt
```

## Running the Application

```bash
shiny run app.py
```

Then open your browser to the URL shown (typically `http://127.0.0.1:8000`)

## Usage

### Protein Analysis
1. Select a cancer type from the dropdown
2. Enter gene names (comma-separated), e.g., `TP53, EGFR, KRAS`
3. Click "Run Analysis"
4. View results table and interactive volcano plot

### Phospho Analysis
1. Select a cancer type
2. Enter phosphosites or genes (comma-separated)
   - Specific sites: `EGFR_Y1068, AKT1_S473`
   - All sites for a gene: `TP53, EGFR`
3. Toggle normalized/unnormalized option
4. Click "Run Analysis"

### Correlation Analysis
1. Select a cancer type
2. Enter genes/phosphosites to correlate (minimum 2 items)
3. Select data type (phospho, proteomics, or both)
4. For 2 items: scatter plot with tumor/normal coloring
5. For multiple items: correlation heatmap

## Data Source

This application uses the [CPTAC Python library](https://github.com/PayneLab/cptac) to access Clinical Proteomic Tumor Analysis Consortium data.

## Based on ProteomicsMCP

This GUI recreates the analysis tools from [ProteomicsMCP](https://github.com/plutzer/ProteomicsMCP) in an interactive web interface with Plotly visualizations.
