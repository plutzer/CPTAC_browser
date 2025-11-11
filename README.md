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

### Option 1: Docker (Recommended)

Docker provides an isolated, reproducible environment for running the application.

#### Build the Docker Image

```bash
docker build -t cptac-browser .
```

#### Run the Container

**Basic usage:**
```bash
docker run -p 8000:8000 cptac-browser
```

**With persistent data cache (recommended):**
```bash
docker run -p 8000:8000 -v cptac-data:/root/.cptac cptac-browser
```

This creates a named volume `cptac-data` that persists CPTAC datasets between container restarts, significantly speeding up subsequent launches.

**Run in background (detached mode):**
```bash
docker run -d -p 8000:8000 -v cptac-data:/root/.cptac --name cptac-browser cptac-browser
```

**Stop the container:**
```bash
docker stop cptac-browser
```

**View logs:**
```bash
docker logs cptac-browser
```

**Remove the container:**
```bash
docker rm cptac-browser
```

Then open your browser to `http://localhost:8000`

#### Docker Compose (Alternative)

For easier management, use Docker Compose:

```bash
# Start the application
docker-compose up -d

# View logs
docker-compose logs -f

# Stop the application
docker-compose down

# Stop and remove volumes (clears cache)
docker-compose down -v
```

### Option 2: Local Python Installation

```bash
pip install -r requirements.txt
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

## Docker Tips

### Persistent Data Storage

The CPTAC library downloads large dataset files on first use. To avoid re-downloading:

1. **Named volume (recommended):** `-v cptac-data:/root/.cptac`
2. **Bind mount to local directory:**
   ```bash
   mkdir -p ./cptac_cache
   docker run -p 8000:8000 -v ./cptac_cache:/root/.cptac cptac-browser
   ```

### Resource Allocation

For large datasets, you may want to increase Docker's memory allocation:

```bash
docker run -p 8000:8000 -m 4g -v cptac-data:/root/.cptac cptac-browser
```

### Troubleshooting

**Port already in use:**
```bash
# Use a different port
docker run -p 8080:8000 cptac-browser
# Then access at http://localhost:8080
```

**Container won't start:**
```bash
# Check logs
docker logs cptac-browser

# Run interactively to see errors
docker run -it -p 8000:8000 cptac-browser
```

**Clear cache and restart:**
```bash
docker-compose down -v
docker-compose up -d
```
