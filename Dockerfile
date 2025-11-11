# CPTAC Browser Dockerfile
# Containerized Python Shiny application for CPTAC data analysis

FROM python:3.11-slim

# Set metadata
LABEL maintainer="CPTAC Browser"
LABEL description="Python Shiny GUI for browsing and analyzing CPTAC proteomics data"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy requirements first for better Docker layer caching
COPY requirements.txt .

# Install Python dependencies
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

# Copy application files
COPY app.py .
COPY README.md .

# Create directory for CPTAC data cache
RUN mkdir -p /root/.cptac

# Expose Shiny default port
EXPOSE 8000

# Set environment variables
ENV PYTHONUNBUFFERED=1
ENV SHINY_HOST=0.0.0.0
ENV SHINY_PORT=8000

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=60s --retries=3 \
    CMD curl -f http://localhost:8000/ || exit 1

# Run the Shiny application
CMD ["shiny", "run", "app.py", "--host", "0.0.0.0", "--port", "8000"]
