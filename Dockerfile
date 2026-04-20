# Legume Root Nodule Symbiosis Atlas — Docker image
#
# Build:
#   docker build -t legume-atlas .
#
# Run (mount your local data + orthologs + celltype overrides):
#   The within-species and cross-species `*_app_slim.rds` files live inside the
#   mounted `within_species_integrated_datasets/` and `app_ready_integration/`
#   folders below, so no extra volume mounts are required after you generate them.
#   docker run --rm -p 3838:3838 \
#       -v "$(pwd)/within_species_integrated_datasets:/srv/atlas/within_species_integrated_datasets:ro" \
#       -v "$(pwd)/app_ready_integration:/srv/atlas/app_ready_integration:ro" \
#       -v "$(pwd)/orthogroups:/srv/atlas/orthogroups:ro" \
#       -v "$(pwd)/annotations:/srv/atlas/annotations:ro" \
#       -v "$(pwd)/metadata:/srv/atlas/metadata:ro" \
#       -v "$(pwd)/celltype_overrides:/srv/atlas/celltype_overrides:ro" \
#       legume-atlas
#
# Override the default orthologs or celltype folder via env vars:
#   -e ATLAS_ORTHOLOGS_PATH=/srv/atlas/orthogroups/my_custom.tsv
#   -e ATLAS_CELLTYPE_OVERRIDES_DIR=/srv/atlas/celltype_overrides
#   -e ATLAS_ACCESS_PASSWORD=changeme     (leave unset locally to disable gate)
#   -e ATLAS_VERSION=0.1.0
#   -e ATLAS_LAST_UPDATED=2026-04-17
# CRAN snapshot pinned to 2026-04-01 for reproducible Docker builds.

FROM rocker/shiny:4.3.2

ENV DEBIAN_FRONTEND=noninteractive \
    R_LIBS_USER=/usr/local/lib/R/site-library \
    ATLAS_APP_DIR=/srv/atlas

RUN apt-get update && apt-get install -y --no-install-recommends \
        curl \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
        libglpk-dev \
        libhdf5-dev \
        libgit2-dev \
        cmake \
        git \
        libudunits2-dev \
        libgdal-dev \
        libgeos-dev \
        libproj-dev \
 && rm -rf /var/lib/apt/lists/*

RUN R -q -e "options(Ncpus = max(1L, parallel::detectCores() - 1L)); \
    install.packages(c( \
        'shiny', 'shinyWidgets', 'DT', 'dplyr', 'tidyr', 'purrr', 'tibble', \
        'ggplot2', 'patchwork', 'svglite', 'plotly', 'uwot', 'viridisLite', \
        'Seurat', 'scCustomize' \
    ), repos = 'https://packagemanager.posit.co/cran/2026-04-01')"

RUN R -q -e "if (!requireNamespace('BiocManager', quietly = TRUE)) \
        install.packages('BiocManager', repos = 'https://packagemanager.posit.co/cran/2026-04-01'); \
    BiocManager::install(c('ComplexHeatmap'), update = FALSE, ask = FALSE)"

WORKDIR ${ATLAS_APP_DIR}

COPY app.R ./app.R
COPY www ./www
COPY scripts ./scripts

# Default data folders. Mount over these with --volume at run time.
RUN mkdir -p within_species_integrated_datasets \
             app_ready_integration \
             orthogroups \
             annotations \
             metadata \
             celltype_overrides

EXPOSE 3838

HEALTHCHECK --interval=30s --timeout=5s --start-period=60s --retries=3 \
    CMD curl -fsS http://localhost:3838/ || exit 1

CMD ["R", "-e", "shiny::runApp('/srv/atlas', host = '0.0.0.0', port = 3838)"]
