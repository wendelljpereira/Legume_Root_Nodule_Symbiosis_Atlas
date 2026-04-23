# Quick Start

## Open The Hosted App

Open the atlas URL provided by the project team. During the pre-publication period, the hosted ShinyApps.io deployment may require a password.

After the app loads, start with the **Overview** page to confirm the atlas version, data-update date, and available species/integration summaries.

## Run The App Locally With Docker

From the repository root:

```bash
docker compose up --build
```

Open <http://localhost:3838>.

Docker is the recommended local route because it keeps R package setup and app execution close to the release environment.

## Run Directly In R

If R and the package dependencies are installed:

```r
shiny::runApp()
```

The app expects the data, metadata, orthology, and annotation folders to be present in the repository layout.

## First Three Things To Try

1. Open **Medicago truncatula** and choose a gene from the gene picker.
2. Click **Generate the expression plots** to render the expression panels.
3. Open **Camex** or **SATURN** and repeat the workflow to inspect ortholog-aware mapping.

```{note}
Adding genes to the selector does not immediately redraw plots. This is intentional: the app separates a staged gene panel from an applied gene panel so users can edit a selection before triggering heavier plotting work.
```
