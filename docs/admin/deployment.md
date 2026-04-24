# Deployment

## ShinyApps.io

The first public deployment target is ShinyApps.io.

Before publishing:

```bash
Rscript scripts/check_bundle_size.R
Rscript scripts/release_check.R --skip-browser
```

Set deployment environment variables in the ShinyApps.io dashboard as needed:

```text
ATLAS_ACCESS_PASSWORD
ATLAS_VERSION
ATLAS_LAST_UPDATED
ATLAS_CITATION
```

Set `ATLAS_ACCESS_PASSWORD` when you need a private, workshop, or restricted-access deployment.

## Bundle Size

The full atlas bundle is close to the practical 5 GB threshold. `scripts/check_bundle_size.R` estimates the size of the app, metadata, annotations, orthogroups, and data folders.

If future releases exceed the limit, use versioned public-read external artifacts for large data files and keep the app interface unchanged.

## Docker Or VM Hosting

Docker, Posit Connect, or a VM with read-only mounted data is the safest deployment route when memory, bundle size, or startup time become limiting.

Docker quick start:

```bash
docker compose up --build
```

The Compose file mounts data, metadata, annotation, and orthology folders read-only by default.
