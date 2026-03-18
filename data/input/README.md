# Input Data Policy

Only small, intentional, version-worthy inputs should live in `data/input`.

Do not commit:

- raw datasets
- large exported simulations
- cached preprocessing outputs
- private or credentialed data

Use these ignored folders for local-only files if needed:

- `data/input/local/`
- `data/input/raw/`
- `data/input/cache/`

If a dataset is too large for normal git hosting, keep it outside the repository and
document how to obtain or regenerate it.
