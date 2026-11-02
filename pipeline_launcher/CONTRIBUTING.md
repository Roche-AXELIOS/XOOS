# Contributing

## Development Setup

Install the package in editable mode with dev dependencies:

```bash
pip install -e ".[dev]"
```

## Testing

Run the full test suite:

```bash
pytest
```

Run with verbose output:

```bash
pytest -v
```

## Code Style

Format and lint all Python files with ruff:

```bash
ruff format src/ tests/
ruff check src/ tests/
```

Run type checking:

```bash
mypy src/
```

## Project Layout

- `src/pipeline_launcher/` — package source
- `tests/` — pytest test suite, mirroring the source layout
- `src/pipeline_launcher/env/` — bundled environment config JSONs
- `src/pipeline_launcher/templates/` — bundled Jinja2 templates (package data)

## Adding a New Executor

1. Create `src/pipeline_launcher/executor/<name>.py` with a class implementing the `Executor` protocol (see `executor/base.py`).
2. Add a corresponding config dataclass in `config/model.py` inheriting from `BaseConfig`.
3. Update `config/loader.py` to parse the new executor type from JSON.
4. Register the new executor in the factory function in `executor/base.py`.
5. Add tests under `tests/executor/test_<name>.py`.

## Adding a New Environment Config

Drop a YAML file into `src/pipeline_launcher/env/`. The file stem becomes the `--env` name. See the README for the config schema.

## Env Config Overrides

Scalar Pydantic model fields in the config model (`str`, `int`, `float`, `bool`, and their `Optional` variants) are automatically overridable from the CLI via `--env-override key=value`. No code changes are needed when adding new scalar fields to config models — they become overridable immediately.

The override mechanism lives in `config/loader.py` (`apply_env_overrides`). It walks the model hierarchy using dot-notation keys and coerces string values to the target field's type. List and dict fields are intentionally excluded.
