# Troubleshooting

## Command-line values do not match the configured defaults

**Symptom:** A parameter uses a different value than expected for the selected
use case.

**Cause:** Each sub-command has its own set of default values tailored for the
associated use case. A value set explicitly on the command line, or in a
`--config` JSON file, takes precedence over the use case default.

**Resolution:** Command-line values always win. When `--config` is supplied,
use case-specific defaults are taken from the JSON file, but any value passed on
the command line overwrites the corresponding config file value. Set the
parameter explicitly on the command line to force a specific value.

## Unexpected variants included or excluded at a threshold boundary

**Symptom:** A variant at the exact threshold value is included when you
expected it to be excluded (or vice versa).

**Cause:** Unless otherwise specified, all threshold cutoffs and minimum values
used for variant filtering and feature computation are inclusive.

**Resolution:** Account for inclusive boundaries when choosing thresholds. For
example, a minimum base quality of 6 includes a base with quality 6 but excludes
one with quality 5.

If you have issues or concerns not documented here, please create an issue on [GitHub](https://github.com/Roche-AXELIOS/XOOS/issues) or reach out to your local Roche affiliate.
