# Contributing

## Setup

```console
$ git clone https://github.com/karine-mue/Silvispect.git
$ cd Silvispect
$ pip install -e ".[dev]"
$ make check          # ruff check, ruff format --check, pytest
```

Silvispect targets Python 3.10+ and has **no runtime dependencies**. That is a
design constraint, not an accident: the tool must run in a field laptop's
system Python without a compiler or a package index. Pull requests adding a
runtime dependency need to argue the case in the description.

## Principles

**Absent data is reported, never imputed.** A missing diameter stays `None` all
the way to a finding. Nodata cells stay `None`, never `-9999`, so a sentinel
can never enter a mean.

**Every number is attributable.** If a metric assumes a form factor, a wood
density or an expansion area, that assumption is a named parameter with a
documented default — not a constant buried in an expression.

**Determinism.** The same input produces the same output, byte for byte.
Anything random is seeded, and ties are broken explicitly.

**Rules are cheap to add and stable forever.** Codes never change meaning; a
changed rule takes a new code.

## Tests

```console
$ pytest                        # the whole suite
$ pytest tests/test_detect.py   # one module
$ pytest -k gap                 # by name
```

Tests are the specification, so a change in behaviour should show up as a
changed test. Some conventions:

- **Test the algorithm against ground truth, not against itself.** Detection
  accuracy is asserted against a synthetic stand whose trees are known
  exactly (`tests/conftest.py`), and curve fitting is checked by recovering the
  parameters the generator used. Do not assert on a value the implementation
  happened to produce.
- **Name the behaviour, not the function.**
  `test_plateau_yields_a_single_top` beats `test_find_treetops_2`.
- **Cover the failure path.** Every `raise` should have a test that provokes it.
- The synthetic stand fixtures are session-scoped; reuse them rather than
  generating new stands in each test.

## Adding an inspection rule

1. Write the rule in `silvispect/inspection.py`, decorated with `@rule`, taking
   an `InspectionContext` and yielding `Finding` objects. Return quietly when
   its inputs are absent — check `ctx.has_field_data` / `ctx.has_raster`.
2. Add any threshold as a field of `InspectionConfig` with a default and a
   docstring entry.
3. Document it in [`docs/findings.md`](docs/findings.md), including *why* the
   threshold has the value it does and what a user should do about a hit.
4. Add a test that fires the rule and one that confirms it stays quiet on clean
   data. The second matters more: a rule that always fires is noise.
5. Add the code to the table in the README.

A worked example is at the end of `docs/findings.md`.

## Style

`ruff` is the arbiter for both linting and formatting; `make format` applies
it. Beyond that:

- Public functions and classes get docstrings with an `Args:` / `Returns:` /
  `Raises:` section where there is something to say.
- Comments explain *why*, not *what*. The comment worth writing is the one
  recording the reason a naive implementation fails — see the morphological
  opening in `canopy.py`.
- Type annotations everywhere, `from __future__ import annotations` at the top.
- Metric units are the convention throughout: metres, centimetres for DBH,
  square metres, hectares.

## Sample data

`data/` holds a committed synthetic plot that the README, the docs and the CI
smoke test all quote. If you need to change it, regenerate rather than edit:

```console
$ make data
```

The generator is seeded, so the output is byte-identical unless the generator
itself changed — CI diffs the result to enforce that. If a change to `synth.py`
legitimately alters the sample, regenerate it in the same commit and update any
figures quoted in the README.
