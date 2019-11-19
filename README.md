# Python Fastjet Bindings

## Installation

Prerequisites:

0. Clone the repository, and initialize the submodule with `git submodule init`.

1. You need to set the environment variables:

    - `FASTJET_ROOT`
    - `CGAL_ROOT`
    - `GMP_ROOT`

    They need to be set to the install locations of the three packages. `alibuild` will set the first two.
    However, `GMP_ROOT` unfortunately needs to be set by hand...

2. Install `poetry`, and install the package via `poetry install` in the repository root directory.

3. The tests can be via `poetry run pytest -l -vv tests/`.

Done. If you need a virtualenv with the package, this can be useful:

```bash
alias poetryShell='source "$(dirname $(poetry run which python))/activate"'
```

## Quick Start

See the tests.

# Alternatives

I was looking for very simple bindings, and control of the API was helpful, so I decided to write my own. But
there are many other excellent binding options:

- Official bindings: Fastjet recently pushed SWIG based python bindings. You need a fairly recent version of
  fastjet, which isn't yet deployed for ALICE.
- [pyjet](https://github.com/scikit-hep/pyjet): It currently (Nov 2019) doesn't have fastjet-contrib bindings.
  They possible could have been added, but it didn't appear as straightforward as I would have liked.
- [heppy](https://github.com/matplo/heppy): I was unable to get this setup :-(


