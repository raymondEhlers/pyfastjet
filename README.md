# Python Fastjet Bindings

**NOTE: This was initial development work that was later in moved into other packages. For current options, see [scikit-hep/fastjet](https://github.com/scikit-hep/fastjet), or for my full analysis framework including jet finding, see [raymondEhlers/mammoth](https://github.com/raymondEhlers/mammoth).

Python bindings for fastjet. fastjet-contrib is not yet supported, but could be added depending on demand.

## Installation

Prerequisites:

0. Clone the repository, and initialize the submodule with `git submodule init`.

1. You need to set the environment variables:

    - `FASTJET`
    - `CGAL_ROOT`
    - `GMP_ROOT`

    They need to be set to the install locations of the three packages. `alibuild` will set the first two.
    (Note that alibuild used to call it `$FASTJET_ROOT`, but now it's known as just `$FASTJET`).
    However, `GMP_ROOT` unfortunately isn't seen as a fastjet dependency, so we have to load it explicitly.
    For convenience with AliBuild, you can use (assuming zsh):

    ```bash
    $ alienv modulecmd zsh load fastjet/latest GMP/latest
    ```

2. Install `poetry`, and install the package via `poetry install` in the repository root directory.

3. The tests can be via `poetry run pytest -l -vv tests/`.

Done. If you need a virtualenv with the package, this can be useful:

```bash
alias poetryShell='source "$(dirname $(poetry run which python))/activate"'
```

## Usage

By way of an example from the fastjet quick start:

```python
import pyfastjet as fj

# Create the PseudoJets
# an event with three particles:   px    py  pz      E
particles = []
particles.append(PseudoJet( 99.0,  0.1,  0, 100.0))
particles.append(PseudoJet(  4.0, -0.1,  0,   5.0))
particles.append(PseudoJet(-99.0,    0,  0,  99.0))

# choose a jet definition
R = 0.7
jet_def = fj.JetDefinition(fj.JetAlgorithm.antikt_algorithm, R)

# run the clustering, extract the jets
cs = fj.ClusterSequence(particles, jet_def)
jets = fj.sorted_by_pt(cs.inclusive_jets())

for jet in cs:
    print(f"jet: {jet}")
    for j, constituent in enumerate(jet.constituents):
        print(f"    constituent {j}: {constituent}")
```

See the tests for more functionality.

# Alternatives

I was looking for very simple bindings, and control of the API was helpful, so I decided to write my own. But
there are many other excellent binding options:

- Official bindings: Fastjet recently pushed SWIG based python bindings. You need a fairly recent version of
  fastjet, which isn't yet deployed for ALICE.
- [pyjet](https://github.com/scikit-hep/pyjet): It currently (Nov 2019) doesn't have fastjet-contrib bindings.
  They possible could have been added, but it didn't appear as straightforward as I would have liked.
- [heppy](https://github.com/matplo/heppy): I was unable to get this setup :-(


