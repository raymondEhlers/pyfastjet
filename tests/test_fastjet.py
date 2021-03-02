#!/usr/bin/env python3

""" Basic tests for the fastjet binding.

.. code-author:: Raymond Ehlers <raymond.ehlers@cern.ch>, ORNL
"""

from typing import List

import awkward as ak
import pandas as pd
import pytest
import numpy as np

import pyfastjet as fj


@pytest.mark.parametrize("list_type", [
    "list", "numpy", "pandas"
], ids = ["PseudoJets", "Numpy arrays", "Pandas DataFrame"])
def test_fastjet(list_type: str) -> None:
    parts: List[List[float]] = []
    # an event with three particles:   px    py  pz      E
    parts.append([ 99.0,  0.1,  0, 100.0])
    parts.append([  4.0, -0.1,  0,   5.0])
    parts.append([-99.0,    0,  0,  99.0])

    # Construct the PseudoJets
    particles_pseudojets = [fj.PseudoJet(*four_momentum) for four_momentum in parts]
    # Construct the numpy array
    particles_numpy = np.array(parts)
    print(f"shape: {particles_numpy.shape}")
    # Select the test case.
    if list_type == "numpy":
        particles = particles_numpy
    elif list_type == "pandas":
        df_particles = pd.DataFrame(parts, columns = ["px", "py", "pz", "E"])
        # Ensure that pz is treated properly.
        df_particles["pz"] = df_particles["pz"].astype(np.float64)
        #print(f"df_particles: {df_particles}")
        #print(f"df dtypes: {df_particles.dtypes}")
        # Must return as a numpy array.
        # NOTE: Do not need to call `np.asarray(...)` or `np.ascontiguousarray(...)` here.
        particles = df_particles.to_numpy()
        #print(f"dtype of particles converted from df: {particles.dtype}")
        #print(f"shape of particles converted from df: {particles.shape}")
        # Sanity check.
        np.testing.assert_allclose(particles, particles_numpy)
    else:
        particles = particles_pseudojets
    print(f"particles: {particles}, particles_numpy: {particles_numpy}, type: {type(particles_numpy)}")

    # choose a jet definition
    R = 0.7
    jet_def = fj.JetDefinition(fj.JetAlgorithm.antikt_algorithm, R)

    # run the clustering, extract the jets
    cs = fj.ClusterSequence(particles, jet_def)
    jets = fj.sorted_by_pt(cs.inclusive_jets())

    # Test direct numpy access
    pt_values = cs.to_numpy()
    print(f"pt_values: type: {type(pt_values)}, values: {pt_values}")

    # Expected jets
    expected_jets = [
        fj.PseudoJet(103.0, 0.0, 0.0, 105.0),
        fj.PseudoJet(-99.0, 0.0, 0.0, 99.0),
    ]
    print(f"measured jets: {jets}")
    print(f"expected jets: {expected_jets}")
    # Check four momentum (they're not really equivalent, so we can't test directly for equivalence).
    assert all([np.isclose(measured.px, expected.px) and np.isclose(measured.py, expected.py) and
                np.isclose(measured.pz, expected.pz) and np.isclose(measured.E, expected.E)
                for measured, expected in zip(jets, expected_jets)])

    # print out some infos
    print(f"Clustering with {jet_def.description()}")

    # print the jets
    print("        pt y phi")
    for i, jet in enumerate(jets):
        print(f"jet {i}: {jet.pt} {jet.rap} {jet.phi}, {jet.px}, {jet.py}, {jet.pz}, {jet.E}")
        for j, constituent in enumerate(jet.constituents):
            print(f"    constituent {j}'s pt: {constituent.pt}")

    print("Trying direct iteration:")
    for jet in cs:
        print(f"jet for cs: {jet}")
        for j, constituent in enumerate(jet.constituents):
            print(f"    constituent {j}: {constituent}")

def test_find_jets() -> None:
    particles = ak.Array({
        "px": [[99.0, 4.0, -99.0], [0.1, -0.1, 0]],
        "py": [[0.1, -0.1, 0], [99.0, 4.0, -99.0]],
        "pz": [[0, 0, 0], [0, 0, 0]],
        "E": [[100.0, 5.0, 99.0], [100.0, 5.0, 99.0]],
    })

    jet_R = 0.7
    jet_defintion = fj.JetDefinition(fj.JetAlgorithm.antikt_algorithm, R = jet_R)
    area_definition = fj.AreaDefinition(fj.AreaType.passive_area, fj.GhostedAreaSpec(1, 1, 0.05))
    settings = fj.JetFinderSettings(jet_definition=jet_defintion, area_definition=area_definition)
    jets, constituents = fj.find_jets(events=particles, settings=settings)

    assert False

if __name__ == "__main__":
    test_fastjet()
