#!/usr/bin/env python3

""" Basic tests for the fastjet binding.

.. code-author:: Raymond Ehlers <raymond.ehlers@cern.ch>, ORNL
"""

import pyfastjet as fj
import numpy as np
from typing import List

import pytest

@pytest.mark.parametrize("use_numpy", [
    False, True
], ids = ["PseudoJets", "Numpy arrays"])
def test_fastjet(use_numpy: bool) -> None:
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
    if use_numpy:
        particles = particles_numpy
    else:
        particles = particles_pseudojets
    print(f"particles: {particles}, particles_numpy: {particles_numpy}, type: {type(particles_numpy)}")

    # choose a jet definition
    R = 0.7
    jet_def = fj.JetDefinition(fj.JetAlgorithm.antikt_algorithm, R)

    # run the clustering, extract the jets
    cs = fj.ClusterSequence(particles, jet_def)
    jets = fj.sorted_by_pt(cs.inclusive_jets())

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

if __name__ == "__main__":
    test_fastjet()
