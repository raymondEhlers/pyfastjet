""" Interface for pyfastjet.

.. codeauthor:: Raymond Ehlers <raymond.ehlers@cern.ch>, ORNL
"""

import awkward1 as ak

# TODO: Shift to all explicit import
# TODO: Add in more of the bindings...
#from pyfastjet._src import _find_jets, AreaType, AreaDefinition, ClusterSequence, ClusterSequenceArea, GhostedAreaSpec, JetFinderSettings, JetDefinition, JetAlgorithm, PseudoJet, sorted_by_pt
from pyfastjet._src import *
# Just import what's explicitly used here...
from pyfastjet._src import _find_jets, JetFinderSettings

def find_jets(events: ak.Array, settings: JetFinderSettings) -> ak.Array:
    """ Find jets according to the events and settings.

    This wrapper is provided so we don't have to remember to pass a layout
    and reconstruct in an `ak.Array`. It's all for convenience.

    Args:
        events: Awkward array containing particles within events. Must contain
            "E", "px", "py", and "pz" columns.
        settings: Jet finder settings for configuring a ClusterSequenceArea (ie. the jet
            JetDefinition and GhostedAreaSpec definitions.
    Returns:
        Jets found according to the settings.
    """
    return ak.Array(_find_jets(events=events.layout, settings=settings))

def my_func() -> str:
    print("I'm in __init__!")
