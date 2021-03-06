#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/AreaDefinition.hh>
#include <fastjet/GhostedAreaSpec.hh>

#include <awkward/array/NumpyArray.h>
#include <awkward/builder/ArrayBuilder.h>
#include <awkward/builder/ArrayBuilderOptions.h>
#include <awkward/kernel-dispatch.h>

namespace py = pybind11;
// Shorthand for literals
using namespace pybind11::literals;

// Convenience for awkward.
namespace ak = awkward;

template<typename T>
class IterableWrapper {
 private:
  T fIterable;
  size_t fIndex;
  size_t fSize;
 public:
  IterableWrapper(T& iterable): fIterable(iterable), fIndex(0), fSize(iterable.size()) {};
  typename T::value_type& operator*() { return fIterable.at(fIndex); }
  IterableWrapper& operator++() {
    fIndex++; return *this;
  }
  const size_t index() const { return fIndex; }
  const size_t size() const { return fSize; }
};

/*template<typename T>
class IterableWrapper {
 private:
  T fIterable;
  typename T::iterator fIterator;
  typename T::iterator fEnd;
 public:
  IterableWrapper(T& iterable): fIterable(iterable), fIterator(iterable.begin()), fEnd(iterable.end()) {};
  typename T::value_type& operator*() { return *fIterator; }
  IterableWrapper& operator++() {
    ++fIterator; return *this;
  }
  typename T::iterator begin() const { return fIterable.begin(); }
  typename T::iterator end() const { return fEnd; }
};*/

class IterableWrapperSentinel {};

/**
 * Determines whether we reached the end of the iterator.
 *
 * @return True if we reached the end of the iterator.
 */
template<typename T>
bool operator==(const IterableWrapper<T> & it, const IterableWrapperSentinel &) {
  return it.index() == it.size();
}

/*template<typename T>
bool operator==(const IterableWrapper<T> & it, const IterableWrapperSentinel &) {
  return it == it.end();
}*/

/**
 * Create PseudoJet objects from a numpy array of px, py, pz, E. Axis 0 is the number of particles,
 * while axis 1 must be the 4 parameters.
 *
 * Note: The aray is required to be c-style, which ensures that it works with other packages. For example,
 *       pandas caused a problem in some cases without that argument.
 *
 * @param[jets] Numpy input array.
 * @returns Vector of PseudoJets.
 */
std::vector<fastjet::PseudoJet> constructPseudojetsFromNumpy(const py::array_t<double, py::array::c_style | py::array::forcecast> & jets)
{
  // Retrieve array and relevant information
  py::buffer_info info = jets.request();
  // I'm not sure which one of these is better.
  //auto inputJets = static_cast<double *>(info.ptr);
  auto inputJets = jets.data();
  std::vector<fastjet::PseudoJet> outputJets;
  // This defines our numpy array shape.
  int nParticles = info.shape[0];
  int nParams = info.shape[1];
  //std::cout << "nParams: " << nParams << ", nParticles: " << nParticles << "\n";

  // Validation.
  if (nParams != 4) {
    throw std::runtime_error("Number of params is not correct. Should be four per particle.");
  }
  // Convert the arrays
  for (size_t i = 0; i < nParticles; ++i) {
    /*std::cout << "i: " << i << " inputs: " << inputJets[i * nParams + 0] << " " << inputJets[i * nParams + 1]
      << " " << inputJets[i * nParams + 2] << " " <<  inputJets[i * nParams + 3] << "\n";*/
    outputJets.push_back(fastjet::PseudoJet(
      inputJets[i * nParams + 0], inputJets[i * nParams + 1],
      inputJets[i * nParams + 2], inputJets[i * nParams + 3]));
  }

  return outputJets;
}

struct JetFinderSettings {
  fastjet::JetDefinition _jetDefinition;
  fastjet::AreaDefinition _areaDefinition;

  const fastjet::JetDefinition JetDefinition() { return _jetDefinition; }
  void SetJetDefinition(fastjet::JetDefinition & def) { _jetDefinition = def; }
  const fastjet::AreaDefinition AreaDefinition() { return _areaDefinition; }
  void SetAreaDefinition(fastjet::AreaDefinition & area) { _areaDefinition = area; }
};

template <typename T>
std::vector<fastjet::PseudoJet> getPseudoJetsFromParticles(const std::shared_ptr<ak::Content>& particles)
{
  std::vector<std::pair<int64_t, T*>> arrays;
  // NOTE: Ordering is important here! It needs to match the PseudoJet constructor.
  std::vector<std::string> columnNames = { "px", "py", "pz", "E" };
  for (const auto& name : columnNames) {
    auto npArray = std::dynamic_pointer_cast<ak::NumpyArray>(particles->getitem_field(name));
    arrays.emplace_back(std::make_pair(npArray->length(), reinterpret_cast<T*>(npArray->data())));
  }

  if (arrays.size() != 4) {
    throw std::runtime_error("Could not find the columns in the passed array.");
  }

  // Explicitly expecting 1D arrays.
  std::vector<fastjet::PseudoJet> pseudoJets;
  for (int64_t i = 0; i < arrays[0].first; i++) {
    pseudoJets.emplace_back(
     fastjet::PseudoJet(arrays[0].second[i], arrays[1].second[i], arrays[2].second[i], arrays[3].second[i]));
    // Keep track of the user index so we can match them later.
    pseudoJets[i].set_user_index(i);
  }

  return pseudoJets;
}

std::tuple<std::shared_ptr<ak::Content>, std::shared_ptr<ak::Content>> findJets(const std::shared_ptr<ak::Content> & arr, JetFinderSettings settings)
{
  ak::ArrayBuilder builder(ak::ArrayBuilderOptions(1024, 2.0));
  ak::ArrayBuilder builderConstituents(ak::ArrayBuilderOptions(1024, 2.0));
  for (int64_t i = 0; i < arr->length(); i++) {
    auto particles = arr->getitem_at(i);
    std::vector<fastjet::PseudoJet> pseudoJets;

    // Make sure that we get the right types.
    auto npArray = std::dynamic_pointer_cast<ak::NumpyArray>(particles->getitem_field("px"));
    auto dtype = npArray->dtype();
    if (dtype == ak::util::dtype::float64) {
      pseudoJets = getPseudoJetsFromParticles<double>(particles);
    }
    else if (dtype == ak::util::dtype::float32) {
      pseudoJets = getPseudoJetsFromParticles<float>(particles);
    }
    else {
      throw std::runtime_error("Can only handle float64");
    }

    // Perform jet finding
    fastjet::ClusterSequenceArea cs(pseudoJets, settings.JetDefinition(), settings.AreaDefinition());
    auto jets = fastjet::sorted_by_pt(cs.inclusive_jets(0));

    // Store the resulting jets.
    builder.beginlist();
    for (auto jet : jets) {
      builder.beginrecord_fast("LorentzVector");
      builder.field_fast("px");
      builder.real(jet.px());
      builder.field_fast("py");
      builder.real(jet.py());
      builder.field_fast("pz");
      builder.real(jet.pz());
      builder.field_fast("E");
      builder.real(jet.E());
      builder.endrecord();
    }
    builder.endlist();
    builderConstituents.beginlist();
    for (auto jet: jets) {
        builderConstituents.beginlist();
        for (auto constituent : jet.constituents()) {
            builderConstituents.real(constituent.user_index());
        }
        builderConstituents.endlist();
    }
    builderConstituents.endlist();
  }

  return std::make_tuple(builder.snapshot(), builderConstituents.snapshot());
}


PYBIND11_MODULE(_src, m) {
  using namespace fastjet;
  /// Jet algorithm definitions
  py::enum_<JetAlgorithm>(m, "JetAlgorithm", py::arithmetic(), "Jet algorithms")
    .value("kt_algorithm", JetAlgorithm::kt_algorithm, "the longitudinally invariant kt algorithm")
    .value("cambridge_algorithm", JetAlgorithm::cambridge_algorithm, "the longitudinally invariant variant of the cambridge algorithm (aka Aachen algoithm).")
    .value("aachen_algorithm", JetAlgorithm::cambridge_algorithm, "the longitudinally invariant variant of the cambridge algorithm (aka Aachen algoithm).")
    .value("cambridge_aachen_algorithm", JetAlgorithm::cambridge_algorithm, "the longitudinally invariant variant of the cambridge algorithm (aka Aachen algoithm).")
    .value("antikt_algorithm", JetAlgorithm::antikt_algorithm, "like the k_t but with distance measures dij = min(1/kti^2,1/ktj^2) Delta R_{ij}^2 / R^2 diB = 1/kti^2")
    .value("genkt_algorithm", JetAlgorithm::genkt_algorithm, "like the k_t but with distance measures dij = min(kti^{2p},ktj^{2p}) Delta R_{ij}^2 / R^2 diB = 1/kti^{2p} where p = extra_param()")
    .value("cambridge_for_passive_algorithm", JetAlgorithm::cambridge_for_passive_algorithm, "a version of cambridge with a special distance measure for particles whose pt is < extra_param(); this is not usually intended for end users, but is instead automatically selected when requesting a passive Cambridge area.")
    .value("genkt_for_passive_algorithm", JetAlgorithm::genkt_for_passive_algorithm, "a version of genkt with a special distance measure for particles whose pt is < extra_param() [relevant for passive areas when p<=0] ***** NB: THERE IS CURRENTLY NO IMPLEMENTATION FOR THIS ALG *******")
    .value("ee_kt_algorithm", JetAlgorithm::ee_kt_algorithm, "the e+e- kt algorithm")
    .value("ee_genkt_algorithm", JetAlgorithm::ee_genkt_algorithm, "the e+e- genkt algorithm (R > 2 and p=1 gives ee_kt)")
    .value("plugin_algorithm", JetAlgorithm::plugin_algorithm, "any plugin algorithm supplied by the user")
    .value("undefined_jet_algorithm", JetAlgorithm::undefined_jet_algorithm, "the value for the jet algorithm in a JetDefinition for which no algorithm has yet been defined")
    .export_values();

  // Recombination scheme definitions
  py::enum_<RecombinationScheme>(m, "RecombinationScheme", py::arithmetic(), "Recombination schemes")
    .value("E_scheme", RecombinationScheme::E_scheme, "summing the 4-momenta")
    .value("pt_scheme", RecombinationScheme::pt_scheme, "pt weighted recombination of y,phi (and summing of pt's) with preprocessing to make things massless by rescaling E=| p|")
    .value("pt2_scheme", RecombinationScheme::pt2_scheme, "pt^2 weighted recombination of y,phi (and summing of pt's) with preprocessing to make things massless by rescaling E=| p|")
    .value("Et_scheme", RecombinationScheme::Et_scheme, "pt weighted recombination of y,phi (and summing of pt's) with preprocessing to make things massless by rescaling | p|->=E")
    .value("Et2_scheme", RecombinationScheme::Et2_scheme, "pt^2 weighted recombination of y,phi (and summing of pt's) with preprocessing to make things massless by rescaling | p|->=E")
    .value("BIpt_scheme", RecombinationScheme::BIpt_scheme, "pt weighted recombination of y,phi (and summing of pt's), with no preprocessing")
    .value("BIpt2_scheme", RecombinationScheme::BIpt2_scheme, "pt^2 weighted recombination of y,phi (and summing of pt's) no preprocessing")
    .value("WTA_pt_scheme", RecombinationScheme::WTA_pt_scheme, "pt-based Winner-Takes-All (WTA) recombination: the result of the recombination has the rapidity, azimuth and mass of the PseudoJet with the larger pt, and a pt equal to the sum of the two pt's")
    .value("WTA_modp_scheme", RecombinationScheme::WTA_modp_scheme, "mod-p-based Winner-Takes-All (WTA) recombination: the result of the recombination gets the 3-vector direction and mass of the PseudoJet with the larger |3-momentum| (modp), and a |3-momentum| equal to the scalar sum of the two |3-momenta|.")
    .value("external_scheme", RecombinationScheme::external_scheme, "for the user's external scheme")
    .export_values();

  // Jet reconstruction strategy definitions
  py::enum_<Strategy>(m, "Strategy", "Jet reconstruction strategies")
    .value("N2MHTLazy9AntiKtSeparateGhosts", Strategy::N2MHTLazy9AntiKtSeparateGhosts, "Like N2MHTLazy9 in a number of respects, but does not calculate ghost-ghost distances and so does not carry out ghost-ghost recombination.")
    .value("N2MHTLazy9", Strategy::N2MHTLazy9, "Only looks into a neighbouring tile for a particle's nearest neighbour (NN) if that particle's in-tile NN is further than the distance to the edge of the neighbouring tile. Uses tiles of size R and a 3x3 tile grid around the particle.")
    .value("N2MHTLazy25", Strategy::N2MHTLazy25, " Similar to N2MHTLazy9, but uses tiles of size R/2 and a 5x5 tile grid around the particle.")
    .value("N2MinHeapTiled", Strategy::N2MinHeapTiled, "faster that N2Tiled above about 500 particles; differs from it by retainig the di(closest j) distances in a MinHeap (sort of priority queue) rather than a simple vector. ")
    .value("N2Tiled", Strategy::N2Tiled, "fastest from about 50..500")
    .value("N2Plain", Strategy::N2Plain, "fastest below 50")
    .value("Best", Strategy::Best, "automatic selection of the best (based on N), including the LazyTiled strategies that are new to FJ3.1")
    .value("NlnN", Strategy::NlnN, "best of the NlnN variants -- best overall for N>10^4. (Does not work for R>=2pi)")
    .value("NlnN3pi", Strategy::NlnN3pi, "legacy N ln N using 3pi coverage of cylinder. (Does not work for R>=2pi)")
    .value("NlnN4pi", Strategy::NlnN4pi, "legacy N ln N using 4pi coverage of cylinder")
    .value("NlnNCam4pi", Strategy::NlnNCam4pi, "Chan's closest pair method (in a variant with 4pi coverage), for use exclusively with the Cambridge algorithm. (Does not work for R>=2pi)")
    .value("NlnNCam2pi2R", Strategy::NlnNCam2pi2R, "Chan's closest pair method (in a variant with 2pi+2R coverage), for use exclusively with the Cambridge algorithm.  (Does not work for R>=2pi)")
    .value("NlnNCam", Strategy::NlnNCam, "Chan's closest pair method (in a variant with 2pi+minimal extra variant), for use exclusively with the Cambridge algorithm. (Does not work for R>=2pi)")
    .value("BestFJ30", Strategy::BestFJ30, "the automatic strategy choice that was being made in FJ 3.0 (restricted to strategies that were present in FJ 3.0)")
    .value("plugin_strategy", Strategy::plugin_strategy, "the plugin has been used...")
    .export_values();

  py::class_<JetDefinition>(m, "JetDefinition", "Jet definition")
    .def(py::init<JetAlgorithm, RecombinationScheme, Strategy>(), "jet_algorithm"_a, "recombination_scheme"_a = E_scheme, "strategy"_a = Best)
    .def(py::init<JetAlgorithm, double, RecombinationScheme, Strategy>(), "jet_algorithm"_a, "R"_a, "recombination_scheme"_a = E_scheme, "strategy"_a = Best)
    .def(py::init<JetAlgorithm, double, double, RecombinationScheme, Strategy>(), "jet_algorithm"_a, "R"_a, "extra_parameter"_a, "recombination_scheme"_a = E_scheme, "strategy"_a = Best)
    // Leave out implementation of an external recombiner because it's not entirely clear how we should handle that in the bindings.
    .def(py::init<JetAlgorithm, double, RecombinationScheme, Strategy, int>(), "jet_algorithm"_a, "R"_a, "recombination_scheme"_a, "strategy"_a, "n_parameters"_a)
    .def_property("recombination_scheme", &JetDefinition::recombination_scheme, &JetDefinition::set_recombination_scheme, "The recombination scheme.")
    .def_property("jet_algorithm", &JetDefinition::jet_algorithm, &JetDefinition::set_jet_algorithm, "The jet algorithm.")
    .def_property_readonly("R", &JetDefinition::R, "The jet resolution parameter.")
    .def_property("extra_parameter", &JetDefinition::extra_param, &JetDefinition::set_extra_param, "A general purpose extra parameter, whose meaning depends on the algorithm, and may often be unused.")
    .def_property_readonly("strategy", &JetDefinition::strategy, "Jet finding strategy")
    .def("description", [](JetDefinition & jet, const bool include_recombiner) { if (include_recombiner == true) { return jet.description(); } return jet.description_no_recombiner(); }, "include_recombiner"_a = true, "A textual description of the current jet definition. The recombiner description is included by default.")
    .def_static("algorithm_description", &JetDefinition::algorithm_description, "jet_algorithm"_a, "A short textual description of the algorithm jet_algorithm")
    .def("__repr__", [](JetDefinition& jetDef){
        std::stringstream s;
        //s << "<JetDefinition jet_algorithm=" << jetDef.jet_algorithm() << ", R=" << jetDef.R() << " at " << &jetDef << ">";
        s << "<JetDefinition R=" << jetDef.R() << " at " << &jetDef << ">";
        return s.str();
      })
    ;

  py::class_<PseudoJet>(m, "PseudoJet")
    .def(py::init<double, double, double, double>(), "px"_a, "py"_a, "pz"_a, "E"_a)
    .def_property_readonly("e", &PseudoJet::E)
    // Alias for `e`
    .def_property_readonly("E", &PseudoJet::E)
    .def_property_readonly("px", &PseudoJet::px)
    .def_property_readonly("py", &PseudoJet::py)
    .def_property_readonly("pz", &PseudoJet::pz)
    .def_property_readonly("phi", &PseudoJet::phi, "phi (in the range 0..2pi)")
    .def_property_readonly("phi_std", &PseudoJet::phi_std, "phi in the range -pi..pi")
    .def_property_readonly("phi_02pi", &PseudoJet::phi_02pi, "phi in the range 0..2pi")
    .def_property_readonly("rap", &PseudoJet::rap, "the rapidity or some large value when the rapidity is infinite.")
    // Alias for `rap`
    .def_property_readonly("rapidity", &PseudoJet::rapidity, "the rapidity or some large value when the rapidity is infinite.")
    .def_property_readonly("pseudorapidity", &PseudoJet::pseudorapidity, "the pseudo-rapidity or some large value when the rapidity is infinite.")
    // Alias for `pseudorapidity`
    .def_property_readonly("eta", &PseudoJet::eta, "the pseudo-rapidity or some large value when the rapidity is infinite.")
    .def_property_readonly("pt", &PseudoJet::pt, "the scalar transverse momentum")
    .def_property_readonly("pt2", &PseudoJet::pt2, "the squared transverse momentum")
    .def_property_readonly("perp", &PseudoJet::perp, "the scalar transverse momentum")
    .def_property_readonly("perp2", &PseudoJet::perp2, "the squared transverse momentum")
    .def_property_readonly("m", &PseudoJet::m, "the squared invariant mass")
    .def_property_readonly("m2", &PseudoJet::m2, "the invariant mass")
    .def_property_readonly("mt", &PseudoJet::mt, "the squared transverse mass = kt^2+m^2")
    .def_property_readonly("mt2", &PseudoJet::mt2, "the transverse mass = sqrt(kt^2+m^2)")
    .def_property_readonly("mperp", &PseudoJet::mperp, "the squared transverse mass = kt^2+m^2")
    .def_property_readonly("mperp2", &PseudoJet::mperp2, "the transverse mass = sqrt(kt^2+m^2)")
    .def_property_readonly("modp", &PseudoJet::modp, "the 3-vector modulus = sqrt(px^2+py^2+pz^2)")
    .def_property_readonly("modp2", &PseudoJet::modp, "the squared 3-vector modulus = sqrt(px^2+py^2+pz^2)")
    .def_property_readonly("et", &PseudoJet::Et, "the transverse energy")
    .def_property_readonly("et2", &PseudoJet::Et2, "the squared transverse energy")
    // Intro ducted in a more recent versions of fastjet.
    /*.def_property_readonly("cos_theta", &PseudoJet::cos_theta, "cosine of the polar angle")
    .def_property_readonly("theta", &PseudoJet::theta, "the polar angle")*/
    .def_property_readonly("beam_distance", &PseudoJet::beam_distance, "distance between this jet and the beam")
    // Methods
    .def("delta_R", &PseudoJet::delta_R, "other"_a)
    .def("boost", &PseudoJet::boost, "rest_frame"_a, "transform this jet (given in the rest frame of the given PseudoJet) into a jet in the lab frame")
    .def("unboost", &PseudoJet::boost, "rest_frame"_a, " transform this jet (given in lab) into a jet in the rest frame of the given PseudoJet")
    // Operations
    .def(py::self + py::self, "other"_a, "Add the given PseudoJets.")
    .def(py::self - py::self, "other"_a, "Subtract the given PseudoJets.")
    .def(py::self += py::self, "other"_a, "Add the given PseudoJet.")
    .def(py::self -= py::self, "other"_a, "Subtract the given PseudoJet.")
    .def(py::self * double(), "val"_a, "Scale by the given value.")
    .def(py::self *= double(), "val"_a, "Scale by the given value.")
    .def(py::self / double(), "val"_a, "Divide by the given value.")
    .def(py::self /= double(), "val"_a, "Divide by the given value.")
    .def(py::self == py::self, "other"_a, "Compare the given PseudoJets. Returns true if the 4 momentum components of the two PseudoJets are identical and all the internal indices (user, cluster_history)")
    .def(py::self != py::self, "other"_a, "Compare the given PseudoJets. Returns true if they are not equal.")
    // Reset
    // Cannot use the templated version directly - we have to instantiate types explicitly.
    .def("reset", (void (PseudoJet::*)(double, double, double, double)) &PseudoJet::reset, "px"_a, "py"_a, "pz"_a, "E"_a, "Reset the PseudoJet to the properties of the given jet px, py, pz, E.")
    .def("reset", &PseudoJet::reset<std::vector<double>>, "jet_four_momentum"_a, "Reset the PseudoJet to the properties of the given jet px, py, py, E.")
    .def("reset", &PseudoJet::reset<PseudoJet>, "new_psueodjet"_a, "Reset the PseudoJet to the properties of the given jet.")
    .def("reset_momentum", (void (PseudoJet::*)(double, double, double, double)) &PseudoJet::reset_momentum, "px"_a, "py"_a, "pz"_a, "E"_a, "Reset the momentum PseudoJet to the momentum of the given jet.")
    .def("reset_momentum", (void (PseudoJet::*)(const PseudoJet &)) &PseudoJet::reset_momentum, "new_psueodjet"_a, "Reset the momentum PseudoJet to the momentum of the given jet.")
    // User index
    .def_property("user_index", &PseudoJet::user_index, &PseudoJet::set_user_index, "The user index allows the user to add simple identifying information to a particle/jet")
    // Area
    .def_property_readonly("area", &PseudoJet::area, "The jet (scalar) area.")
    .def_property_readonly("area_error", &PseudoJet::area_error, "The error (uncertainty) associated with the determination of the area of this jet.")
    .def_property_readonly("area_4vector", &PseudoJet::area_4vector, "The jet 4-vector area.")
    .def_property_readonly("is_pure_ghost", &PseudoJet::is_pure_ghost, "True if this jet is made exclusively of ghosts.")
    // Constituents
    .def_property_readonly("constituents", &PseudoJet::constituents, "The constituents.")
    .def("__getitem__", [](const PseudoJet &jet, size_t i) {
      auto constituents = jet.constituents();
      if (i >= constituents.size()) throw py::index_error();
      return constituents[i];
    }, "i"_a, R"pbdoc(
      Retrieve an individual constituent at a given index.

      Args:
        i: Index of the constituent to retrieve.
      Returns:
        Constituent at that index.
    )pbdoc")
    .def("__iter__",
      [](const PseudoJet &jet) {
        auto constituents = jet.constituents();
        return py::make_iterator(IterableWrapper<std::vector<PseudoJet>>(constituents), IterableWrapperSentinel());
      }, py::keep_alive<0, 1>(), R"pbdoc(
        Iterates over the constituents in the PseudoJet.
      )pbdoc")
    // Cluster sequence related
    .def_property_readonly("has_associated_cluster_sequence", &PseudoJet::has_associated_cluster_sequence, "true if this PseudoJet has an associated ClusterSequence")
    .def_property_readonly("has_associated_cs", &PseudoJet::has_associated_cluster_sequence, "true if this PseudoJet has an associated ClusterSequence")
    .def_property_readonly("has_valid_cluster_sequence", &PseudoJet::has_valid_cluster_sequence, "true if this PseudoJet has an associated and still valid(ated) ClusterSequence.")
    .def_property_readonly("has_valid_cs", &PseudoJet::has_valid_cluster_sequence, "true if this PseudoJet has an associated and still valid(ated) ClusterSequence.")
    .def("has_partner", &PseudoJet::has_partner, "partner"_a, "Check if it has been recombined with another PseudoJet in which case, return its partner through the argument. Otherwise, 'partner' is set to 0.")
    .def("has_child", &PseudoJet::has_child, "child"_a, "Check if it has been recombined with another PseudoJet in which case, return its child through the argument. Otherwise, 'child' is set to 0.")
    .def("parents", [](PseudoJet & jet) {
        PseudoJet parent1, parent2;
        jet.has_parents(parent1, parent2);
        // We'd like to return None if the parents aren't available. However, we can't handle this nicely without c++17.
        // Instead, we just return both the parents. The user has to check...
        //return (parent1.pt() != 0 || parent2.pt() != 0) ? std::optional<std::tuple<PseudoJet, PseudoJet>>{std::make_tuple(parent1, parent2)} : py::cast<py::none>(Py_None);
        return std::make_tuple(parent1, parent2);
      }, "Return the parents of the PseudoJet if it is the product of a recombination. If not, returns two empty PseudoJets.")
    .def("contains", &PseudoJet::contains, "other"_a, " Check if the current PseudoJet contains the one passed as argument.")
    .def("is_inside", &PseudoJet::is_inside, "other"_a, "Check if the current PseudoJet is contained the one passed as argument.")
    .def("__repr__", [](PseudoJet & jet){
        std::stringstream s;
        //s << "<PseudoJet px=" << jet.px() << ", py=" << jet.py() << ", pz=" << jet.pz() << ", E=" << jet.E();
        s << "<PseudoJet pt=" << jet.pt() << ", eta=" << jet.eta() << ", phi=" << jet.phi() << ", m=" << jet.m() << " at " << &jet << ">";
        return s.str();
      });

  // Helper functions
  m.def("dot_product", &dot_product, "jet_1"_a, "jet_2"_a, "Returns the 4-vector dot product of a and b");
  m.def("have_same_momentum", &have_same_momentum, "jet_1"_a, "jet_2"_a, "Returns true if the momenta of the two input jets are identical");
  m.def("sorted_by_pt", &sorted_by_pt, "jets"_a, "Return a vector of jets sorted into decreasing transverse momentum");
  m.def("sorted_by_pz", &sorted_by_pz, "jets"_a, "Return a vector of jets sorted into increasing pz");
  m.def("sorted_by_rapidity", &sorted_by_rapidity, "jets"_a, "Return a vector of jets sorted into increasing rapidity");
  m.def("sorted_by_E", &sorted_by_E, "jets"_a, "Return a vector of jets sorted into decreasing energy");

  py::class_<ClusterSequence>(m, "ClusterSequence")
    .def(py::init<const std::vector<PseudoJet> &, const JetDefinition &, const bool &>(), "pseudojets"_a, "jet_definition"_a, "write_out_combinations"_a = false, "Create a ClusterSequence, starting from the supplied set of PseudoJets and clustering them with jet definition specified by jet_definition (which also specifies the clustering strategy)")
    // numpy constructor.
    .def(py::init([](const py::array_t<double> & pseudojets, const JetDefinition & jetDef, const bool & writeOutCombination){ auto convertedPseudojets = constructPseudojetsFromNumpy(pseudojets); return ClusterSequence(convertedPseudojets, jetDef, writeOutCombination); }), "pseudojets"_a, "jet_definition"_a, "write_out_combinations"_a = false, "Create a ClusterSequence, starting from the supplied set of PseudoJets and clustering them with jet definition specified by jet_definition (which also specifies the clustering strategy)")
    .def("inclusive_jets", &ClusterSequence::inclusive_jets, "pt_min"_a = 0., "Return a vector of all jets (in the sense of the inclusive algorithm) with pt >= ptmin. Time taken should be of the order of the number of jets returned.")
    .def("__getitem__", [](const ClusterSequence &cs, size_t i) {
      auto inclusive_jets = cs.inclusive_jets();
      if (i >= inclusive_jets.size()) throw py::index_error();
      return inclusive_jets[i];
    }, "i"_a, R"pbdoc(
      Retrieve an individual jet at a given index.

      Args:
        i: Index of the jet to retrieve.
      Returns:
        Jet at that index.
    )pbdoc")
    .def("__iter__",
      [](const ClusterSequence &cs) {
        auto jets = cs.inclusive_jets();
        return py::make_iterator(IterableWrapper<std::vector<PseudoJet>>(jets), IterableWrapperSentinel());
      }, py::keep_alive<0, 1>(), R"pbdoc(
        Finds the include jets and iterates over them.
      )pbdoc")
    .def("__call__",
      [](const ClusterSequence &cs, double min_pt = 0) {
        auto jets = cs.inclusive_jets(min_pt);
        return py::make_iterator(IterableWrapper<std::vector<PseudoJet>>(jets), IterableWrapperSentinel());
      }, py::keep_alive<0, 1>(), "min_pt"_a = 0, R"pbdoc(
        Retrieves the inclusive jets.

        Args:
          min_pt: Minimum jet pt to include. Default: 0.
        Returns:
          List of inclusive jets.
      )pbdoc")
    .def("to_numpy",
      [](const ClusterSequence &cs, double min_pt = 0) {
        auto jets = cs.inclusive_jets(min_pt);
        // Don't specify the size if using push_back.
        std::vector<double> pt, eta, phi, m;
        for (const auto & jet : jets) {
          pt.push_back(jet.pt());
          eta.push_back(jet.eta());
          phi.push_back(jet.phi());
          m.push_back(jet.m());
        }
        return std::make_tuple(
            py::array(py::cast(pt)),
            py::array(py::cast(eta)),
            py::array(py::cast(phi)),
            py::array(py::cast(m))
          );
      }, "min_pt"_a = 0, R"pbdoc(
        Retrieves the inclusive jets and converts them to numpy arrays.

        Args:
          min_pt: Minimum jet pt to include. Default: 0.
        Returns:
          pt, eta, phi, m of inclusive jets.
      )pbdoc");

  py::class_<ClusterSequenceArea, ClusterSequence>(m, "ClusterSequenceArea")
    .def(py::init<const std::vector<PseudoJet> &, const JetDefinition &, const AreaDefinition &>(), "pseudojets"_a, "jet_definition"_a, "area_definition"_a, "Create a ClusterSequenceArea, starting from the supplied set of PseudoJets and clustering them with jet definition specified by jet_definition (which also specifies the clustering strategy)")
    .def(py::init<const std::vector<PseudoJet> &, const JetDefinition &, const GhostedAreaSpec &>(), "pseudojets"_a, "jet_definition"_a, "ghost_spec"_a, "Create a ClusterSequenceArea, starting from the supplied set of PseudoJets and clustering them with jet definition specified by jet_definition (which also specifies the clustering strategy)");

  py::enum_<AreaType>(m, "AreaType", py::arithmetic(), "Jet area definitions")
    .value("invalid_area", AreaType::invalid_area, "Invalid jet area")
    .value("active_area", AreaType::active_area, "Active jet area")
    .value("active_area_explicit_ghosts", AreaType::active_area_explicit_ghosts, "Active jet area with explicit ghosts")
    .value("one_ghost_passive_area", AreaType::one_ghost_passive_area, "One ghost with passive area")
    .value("passive_area", AreaType::passive_area, "Passive area")
    .value("voronoi_area", AreaType::voronoi_area, "Voronoi area")
    .export_values();

  py::class_<GhostedAreaSpec>(m, "GhostedAreaSpec", "Ghost area specification")
    .def(py::init<double, int, double, double, double, double>(),
        "max_rapidity"_a,
        "repeat_in"_a = gas::def_repeat,
        "ghost_area_in"_a = gas::def_ghost_area,
        "grid_scatter_in"_a = gas::def_grid_scatter,
        "pt_scatter_in"_a = gas::def_pt_scatter,
        "mean_ghost_pt_in"_a = gas::def_mean_ghost_pt)
    .def_property("max_rapidity", &GhostedAreaSpec::ghost_maxrap, &GhostedAreaSpec::set_ghost_maxrap, "Max ghost rapidity.")
    .def_property("ghost_area", &GhostedAreaSpec::ghost_area, &GhostedAreaSpec::set_ghost_area, "Ghost area.")
    .def_property("grid_scatter", &GhostedAreaSpec::grid_scatter, &GhostedAreaSpec::set_grid_scatter, "Grid scatter.")
    .def_property("pt_scatter", &GhostedAreaSpec::pt_scatter, &GhostedAreaSpec::set_pt_scatter, "pt scatter.")
    .def_property("mean_ghost_pt", &GhostedAreaSpec::mean_ghost_pt, &GhostedAreaSpec::set_mean_ghost_pt, "Mean ghost pt.");

  py::class_<AreaDefinition>(m, "AreaDefinition", "Area definition")
    .def(py::init<AreaType, const GhostedAreaSpec &>(), "area_type"_a, "ghost_spec"_a)
    .def(py::init<AreaType>(), "area_type"_a)
    .def_property_readonly("area_type", &AreaDefinition::area_type)
    .def_property_readonly("ghost_spec", (const GhostedAreaSpec & (AreaDefinition::*)() const) &AreaDefinition::ghost_spec)
    .def("__repr__", [](AreaDefinition& areaDef){
        std::stringstream s;
        //s << "<AreaDefinition area_type=" << areaDef.area_type() << ", ghost_spec=" << areaDef.ghost_spec() << " at " << &areaDef << ">";
        s << "<AreaDefinition at " << &areaDef << ">";
        return s.str();
      })
    ;

  // Awkward array
  // Ensure dependencies are loaded.
  py::module::import("awkward");

  py::class_<JetFinderSettings>(m, "JetFinderSettings", "Encompasses jet finder settings")
    .def(py::init<const JetDefinition, const AreaDefinition>(), "jet_definition"_a, "area_definition"_a)
    .def_property("jet_definition", &JetFinderSettings::JetDefinition, &JetFinderSettings::SetJetDefinition)
    .def_property("area_definition", &JetFinderSettings::AreaDefinition, &JetFinderSettings::SetAreaDefinition)
    .def("__repr__", [](JetFinderSettings & settings){
        std::stringstream s;
        //s << "<JetFinderSettings jet_definition R=" << settings.JetDefinition() << ", area_definition=" << settings.AreaDefinition()  << " at " << &settings << ">";
        s << "<JetFinderSettings jet_def R=" << settings.JetDefinition().R() << " at " << &settings << ">";
        return s.str();
      })
    ;

  m.def("_find_jets", &findJets, "events"_a, "settings"_a, "Find jets for the given events according to the provided settings.");

  // TODO: fastjet-contrib bindings. Look at a substructure analysis.
  // Constituent subtractor
  // Recursive tools
  // Lund plane
}
