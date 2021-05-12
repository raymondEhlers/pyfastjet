
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>
#include <iostream>

#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/AreaDefinition.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/tools/JetMedianBackgroundEstimator.hh>
#include <fastjet/tools/Subtractor.hh>
#include <fastjet/contrib/ConstituentSubtractor.hh>

/**
 * Common jet finding approach: Find jets using the anti-kt algorithm.
 */
template <typename T>
std::tuple<std::shared_ptr<fastjet::ClusterSequenceArea>, std::vector<fastjet::PseudoJet>> singleAntiKtJetFinding(
  double jetR,
  // Particle four vectors
  std::vector<T> px,
  std::vector<T> py,
  std::vector<T> pz,
  std::vector<T> E)
{
  // NOTE: Not all users will need to use an area definition.
  //       Note that using a jet area def requires using a `ClusterSequenceArea`.
  fastjet::AreaDefinition areaDefinition = fastjet::AreaDefinition(fastjet::AreaType::passive_area, fastjet::GhostedAreaSpec(1, 1, 0.05));
  // I'm assuming anti-kt here from simplicity.
  fastjet::JetDefinition jetDefinition(fastjet::JetAlgorithm::antikt_algorithm, jetR);

  // Convert inputs to PseudoJet objects
  std::vector<fastjet::PseudoJet> inputPseudoJets;
  for (std::size_t i = 0; i < px.size(); i++) {
    inputPseudoJets.emplace_back(
      fastjet::PseudoJet(
        px.at(i),
        py.at(i),
        pz.at(i),
        E.at(i)
      )
    );
  }

  std::shared_ptr<fastjet::ClusterSequenceArea> cs = std::make_shared<fastjet::ClusterSequenceArea>(inputPseudoJets, jetDefinition, areaDefinition);
  auto jets = fastjet::sorted_by_E(cs->inclusive_jets(0));

  // Now I use my jets for whatever I want.
  // I don't necessarily care about the cs itself, but it needs to stay alive so I can access jet constituents, etc. (detail of fastjet)
  return std::make_tuple(cs, jets);
}


/**
 * Common jet finding approach, but this time taking advantage of background subtraction.
 *
 * Supports either rho subtraction, or event-wise constituent subtraction.
 */
template <typename T>
std::tuple<std::shared_ptr<fastjet::ClusterSequenceArea>, std::vector<fastjet::PseudoJet>> jetFindingWithBackgroundSubtraction(
  double jetR,
  // Particle four vectors
  std::vector<T> px,
  std::vector<T> py,
  std::vector<T> pz,
  std::vector<T> E,
  bool useConstituentSubtraction)
{
  // General settings
  double etaMin = -0.9;
  double etaMax = 0.9;
  // So-called "Ghost" settings
  // They're often the same for both signal and background, so we define it here first.
  // We don't always use all of these options, but tried to be more complete.
  double ghostEtaMin = etaMin;
  double ghostEtaMax = etaMax;
  double ghostArea = 0.005;
  int ghostRepeatN = 1;
  double ghostktMean = 1e-100;
  double gridScatter = 1.0;
  double ktScatter = 0.1;
  fastjet::GhostedAreaSpec ghostAreaSpec(ghostEtaMax, ghostRepeatN, ghostArea, gridScatter, ktScatter, ghostktMean);

  // Background settings
  // These should all be settable, but I wanted to keep the signature simple, so I just define them here with some reasonable(ish) defaults.
  double backgroundJetR = 0.2;
  double backgroundJetEtaMin = etaMin;
  double backgroundJetEtaMax = etaMax;
  double backgroundJetPhiMin = 0;
  double backgroundJetPhiMax = 2 * M_PI;
  // Fastjet background settings
  fastjet::JetAlgorithm backgroundJetAlgorithm(fastjet::kt_algorithm);
  fastjet::RecombinationScheme backgroundRecombinationScheme(fastjet::E_scheme);
  fastjet::Strategy backgroundStrategy(fastjet::Best);
  fastjet::AreaType backgroundAreaType(fastjet::active_area);
  // Derived fastjet settings
  fastjet::JetDefinition backgroundJetDefinition(backgroundJetAlgorithm, backgroundJetR, backgroundRecombinationScheme, backgroundStrategy);
  fastjet::AreaDefinition backgroundAreaDefinition(backgroundAreaType, ghostAreaSpec);
  fastjet::Selector selRho = fastjet::SelectorRapRange(backgroundJetEtaMin, backgroundJetEtaMax) && fastjet::SelectorPhiRange(backgroundJetPhiMin, backgroundJetPhiMax) && !fastjet::SelectorNHardest(2);
  // Constituent subtraction options (if used)
  double constituentSubAlpha = 1.0;
  double constituentSubRMax = 0.25;

  // Finally, define the background estimator
  // This is needed for all background subtraction cases.
  fastjet::JetMedianBackgroundEstimator backgroundEstimator(selRho, backgroundJetDefinition, backgroundAreaDefinition);

  // Signal jet settings
  // Again, these should all be settable, but I wanted to keep the signature simple, so I just define them here with some reasonable(ish) defaults.
  double jetPtMin = 0;
  double jetPtMax = 1000;
  // Would often set as abs(eta - R), but should be configurable.
  double jetEtaMin = etaMin + jetR;
  double jetEtaMax = etaMax - jetR;
  double jetPhiMin = 0;
  double jetPhiMax = 2 * M_PI;
  // Fastjet settings
  fastjet::JetAlgorithm jetAlgorithm(fastjet::antikt_algorithm);
  fastjet::RecombinationScheme recombinationScheme(fastjet::E_scheme);
  fastjet::Strategy strategy(fastjet::Best);
  fastjet::AreaType areaType(fastjet::active_area);
  // Derived fastjet settings
  fastjet::JetDefinition jetDefinition(jetAlgorithm, jetR, recombinationScheme, strategy);
  fastjet::AreaDefinition areaDefinition(areaType, ghostAreaSpec);
  fastjet::Selector selJets = fastjet::SelectorPtRange(jetPtMin, jetPtMax) && fastjet::SelectorEtaRange(jetEtaMin, jetEtaMax) && fastjet::SelectorPhiRange(jetPhiMin, jetPhiMax);

  // Convert inputs to PseudoJet objects
  std::vector<fastjet::PseudoJet> inputPseudoJets;
  for (std::size_t i = 0; i < px.size(); i++) {
    inputPseudoJets.emplace_back(
      fastjet::PseudoJet(
        px.at(i),
        py.at(i),
        pz.at(i),
        E.at(i)
      )
    );
  }

  // Setup the background estimator to be able to make the estimation.
  backgroundEstimator.set_particles(inputPseudoJets);

  // Now, deal with applying the background subtraction.
  // The subtractor will subtract the background from jets. It's not used in the case of constituent subtraction.
  std::shared_ptr<fastjet::Subtractor> subtractor = nullptr;
  // The constituent subtraction (here, it's implemented as event-wise subtraction, but that doesn't matter) takes
  // a different approach to background subtraction. It's used here to illustrate a different work flow.
  std::shared_ptr<fastjet::contrib::ConstituentSubtractor> constituentSubtraction = nullptr;
  // Now, set them up as necessary.
  if (!useConstituentSubtraction) {
    subtractor = std::make_shared<fastjet::Subtractor>(&backgroundEstimator);
  }
  else {
    constituentSubtraction = std::make_shared<fastjet::contrib::ConstituentSubtractor>(&backgroundEstimator);
    constituentSubtraction->set_distance_type(fastjet::contrib::ConstituentSubtractor::deltaR);
    constituentSubtraction->set_max_distance(constituentSubRMax);
    constituentSubtraction->set_alpha(constituentSubAlpha);
    constituentSubtraction->set_ghost_area(ghostArea);
    constituentSubtraction->set_max_eta(backgroundJetEtaMax);
    constituentSubtraction->set_background_estimator(&backgroundEstimator);
  }

  // For constituent subtraction, we subtract the input particles
  if (useConstituentSubtraction) {
    inputPseudoJets = constituentSubtraction->subtract_event(inputPseudoJets);
  }

  // Perform jet finding on signal
  std::shared_ptr<fastjet::ClusterSequenceArea> cs = std::make_shared<fastjet::ClusterSequenceArea>(inputPseudoJets, jetDefinition, areaDefinition);
  // Get the jets, implicitly grabbing jets down to pt of 0.
  auto jets = cs->inclusive_jets(0);
  // Apply the subtractor when appropriate
  if (!useConstituentSubtraction) {
    jets = (*subtractor)(jets);
  }
  // It's also not uncommon to apply a sorting by E or pt.
  jets = fastjet::sorted_by_E(jets);

  // Now I use my jets for whatever I want.
  // I don't necessarily care about the cs itself, but it needs to stay alive so I can access jet constituents, etc. (detail of fastjet)
  return std::make_tuple(cs, jets);
}

int main() {
  // Inputs to functions are vectors of px, py, pz, and E.
  // Would be better to get some larger tests inputs.
  // Could use power law of ~-4 for pt, flat in eta from -1 to 1, and full 0-2pi in phi
  std::cout << "Simple\n";
  singleAntiKtJetFinding(0.5, std::vector<double>{1.0},
    std::vector<double>{1.0},std::vector<double>{1.0},std::vector<double>{1.0});
  std::cout << "Rho subtraction\n";
  jetFindingWithBackgroundSubtraction(0.5, std::vector<double>{1.0},
    std::vector<double>{1.0},std::vector<double>{1.0},std::vector<double>{1.0}, false);
  std::cout << "Constituent subtraction (event-wise)\n";
  jetFindingWithBackgroundSubtraction(0.5, std::vector<double>{1.0},
    std::vector<double>{1.0},std::vector<double>{1.0},std::vector<double>{1.0}, true);
  std::cout << "Done!\n";
}
