#ifndef MMT_ANALYSIS_H
#define MMT_ANALYSIS_H

#include "../entities/Universe.h"

namespace pylimer_tools::calc::mmt {

/*
Compute the stoichiometric inbalance
( nr. of bonds formable of crossLinker / nr. of formable bonds of precursor
)

NOTE:
  - if your system has a non-integer number of possible bonds (e.g. one site
unbonded), this will not be rounded/respected in any way.

Arguments:
  - network: The poylmer network to do the computation for
  - crossLinkerType: The type of the junctions/crossLinkers to select them
in the network
  - strandLength: The length of the network strands (in nr. of beads).
          Used to infer the number of precursor strands.
          If `None`: will use average length of each connected system when
ignoring the crossLinkers.
  - functionalityPerType: a dictionary with key: type, and value:
functionality of this atom type. If `None`: will use max functionality per
type.

Returns:
  - r (float): The stoichiometric inbalance
*/
double
computeStoichiometricInbalance(pylimer_tools::entities::Universe network,
                               const int crossLinkerType,
                               const int strandLength,
                               std::map<int, int> functionalityPerType)
{

  if (network.getNrOfAtoms() == 0) {
    return 0;
  }

  std::vector<int> allTypes = network.getPropertyValues<int>("type");

  int crossLinkerFormableBonds = 0;
  int otherFormableBonds = 0;

  // count how often every type occurs
  for (int type : allTypes) {
    if (type == crossLinkerType) {
      crossLinkerFormableBonds += functionalityPerType[crossLinkerType];
    } else {
      otherFormableBonds += functionalityPerType[type];
    }
  }

  // division by 2 is implicit
  return crossLinkerFormableBonds / (otherFormableBonds / strandLength);
}

/**
Compute the extent of reaction
(nr. of formed bonds in reaction / max. nr. of bonds formable)
NOTE: if your system has a non-integer number of possible bonds (e.g. one
site unbonded), this will not be rounded/respected in any way.

Arguments:
  - network: The poylmer network to do the computation for
  - functionalityPerType: a dictionary with key: type, and value:
functionality of this atom type. If None: will use max functionality per
type.

Returns:
  - p (float): The extent of reaction
*/
double
computeExtentOfReaction(pylimer_tools::entities::Universe network,
                        std::map<int, int> functionalityPerType)
{
  if (network.getNrOfAtoms() == 0) {
    return 1;
  }

  std::vector<int> allTypes = network.getPropertyValues<int>("type");

  int maxFormableBonds = 0;

  // count how often every type occurs
  for (int type : allTypes) {
    maxFormableBonds += functionalityPerType[type];
  }

  return network.getNrOfBonds() * 2.0 / static_cast<double>(maxFormableBonds);
}

double
computeExtentOfReaction(pylimer_tools::entities::Universe network)
{
  return computeExtentOfReaction(network,
                                 network.determineFunctionalityPerType());
}

/**
Compute the gelation point $p_{gel}$ as theoretically predicted
(gelation point = critical extent of reaction for gelation)

Source:
  - https://www.sciencedirect.com/science/article/pii/003238618990253X

Arguments:
  - r (double): The stoichiometric inbalance of reactants (see:
#computeStoichiometricInbalance)
  - f (int): functionality of the crossLinkers
  - g (int): functionality of the precursor polymer

Returns:
  - p_gel: critical extent of reaction for gelation
*/
double
predictGelationPoint(const double r, const int f, const int g = 2)
{
  return 1 / (r * (static_cast<double>(f) - 1) * (static_cast<double>(g) - 1));
}
}

#endif
