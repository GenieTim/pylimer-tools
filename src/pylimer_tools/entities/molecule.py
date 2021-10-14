
import warnings
from enum import Enum
from typing import Iterable

import igraph
import numpy as np
from pylimer_tools.entities._graphDecorator import GraphDecorator
from pylimer_tools.entities.atom import Atom


class Molecule(GraphDecorator, Iterable):

    class MoleculeType(Enum):
        UNDEFINED = 0
        NETWORK_STRAND = 1
        PRIMARY_LOOP = 2
        DANGLING_CHAIN = 3
        FREE_CHAIN = 4

    moleculeType = MoleculeType.UNDEFINED

    def __init__(self, molecule_graph, chainType=MoleculeType.UNDEFINED):
        self.underlying_graph = molecule_graph
        self.underlying_graph.simplify()
        self.moleculeType = chainType

    def decomposeFurther(self, splitAtomType):
        """
        Split this molecule into smaller molecules by ignoring all atoms with a given type.

        Arguments:
          - splitAtomType: the type of the atom to omit

        Returns:
          - subMolecules (list): a list of Molecule objects        
        """
        verticesToRemove = self.underlying_graph.vs.select(
            type_eq=splitAtomType)
        graphCopy = self.underlying_graph.copy()
        verticesToRemoveIndices = [v.index for v in verticesToRemove]
        graphCopy.delete_vertices(verticesToRemoveIndices)

        subMolecules = []
        subgraphs = graphCopy.decompose()
        for subgraph in subgraphs:
            subMolecules.append(
                Molecule(subgraph, Molecule.MoleculeType.FREE_CHAIN))
        return subMolecules

    def computeEndToEndDistance(self):
        """
        Compute the end-to-end distance of this molecule/chain.

        Returns:
          - $R_{ee}$ (float): the end-to-end distance, 
              `None` if the molecule/chain does not have two distinct ends.
        """
        if (self.underlying_graph.vcount() < 2):
            # warnings.warn("Underlying_graph too small for end-to-end distance")
            return None
        endNodes = self.underlying_graph.vs.select(_degree_eq=1)
        if (len(endNodes) != 2):
            # warnings.warn("Underlying_graph has {} instead of 2 end nodes".format(len(endNodes)))
            return None

        atom1 = endNodes[0]["atom"]
        atom2 = endNodes[1]["atom"]
        assert(isinstance(atom1, Atom) and isinstance(atom2, Atom))
        return atom1.computeDistanceTo(atom2)

    def computeBondLengths(self):
        """
        Calculate the bond lengths

        Returns:
          - a np.array of all bond lengths
        """
        lastAtom = None
        Rs = []
        for atom in self:
            if (lastAtom is not None):
                Rs.append(atom.computeDistanceTo(lastAtom))
            lastAtom = atom

        return np.array(Rs)

    def getLength(self):
        """
        Query the length of the molecule/chain (the nr. of atoms)

        Returns:
          - len (int): the nr. of nodes (atoms) in this molecule
        """
        return self.underlying_graph.vcount()

    def getType(self) -> MoleculeType:
        """
        Query the type of the molecule/chain (loop, dangling, etc.)

        Returns:
          - moleculeType (MoleculeType): one of the enum's Molecule.MoleculeType...
        """
        return self.moleculeType

    def __iter__(self):
        """
        Decorator of the depth-first iterator of the underlying graph.

        Yields:
          - atom (Atom): the atom of the current iteration.
        """
        endNodes = self.underlying_graph.vs.select(_degree_eq=1)
        if (self.getLength() < 1):
            return
        if (len(endNodes) < 1 and self.getLength() == 1):
            endNodes = self.underlying_graph.vs.select(_degree_eq=0)
        subIterator = self.underlying_graph.dfsiter(
            endNodes[0], mode="all", advanced=False)
        for node in subIterator:
            yield node["atom"]
