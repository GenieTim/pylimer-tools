from __future__ import annotations

import igraph
import numpy as np
import pandas as pd
from pylimer_tools.entities._graphDecorator import GraphDecorator
from pylimer_tools.entities.atom import Atom
from pylimer_tools.entities.molecule import Molecule


class Universum(GraphDecorator):

    boxSizes: list

    def __init__(self, boxSizes: list):
        """
        Instantiate this Universe (Collection of Molecules)

        Arguments:
          - boxSizes: a list containing the box lengths (x, y, z)

        Returns:
          - self (pylimer_tools.entities.Universum): the new Universum object
        """
        self.underlying_graph = igraph.Graph(directed=False)
        self.boxSizes = boxSizes

    def addAtomBondData(self, atomData: pd.DataFrame, bondData: pd.DataFrame) -> Universum:
        """
        Add atoms and bonds to the underlying graph.

        Arguments:
          - atomData: the dataframe containing atoms with their positions, id, type etc.
          - bondData: the dataframe containing two columns: one indicating where the bond originates, one where it goes. Direction irrelevant.

        Returns:
          - self: the Universum object for a fluent interface.
        """
        assert("id" in atomData.columns)
        assert("type" in atomData.columns)
        # first, create atoms and add them to the graph
        for index, row in atomData.iterrows():
            self.underlying_graph.add_vertices(
                ["atom{}".format(row['id'])],
                {
                    "type": row["type"],
                    "atom": Atom(row, boxSizes=self.boxSizes)
                })
        # then, follow up with the bonds.
        assert(len(bondData.columns) == 2)
        bondNames = bondData.applymap(lambda x: "atom{}".format(x))
        bondsArray = bondNames.to_numpy()
        self.underlying_graph.add_edges(bondsArray)
        self.underlying_graph.simplify()

        return self

    def getMolecules(self, ignoreAtomType=None) -> list[Molecule]:
        """
        Decompose the Universe into molecules, which could be either chains, networks, or even lonely atoms.

        Arguments:
          - ignoreAtomType: the atom type to ignore/omit from the molecules

        Returns:
          - molecules (list): a list of Molecule objects
        """
        molecules = []
        subgraphs = self.underlying_graph.decompose()

        for subgraph in subgraphs:
            if (ignoreAtomType is None):
                molecules.append(Molecule(subgraph))
            else:
                moleculeToSplit = Molecule(subgraph)
                molecules.extend(
                    moleculeToSplit.decomposeFurther(ignoreAtomType))

        return molecules

    def getAtom(self, atomId: int) -> Atom:
        """
        Find an atom by its ID

        Arguments:
          - atomId: the ID of the atom

        Returns:
          - atom (pylimer_tools.entities.Atom): the Atom object or None if it is not found
        """
        try:
            vertex = self.underlying_graph.vs.find(
                name="atom{}".format(atomId))
            if (vertex is not None):
                return vertex["atom"]
        except (ValueError, KeyError):
            pass
        return None

    def getAtomsWithType(self, atomType) -> list[Atom]:
        """
        Find an atom by its type

        Arguments:
          - atomType: the type of the atom

        Returns:
          - atoms (list<pylimer_tools.entities.Atom>): the Atom objects or None if it is not found
        """
        try:
            vertices = self.underlying_graph.vs.select(type_eq=atomType)
            if (vertices is not None):
                return [v["atom"] for v in vertices]
        except (ValueError, KeyError):
            pass
        return None

    def getVolume(self):
        """
        Get this object's volume

        Returns:
          - volume (float): the volume of the box
        """
        return np.prod(self.boxSizes)

    def getSize(self):
        """
        Get the number of atoms in this universe

        Returns:
          - nr (int): the number of atoms (nodes)
        """
        return self.underlying_graph.vcount()

    def setBoxSizes(self, boxSizes: list) -> Universum:
        """
        Re-set this Universe's size.

        Arguments:
          - boxSizes: a list containing the box lengths (x, y, z)

        Returns:
          - self (pylimer_tools.entities.Universum): the Universum object for a fluent interface.
        """
        self.boxSizes = boxSizes
        return self

    def reset(self) -> Universum:
        """
        Reset this Universe to be empty again.

        Returns:
          - self (pylimer_tools.entities.Universum): the Universum object for a fluent interface.
        """
        self.underlying_graph = igraph.Graph(directed=False)
        return self
