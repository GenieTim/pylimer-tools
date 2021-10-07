import igraph

class GraphDecorator(object):

    underlying_graph: igraph.Graph

    def getUnderlyingGraph(self):
        """
        Return the underlying graph

        Returns:
          - self.underlying_graph (igraph.Graph): the igraph.Graph object
        """
        return self.underlying_graph
