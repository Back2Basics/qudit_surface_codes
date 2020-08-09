import matplotlib.pyplot as plt
import numpy as np
from multiset import Multiset
from networkx import MultiGraph
from networkx import nx
from networkx.algorithms import bipartite
from sympy.combinatorics import Permutation
from typing import Tuple, List


class HybridSurfaceCodeGraph(MultiGraph):
    def __init__(self, sigma: Tuple[Tuple[int]], alpha: Tuple[Tuple[int]]):
        super().__init__()
        self.sigma = sigma  # should include singletons corresponding to fixed points
        self.alpha = alpha  # should include singletons corresponding to fixed points
        f = self.compute_phi()
        self.phi = self.permlist_to_tuple(f)
        self.build_node_info()  # print dictionary for [sigma, alpha, phi]
        self.node_dict = self.sigma_dict, self.alpha_dict, self.phi_dict

        self.node_info = ["sigma:", self.sigma_dict,
                          "alpha:", self.alpha_dict,
                          "phi:", self.phi_dict]

        self.code_graph = nx.MultiGraph()

        # Create black nodes for each cycle in sigma along with white nodes
        # representing "half edges" around the black nodes
        for cycle in self.sigma:
            self.code_graph.add_node(cycle, bipartite=1)
            for node in cycle:
                self.code_graph.add_node(node, bipartite=0)
                self.code_graph.add_edge(cycle, node)

        # Create black nodes for each cycle in phi along with white nodes
        # representing "half edges" around the black nodes
        for cycle in self.phi:
            self.code_graph.add_node(cycle, bipartite=1)
            for node in cycle:
                self.code_graph.add_edge(cycle, node)

        # Create nodes for each cycle in alpha then
        # glue the nodes corresponding to a the pairs
        for pair in self.alpha:
            self.code_graph.add_node(pair)
            self.code_graph = nx.contracted_nodes(self.code_graph, pair[0], pair[1], self_loops=True)
            # Now contract pair with pair[0] to make sure edges (white nodes) are labeled
            # by the pairs in alpha to keep track of the gluing from the previous step
            self.code_graph = nx.contracted_nodes(self.code_graph, pair, pair[0], self_loops=True)

        # Define the white and black nodes. White correspond to edges labeled by
        # cycles in alpha. Black correspond to nodes labeled by cycles in sigma
        # (vertices) and phi (faces)
        self.black_nodes, self.white_nodes = bipartite.sets(self.code_graph)
        self.vertex_basis()
        self.edge_basis()
        self.face_basis()
        self.d_2()
        self.d_1()

    def permlist_to_tuple(self, perms):
        """
        convert list of lists to tuple of tuples in order to have two level iterables
        that are hashable for the dictionaries used later
        """
        return tuple(tuple(perm) for perm in perms)

    def compute_phi(self):
        """compute the list of lists full cyclic form of phi (faces of dessin [sigma, alpha, phi])"""
        s = Permutation(self.sigma)
        a = Permutation(self.alpha)
        f = ~(a * s)
        f = f.full_cyclic_form  # prints permutation as a list of lists including all singletons (fixed points)
        return f

    def build_node_info(self):
        count = -1

        self.sigma_dict = dict()
        for count, cycle in enumerate(self.sigma):
            self.sigma_dict[cycle] = count

        self.phi_dict = dict()
        for count, cycle in enumerate(self.phi, start=count + 1):
            self.phi_dict[cycle] = count

        self.alpha_dict = dict()
        for count, pair in enumerate(self.alpha, start=count + 1):
            self.alpha_dict[pair] = count
        return tuple([self.sigma_dict, self.alpha_dict, self.phi_dict])

    def boundary_1(self, edge):
        """
        compute boundary of a single edge given by a white node (cycle in alpha)
        """
        boundary1 = [node for node in self.code_graph.neighbors(edge) if node in self.sigma_dict]
        return boundary1

    def del_1(self, edges: List[Tuple[int]]):
        """
        boundary of a list of edges, i.e. an arbitrary 1-chain over Z/2Z
        """
        boundary_list = Multiset()
        for edge in edges:
            boundary_list.update(Multiset(self.boundary_1(edge)))
        return boundary_list

    def boundary_2(self, face):
        """
        compute boundary of a single face
        """
        boundary2 = list(self.code_graph.neighbors(face))
        return boundary2

    def del_2(self, faces: List[Tuple[int]]):
        """
        boundary of a list of faces, i.e. an arbitrary 2-chain over Z/2Z
        """
        boundary_list = Multiset()
        for face in faces:
            boundary_list.update(Multiset(self.boundary_2(face)))
        return boundary_list

    def coboundary_1(self, star):
        """
        compute coboundary of a single star
        """
        coboundary = self.code_graph.neighbors(star)
        return coboundary

    def delta_1(self, stars: List[Tuple[int]]):
        """
        coboundary of a list of stars, i.e. an arbitrary 0-cochain over Z/2Z
        """
        coboundary_list = Multiset()
        for star in stars:
            coboundary_list.update(Multiset(self.boundary_1(star)))
        return coboundary_list

    def coboundary_2(self, edge):
        """
        compute coboundary of a single edge given by a white node (cycle in alpha)
        """
        coboundary2 = [node for node in self.code_graph.neighbors(edge) if node in self.phi_dict]
        return coboundary2

    def delta_2(self, edges: List[Tuple[int]]):
        """
        coboundary of a list of edges, i.e. an arbitrary 1-cochain over Z/2Z
        given by a list of cycles in alpha
        """
        coboundary_list = Multiset()
        for edge in edges:
            coboundary_list.update(Multiset(self.boundary_1(edge)))
        return coboundary_list

    def euler_characteristic(self):
        """
        Compute the Euler characteristic of the surface in which the graph is embedded
        """
        chi = len(self.phi) - len(self.alpha) + len(self.sigma)
        return (chi)

    def genus(self):
        """
        Compute the genus of the surface in which the graph is embedded
        """
        g = int(-(len(self.phi) - len(self.alpha) + len(self.sigma) - 2) / 2)
        return (g)

    def vertex_basis(self):
        self.v_basis_dict = dict()
        self.v_dict = dict()
        A = np.eye(len(self.sigma), dtype=np.uint8)
        for count, cycle in enumerate(self.sigma):
            self.v_dict[cycle] = count
            self.v_basis_dict[cycle] = A[:, count]

    def edge_basis(self):
        self.e_basis_dict = dict()
        self.e_dict = dict()
        B = np.eye(len(self.alpha), dtype=np.uint8)
        for count, cycle in enumerate(self.alpha):
            self.e_dict[cycle] = count
            self.e_basis_dict[cycle] = B[count, :].reshape(1, -1)

    def face_basis(self):
        self.f_basis_dict = dict()
        self.f_dict = dict()
        C = np.eye(len(self.phi), dtype=np.uint8)
        for count, cycle in enumerate(self.phi):
            self.f_dict[cycle] = count
            self.f_basis_dict[cycle] = C[count, :].reshape(1, -1)

    def d_2(self):
        """
        Create the matrix D2 for the differential d_2 in the chain complex
        that maps faces to their boundary edges. This can be used in the homology
        class as one of the input matrices.
        :return: matrix D2 of the differential d_2
        """
        self.D2 = np.zeros(len(self.e_basis_dict), dtype=np.uint8)
        for cycle in self.phi:
            bd = self.boundary_2(cycle)
            if bd == []:
                image = np.zeros(len(self.e_basis_dict))
            else:
                image = sum([self.e_basis_dict[edge] for edge in bd])
            self.D2 = np.vstack((self.D2, image))
        self.D2 = np.array(self.D2[1:, :]).T

    def d_1(self):
        """
        Create the matrix D1 for the differential d_1 in the chain complex
        that maps edges to their boundary vertices. This can be used in the homology
        class as one of the input matrices
        :return: matrix D1 of the differential d_1
        """
        self.D1 = np.zeros(len(self.v_basis_dict), dtype=np.uint8)
        for cycle in self.alpha:
            bd = self.boundary_1(cycle)
            if bd == []:
                image = np.zeros(len(self.v_basis_dict))
            else:
                image = sum([self.v_basis_dict[vertex] for vertex in bd])
            self.D1 = np.vstack((self.D1, image))
        self.D1 = np.array(self.D1[1:, :]).T

    def draw(self, node_type='', layout=''):
        """
        Draw graph with vertices, edges, and faces labeled by colored nodes and their integer indices
        corresponding to the qubit indices for the surface code
        """
        if not node_type in ['cycles', 'dict']:
            raise ValueError('node_type can be "cycles" or "dict"')

        if layout == 'spring':
            pos = nx.spring_layout(self.code_graph)
        if layout == 'spectral':
            pos = nx.spectral_layout(self.code_graph)
        if layout == 'planar':
            pos = nx.planar_layout(self.code_graph)
        if layout == 'shell':
            pos = nx.shell_layout(self.code_graph)
        if layout == 'circular':
            pos = nx.circular_layout(self.code_graph)
        if layout == 'spiral':
            pos = nx.spiral_layout(self.code_graph)
        if layout == 'random':
            pos = nx.random_layout(self.code_graph)
        # white nodes
        nx.draw_networkx_nodes(self.code_graph, pos,
                               nodelist=list(self.alpha),
                               node_color='c',
                               node_size=500,
                               alpha=0.3)
        # vertex nodes
        nx.draw_networkx_nodes(self.code_graph, pos,
                               nodelist=list(self.sigma),
                               node_color='b',
                               node_size=500,
                               alpha=0.6)
        # face nodes
        nx.draw_networkx_nodes(self.code_graph, pos,
                               nodelist=list(self.phi),
                               node_color='r',
                               node_size=500,
                               alpha=0.6)
        # edges
        nx.draw_networkx_edges(self.code_graph, pos, width=1.0, alpha=0.5)

        labels = {}

        if node_type == 'cycles':
            '''
            label nodes the cycles of sigma, alpha, and phi
            '''
            for node in self.alpha_dict:
                # stuff = self.alpha_dict[node]
                labels[node] = f'$e$({node})'
            for node in self.sigma_dict:
                # something = self.sigma_dict[node]
                labels[node] = f'$v$({node})'
            for node in self.phi_dict:
                # something2 = self.phi_dict[node]
                labels[node] = f'$f$({node})'
            nx.draw_networkx_labels(self.code_graph, pos, labels, font_size=12)

        if node_type == 'dict':
            '''
            label nodes with v, e, f and indices given by node_dict corresponding to
            qubit indices of surface code
            '''

            for node in self.alpha_dict:
                # stuff = self.alpha_dict[node]
                labels[node] = f'$e$({self.alpha_dict[node]})'
            for node in self.sigma_dict:
                # something = self.sigma_dict[node]
                labels[node] = f'$v$({self.sigma_dict[node]})'
            for node in self.phi_dict:
                # something2 = self.phi_dict[node]
                labels[node] = f'$f$({self.phi_dict[node]})'
            nx.draw_networkx_labels(self.code_graph, pos, labels, font_size=12)

        # plt.axis('off')
        # plt.savefig("labels_and_colors.png") # save as png
        plt.show()  # display
