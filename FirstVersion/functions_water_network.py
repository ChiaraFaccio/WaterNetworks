# Chiara Faccio, Scuola Normale Superiore, Pisa
# chiara.faccio@sns.it
# May 2022

import numpy as np
import numpy.linalg
from numpy.linalg import norm
from numpy import ones, zeros, arange, inf
import networkx as nx
import sys, os, math
from math import exp, log
from tkinter import Tk, filedialog
import matplotlib.pyplot as plt
import scipy
import scipy.io
from scipy.sparse import csr_matrix, find
from scipy.spatial import distance_matrix
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply
from scipy.sparse import identity, diags, bmat, kron
from mpl_toolkits.mplot3d import Axes3D
import tracemalloc
import time
from scipy.sparse.linalg import spsolve
from scipy.spatial.distance import pdist



def exponential_quadrature(A, u, v, tol=1e-7, maxit=50):

    """
    
    This functions is from the Python package NetworkSNS  https://github.com/alb95/NetworkSNS/blob/main/networksns/centrality_measures.py of Alberto Bucci
    
    BSD 2-Clause License

    Copyright (c) 2021, Alberto Bucci
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    1. Redistributions of source code must retain the above copyright notice, this
       list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above copyright notice,
       this list of conditions and the following disclaimer in the documentation
       and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
    FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.   
    
    
    Author: Alberto Bucci
    
    Computes :math:`q=u^T e^A v`. For the computation polarization rule and Lanczos iteration are used [1]_.

    Parameters
    __________
    A: array_like
        sparse/dense symmetric matrix
    u: array
        vector
    v: array
        vector
    tol: float,optional
        tolerance for convergence, relative accuracy, default: 1e^-7.
    maxit: integer, optional
        maximum number of Lanczos iterations, default: 50.


    :return: **q**: (float)
        value of the quadrature :math:`u^Te^Av`.

    Examples
    ________
    .. code:: python

     >>>  from sobigdatainit.sobigdata_pkg import sobigdata as sbd
     >>>  import numpy as np

    Create symmetric matrix :math:`A`

    .. code:: python

     >>>    A = np.arange(0, 9, 1)
     >>>    A = A.reshape(3, 3)
     >>>    A = A + A.transpose()
            array([[ 0,  4,  8],
                   [ 4,  8, 12],
                   [ 8, 12, 16]])

    Create vectors :math:`u` and :math:`v`

    .. code:: python

     >>>    u = np.arange(0, 3)
            array([0, 1, 2])
     >>>    v = np.array([2,5,2])
            array([2, 5, 2])

    Compute :math:`q=u^T e^A v`

     >>>    q = sbd.exponential_quadrature(A, u, v)

    References
    ----------
    .. [1] G. H. Golub, and G. Meurant (2010)
           "Matrices, Moments and Quadrature with Applications",
           Princeton University Press, Princeton, NJ.
    """

    quadrature_sum = exponential_symmetric_quadrature(A, u+v, tol, maxit)
    quadrature_difference = exponential_symmetric_quadrature(A, u-v, tol, maxit)

    # polarization formula
    quadrature = (quadrature_sum - quadrature_difference)/4
    return quadrature


def exponential_symmetric_quadrature(A, u, tol=1e-7, maxit=50):

    """
    
    This functions is from the Python package NetworkSNS  https://github.com/alb95/NetworkSNS/blob/main/networksns/centrality_measures.py of Alberto Bucci
    
    BSD 2-Clause License

    Copyright (c) 2021, Alberto Bucci
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    1. Redistributions of source code must retain the above copyright notice, this
       list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above copyright notice,
       this list of conditions and the following disclaimer in the documentation
       and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
    FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    
    
    Author : Alberto Bucci
    
    Computes :math:`q=u^Te^Au`.
    The computation is done by means of Lanczos method according to [1]_.
    Parameters
    __________
    A: array_like
        sparse/dense symmetric matrix.
    u: array
        vector.
    tol: float,optional
        tolerance for convergence, relative accuracy, default: 1e-7.
    maxit: integer, optional
     maximum number of Lanczos iterations, default: 50.
    :return: **q**: (float)
     value of the quadratic form :math:`u^Te^Au`.
    Examples
    ________
    .. code:: python
     >>>  from networksns import centrality_measures as cm
     >>>  import numpy as np
    Create symmetric matrix :math:`A`
    .. code:: python
     >>>    A = np.arange(0, 9, 1)
     >>>    A = A.reshape(3, 3)
     >>>    A = A + A.transpose()
            array([[ 0,  4,  8],
                   [ 4,  8, 12],
                   [ 8, 12, 16]])
    Create vector :math:`u`
        .. code:: python
    >>>     u = np.arange(0, 3)
            array([0, 1, 2])
    Compute :math:`q=u^T e^A u`.
     >>>    q = cm.exponential_symmetric_quadrature(A, u)
    References
    ----------
    .. [1] G. H. Golub, and G. Meurant (2010)
           "Matrices, Moments and Quadrature with Applications",
           Princeton University Press, Princeton, NJ.
    """
    quadrature = 1
    quadrature_old = 2
    old_err = 1
    err = 1
    omega = []
    gamma = []
    u_norm = norm(u)
    if u_norm == 0:
        return 0
    else:
        x_0 = u/u_norm
        #  computing Lanczos matrix J
        omega.append(x_0.dot(A.dot(x_0)))
        r = A.dot(x_0) - omega[0] * x_0
        r_norm = norm(r)
        if r_norm == 0:
            return exp(omega[0]) * u_norm ** 2
        gamma.append(r_norm)
        x_1 = r / r_norm
        it = 1  # iterations
        while err > tol and old_err > tol and it < maxit:
            z = A.dot(x_1) - gamma[it - 1] * x_0
            omega.append(x_1.dot(z))
            x_2 = z - omega[it] * x_1
            if norm(x_2) == 0:
                gamma.append(0)  # variable used only to avoid deletion in line: eta = gamma[:-1]
                eta = gamma[:-1]
                J = diags(omega)
                J = J + diags(eta, 1)
                J = J + diags(eta, -1)
                e_1 = np.zeros(len(omega))
                e_1[0] = 1
                quadrature = e_1.dot(expm_multiply(J, e_1)) * (norm(u)) ** 2
                break
            gamma.append(norm(x_2))
            x_0 = x_1
            x_1 = x_2 / gamma[it]
            eta = gamma[:-1]
            J = diags(omega)
            J = J + diags(eta, 1)
            J = J + diags(eta, -1)
            e_1 = np.zeros(len(omega))
            e_1[0] = 1
            quadrature_very_old = quadrature_old
            quadrature_old = quadrature
            quadrature = e_1.dot(expm_multiply(J, e_1))*(norm(u))**2
            old_err = err
            err = max(abs((quadrature_old-quadrature))/quadrature_old,
                      abs((quadrature_very_old-quadrature_old))/quadrature_very_old)
            it = it+1
    return quadrature
    
    
def subgraph_centrality_Bucci(G, t=1, tol=1e-7, maxit=50):

    """
    
    This functions is from the Python package NetworkSNS  https://github.com/alb95/NetworkSNS/blob/main/networksns/centrality_measures.py of Alberto Bucci
    
    BSD 2-Clause License

    Copyright (c) 2021, Alberto Bucci
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    1. Redistributions of source code must retain the above copyright notice, this
       list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above copyright notice,
       this list of conditions and the following disclaimer in the documentation
       and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
    FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    
    
    Author : Alberto Bucci
    
    Computes the subgraph centrality of all the nodes in graph :math:`G`.
    The subgraph centrality of the :math:`i^{th}` node is given by :math:`[e^{tA}]_{ii}=e_i^T (e^{tA})e_i`,
    where :math:`e_i` and :math:`A` denote respectively the :math:`i^{th}` vector of the canonical basis and the adjacency matrix of the graph [1]_.
    Parameters
    __________
    G: Graph or DiGraph object
        a graph.
    t: scalar, optional
     when exponentiating multiply the adjacency matrix by t, default: 1.
    tol: float,optional
     tolerance for convergence, relative accuracy, default: 1e-7.
    maxit: integer, optional
     maximum number of Lanczos iterations, default: 50.
    :return: **sc** (dict)  subgraph centrality of all nodes in :math:`G`.
    Examples
    ________
    .. code:: python
     >>>  from networksns import centrality_measures as cm
     >>>  import networkx as nx
    Create graph :math:`G`.
    .. code:: python
     >>>    G = nx.Graph()
     >>>    G.add_edge(1, 2)
     >>>    G.add_edge(2, 3)
            EdgeView([(1, 2), (2, 3)])
    Compute the subgraph centrality.
     >>>    sc = cm.subgraph_centrality(G)
    References
    ----------
    .. [1] Ernesto Estrada and Juan A. Rodríguez-Velázquez (2005)
           Subgraph centrality in complex networks,
           Physical Review, Volume 71, Issue 5,
           https://doi.org/10.1103/PhysRevE.71.056103
    """
    
    n = G.number_of_nodes()
    node_list = list(G.nodes)
    subgraph_centrality = np.zeros(n)
    for i in range(n):
        subgraph_centrality[i] = node_subgraph_centrality(G, node_list[i], t, tol, maxit)
    centrality = dict(zip(node_list, subgraph_centrality))
    return centrality
    

def node_subgraph_centrality(G, u, t=1, tol=1e-7, maxit=50):

    """
    
    This functions is from the Python package NetworkSNS  https://github.com/alb95/NetworkSNS/blob/main/networksns/centrality_measures.py of Alberto Bucci
    
    BSD 2-Clause License

    Copyright (c) 2021, Alberto Bucci
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    1. Redistributions of source code must retain the above copyright notice, this
       list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above copyright notice,
       this list of conditions and the following disclaimer in the documentation
       and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
    FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    
    
    Author : Alberto Bucci
    
    Computes the subgraph centrality of node :math:`u`.
    If node :math:`u` is the :math:`i^{th}` node of the graph, the subgraph centrality of node :math:`u` is given by :math:`[e^{tA}]_{ii}=e_i^T (e^{tA})e_i`,
    where :math:`e_i` and :math:`A` denote respectively the :math:`i^{th}` vector of the canonical basis and the adjacency matrix of the graph [1]_.
    Parameters
    __________
    G: Graph object
        a graph.
    u: node_id
        node in G.
    t: scalar, optional
     when exponentiating multiply the adjacency matrix by t, default: 1.
    tol: float,optional
     tolerance for convergence, relative accuracy, default: 1e-7.
    maxit: integer, optional
     maximum number of Lanczos iterations, default: 50.
     :return: **sc_u** (float) subgraph centrality of node :math:`u`.
    Examples
    ________
    .. code:: python
     >>>  from networksns import centrality_measures as cm
     >>>  import networkx as nx
    Create graph :math:`G`.
    .. code:: python
     >>>    G = nx.Graph()
     >>>    G.add_edge(1, 'u')
     >>>    G.add_edge('u', 2)
            EdgeView([(1, 'u'), ('u', 2)])
    Compute the node total communicability
     >>>    sc_u = cm.node_subgraph_centrality(G, 'u')
    References
    ----------
    .. [1] Ernesto Estrada and Juan A. Rodríguez-Velázquez (2005)
           Subgraph centrality in complex networks,
           Physical Review, Volume 71, Issue 5,
           https://doi.org/10.1103/PhysRevE.71.056103
    """
    
    n = G.number_of_nodes()
    node_list = G.nodes
    enumerated_nodes = dict(zip(node_list, np.arange(n)))
    node_position = enumerated_nodes[u]
    e_node = np.zeros(n)
    e_node[node_position] = 1
    Adj = nx.adjacency_matrix(G)
    if t != 1:
        Adj = Adj * t
    subgraph_centrality = exponential_symmetric_quadrature(Adj, e_node, tol, maxit)
    return subgraph_centrality


def total_communicability(G, t=None):

    """
    
    This functions is from the Python package NetworkSNS  https://github.com/alb95/NetworkSNS/blob/main/networksns/centrality_measures.py of Alberto Bucci
    
    BSD 2-Clause License

    Copyright (c) 2021, Alberto Bucci
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    1. Redistributions of source code must retain the above copyright notice, this
       list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above copyright notice,
       this list of conditions and the following disclaimer in the documentation
       and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
    FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    
    
    Author : Alberto Bucci
    
    Computes the total communicability of :math:`G`.

    Total communicability is defined as the row sum of the exponential of the adjacency matrix [1]_, so denoting with
    :math:`A` the adjacency matrix of graph :math:`G` and with :math:`\\mathbf{1}` the vector of all ones,
    we have :math:`tc= e^A \\mathbf{1}`.

    Parameters
    __________
    G: Graph or DiGraph object
        a graph.
    t: scalar, optional
        exponentiate multiply the adjacency matrix by :math:`t`, default None.


    :return: **tc** (dict)  total communicability of graph :math:`G`.

    Examples
    ________
    .. code:: python

     >>>  from sobigdatainit.sobigdata_pkg import sobigdata as sbd
     >>>  import networkx as nx

    Create graph :math:`G`

    .. code:: python

     >>>    G = nx.Graph()
     >>>    G.add_edge(1, 2)
     >>>    G.add_edge(2, 3)
            EdgeView([(1, 2), (2, 3)])

    Compute :math:`tc= e^A \\mathbf{1}`

     >>>    tc = sbd.total_communicability(G)


    References
    ----------
    .. [1] Michele Benzi, Christine Klymko (2013)
           Total communicability as a centrality measure,
           Journal of Complex Networks, Volume 1, Issue 2, Pages 124–149,
           https://doi.org/10.1093/comnet/cnt007
    """

    n = G.number_of_nodes()
    node_list = G.nodes
    e = ones(n)  # vector of all ones
    Adj = nx.adjacency_matrix(G)
    if t != 1:
        Adj = Adj*t
    tot_communicability = expm_multiply(Adj, e)
    centrality = dict(zip(node_list, tot_communicability))
    return centrality


def total_directed_communicability(G, t=None):

    """
    
    This functions is from the Python package NetworkSNS  https://github.com/alb95/NetworkSNS/blob/main/networksns/centrality_measures.py of Alberto Bucci
    
    BSD 2-Clause License

    Copyright (c) 2021, Alberto Bucci
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    1. Redistributions of source code must retain the above copyright notice, this
       list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above copyright notice,
       this list of conditions and the following disclaimer in the documentation
       and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
    FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.   
    
    
    Author: Alberto Bucci
    
    Computes broadcast and receive communicability of a directed graph :math:`G`.

    Let :math:`A` be the adjacency matrix of :math:`G`, then\
     :math:`\\mathcal{A}=\\begin{pmatrix} 0 & A \\\\ A^T & 0 \\end{pmatrix}`
    is the adjacency matrix of the associated undirected bipartite graph [1]_.

    Broadcast communicability  (hub centrality) is given by :math:`bc = (e^{\\mathcal{A}}\\mathbf{1})_{1:n}`.


    Receive communicability (authority centrality) is given by :math:`rc = (e^{\\mathcal{A}}\\mathbf{1})_{n+1:2n}`.

    Parameters
    __________

    G: DiGraph object
        a directed graph.
    t: scalar, optional
     when exponentiate multiply the adjacency matrix by t, default None.

    Returns
    ________

    bc: dict
     broadcast communicability
    rc: dict
     receive communicability

    Examples
    ________

    .. code:: python

     >>>  from sobigdatainit.sobigdata_pkg import sobigdata as sbd
     >>>  import networkx as nx

    Create graph :math:`G`

    .. code:: python

     >>>    G = nx.DiGraph()
     >>>    G.add_edge(1, 2)
     >>>    G.add_edge(1, 3)
     >>>    G.add_edge(2, 3)
     >>>    G.add_edge(3, 1)
            OutEdgeView([(1, 2), (1, 3), (2, 3), (3, 1)])

    Compute broadcast and receive communicability

     >>>    bc, rc = sbd.total_directed_communicability(G)

    References
    ----------
    .. [1] Michele Benzi, Paola Boito (2020),
           Matrix functions in network analysis,
           GAMM‐Mitteilungen, Volume 43, Issue 3,
           https://onlinelibrary.wiley.com/doi/abs/10.1002/gamm.202000012
    """

    n = G.number_of_nodes()
    node_list = G.nodes
    e = ones(2*n)  # vector of all ones
    Adj = nx.adjacency_matrix(G, nodelist = range(0,n))
    Bip_Adj = bmat([[None, Adj], [Adj.transpose(), None]])
    if t is not None:
        Bip_Adj = Bip_Adj*t
    tot_communicability = expm_multiply(Bip_Adj, e)
    bc = dict(list(zip(node_list, tot_communicability[0:n])))  # hub centrality
    rc = dict(list(zip(node_list, tot_communicability[n:2*n])))  # authority centrality
    return bc, rc


def total_network_communicability(G, t=None, tol=1e-7, maxit=50):

    """
    
    This functions is from the Python package NetworkSNS  https://github.com/alb95/NetworkSNS/blob/main/networksns/centrality_measures.py of Alberto Bucci
    
    BSD 2-Clause License

    Copyright (c) 2021, Alberto Bucci
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    1. Redistributions of source code must retain the above copyright notice, this
       list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above copyright notice,
       this list of conditions and the following disclaimer in the documentation
       and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
    FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.   
    
    Author : Alberto Bucci
    
    Computes the total network communicability of :math:`G`.

    Total network communicability is defined as the sum over all elements of the exponential of the adjacency matrix\
     [1]_, so denoting with
    :math:`A` the adjacency matrix of graph :math:`G` and with :math:`\\mathbf{1}` the vector of all ones,
    we have :math:`tnc = \\mathbf{1}^T e^A \\mathbf{1}`.

    Parameters
    __________
    G: Graph object
        an undirected graph.
    t: scalar, optional
        exponentiate multiply the adjacency matrix by t, default None.
    tol: float,optional
     tolerance for convergence, relative accuracy, default: 1e^7.
    maxit: integer, optional
     maximum number of Lanczos iterations, default: 50


    :return: **tnc** (float)  total network communicability of graph :math:`G`.

    Examples
    ________
    .. code:: python

     >>>  from sobigdatainit.sobigdata_pkg import sobigdata as sbd
     >>>  import networkx as nx

    Create graph :math:`G`

    .. code:: python

     >>>    G = nx.Graph()
     >>>    G.add_edge(1, 'u')
     >>>    G.add_edge('u', 2)
            EdgeView([(1, 'u'), ('u', 2)])

    Compute the total network communicability

     >>>    tnc = sbd.total_network_communicability(G)

    References
    ----------
    .. [1] Michele Benzi, Christine Klymko (2013)
           Total communicability as a centrality measure,
           Journal of Complex Networks, Volume 1, Issue 2, Pages 124–149,
           https://doi.org/10.1093/comnet/cnt007
    """

    n = G.number_of_nodes()
    Adj = nx.adjacency_matrix(G, nodelist = range(0,n))
    if t is not None:
        Adj = Adj*t
    e = ones(n)
    tot_net_communicability = exponential_symmetric_quadrature(Adj, e, tol, maxit)
    return tot_net_communicability


def divide_file():

    '''
    Splits the file.gro into the various frames, and saves in the folder "frames"
    '''
    
    print("Choose the file.gro : ")
    root = Tk()
    root.withdraw()
    fname = filedialog.askopenfilename()
    root.destroy()

    print("You have chosen the file {} ".format(fname))
    name_file = os.path.splitext(os.path.basename(fname))[0]

    try: 
        
        smallfile = None
        path = sys.path[0]
        path = os.path.join(path, str(name_file) )
        
        try:
            os.mkdir(path)
        except FileExistsError:
            print('Attention, this file has already been split and its folder already exists.')
            sys.exit(1)
            
        path2 = os.path.join(path, 'frames' )
        os.mkdir(path2)
        
        with open(fname, 'r') as bigfile:
            bigfile.seek(0)
            
            bigfile.readline()
            n = bigfile.readline()
            n = int(int(n[1:]))
        
        lines_per_file = n + 3        
        
        with open(fname) as bigfile:
            bigfile.seek(0)
            for lineno, line in enumerate(bigfile):
                if lineno % lines_per_file == 0:
                    if smallfile:
                        smallfile.close()
                    small_filename = '{}.gro'.format(lineno + lines_per_file)
                    smallfile = open(os.path.join(path2, small_filename), "w")
                smallfile.write(line)
            if smallfile:
                smallfile.close()
    except FileNotFoundError:
        print('ERROR: File {} is not present'.format(fname))
        sys.exit(1)
        
    return path2, name_file
    
    
def extract_number_nodes(file, name_node):
    
    '''
    Returns the number of nodes in the central box
    '''
    
    count = 0
    with open(file, 'r') as data:
        data.seek(0)    
        for line in data:
            if name_node in line:
                count += 1
                
    return count


def extract_data(name_atom, path, file):

    '''
    Extracts the coordinates of 'name_atom' and the values of the box for the frame 'file'
    '''
    
    try: 
        # estraggo in una matrice coord le coordinate degli atomi che mi interessano, e il vettore del box 
        coord = []
        count = 0
        with open(os.path.join(path, file), 'r') as data:
            data.seek(0)    
            for line in data:
                if name_atom not in line:
                    count += 1
                    if str('water') not in line:
                        if count > 2:
                            try:
                                box = [float(x) for x in line.split()]
                            except:
                                continue
                    continue
                else:
                    count += 1
                    xyz = np.array((float(line[20:28]), float(line[28:36]), float(line[36:44])))
                    coord.append(xyz)
        coord = np.asmatrix(coord)
    except FileNotFoundError:
        print('ERROR: File {} is not present or name atom {} is wrong'.format(file, name_atom))
        sys.exit(1)
        
    return {'coord': coord, 'box': box}
   
   
def impose_boundary_condition_27_box_given_cood(coord, box):

    '''
    Imposes the periodic boundary conditions given the coordinates 
    '''
    
    box_plus_x = coord.copy()
    box_plus_x[:,0] = box_plus_x[:,0] + box[0]
    
    box_plus_x_plus_y = coord.copy()
    box_plus_x_plus_y[:,0] = box_plus_x_plus_y[:,0] + box[0]
    box_plus_x_plus_y[:,1] = box_plus_x_plus_y[:,1] + box[1]
    
    box_plus_x_plus_z = coord.copy()
    box_plus_x_plus_z[:,0] = box_plus_x_plus_z[:,0] + box[0]
    box_plus_x_plus_z[:,2] = box_plus_x_plus_z[:,2] + box[2]
    
    box_plus_x_minus_y = coord.copy()
    box_plus_x_minus_y[:,0] = box_plus_x_minus_y[:,0] + box[0]
    box_plus_x_minus_y[:,1] = box_plus_x_minus_y[:,1] - box[1]
    
    box_plus_x_minus_z = coord.copy()
    box_plus_x_minus_z[:,0] = box_plus_x_minus_z[:,0] + box[0]
    box_plus_x_minus_z[:,2] = box_plus_x_minus_z[:,2] - box[2]
    
    box_plus_x_plus_y_plus_z = coord.copy()
    box_plus_x_plus_y_plus_z[:,0] = box_plus_x_plus_y_plus_z[:,0] + box[0]
    box_plus_x_plus_y_plus_z[:,1] = box_plus_x_plus_y_plus_z[:,1] + box[1]
    box_plus_x_plus_y_plus_z[:,2] = box_plus_x_plus_y_plus_z[:,2] + box[2]
    
    box_plus_x_minus_y_plus_z = coord.copy()
    box_plus_x_minus_y_plus_z[:,0] = box_plus_x_minus_y_plus_z[:,0] + box[0]
    box_plus_x_minus_y_plus_z[:,1] = box_plus_x_minus_y_plus_z[:,1] - box[1]
    box_plus_x_minus_y_plus_z[:,2] = box_plus_x_minus_y_plus_z[:,2] + box[2]
    
    box_plus_x_plus_y_minus_z = coord.copy()
    box_plus_x_plus_y_minus_z[:,0] = box_plus_x_plus_y_minus_z[:,0] + box[0]
    box_plus_x_plus_y_minus_z[:,1] = box_plus_x_plus_y_minus_z[:,1] + box[1]
    box_plus_x_plus_y_minus_z[:,2] = box_plus_x_plus_y_minus_z[:,2] - box[2]
    
    box_plus_x_minus_y_minus_z = coord.copy()
    box_plus_x_minus_y_minus_z[:,0] = box_plus_x_minus_y_minus_z[:,0] + box[0]
    box_plus_x_minus_y_minus_z[:,1] = box_plus_x_minus_y_minus_z[:,1] - box[1]
    box_plus_x_minus_y_minus_z[:,2] = box_plus_x_minus_y_minus_z[:,2] - box[2]   
    
    box_minus_x = coord.copy()
    box_minus_x[:,0] = box_minus_x[:,0] - box[0]
    
    box_minus_x_plus_y = coord.copy()
    box_minus_x_plus_y[:,0] = box_minus_x_plus_y[:,0] - box[0]
    box_minus_x_plus_y[:,1] = box_minus_x_plus_y[:,1] + box[1]
    
    box_minus_x_plus_z = coord.copy()
    box_minus_x_plus_z[:,0] = box_minus_x_plus_z[:,0] - box[0]
    box_minus_x_plus_z[:,2] = box_minus_x_plus_z[:,2] + box[2]
    
    box_minus_x_minus_y = coord.copy()
    box_minus_x_minus_y[:,0] = box_minus_x_minus_y[:,0] - box[0]
    box_minus_x_minus_y[:,1] = box_minus_x_minus_y[:,1] - box[1]
    
    box_minus_x_minus_z = coord.copy()
    box_minus_x_minus_z[:,0] = box_minus_x_minus_z[:,0] - box[0]
    box_minus_x_minus_z[:,2] = box_minus_x_minus_z[:,2] - box[2]
    
    box_minus_x_plus_y_plus_z = coord.copy()
    box_minus_x_plus_y_plus_z[:,0] = box_minus_x_plus_y_plus_z[:,0] - box[0]
    box_minus_x_plus_y_plus_z[:,1] = box_minus_x_plus_y_plus_z[:,1] + box[1]
    box_minus_x_plus_y_plus_z[:,2] = box_minus_x_plus_y_plus_z[:,2] + box[2]
    
    box_minus_x_minus_y_plus_z = coord.copy()
    box_minus_x_minus_y_plus_z[:,0] = box_minus_x_minus_y_plus_z[:,0] - box[0]
    box_minus_x_minus_y_plus_z[:,1] = box_minus_x_minus_y_plus_z[:,1] - box[1]
    box_minus_x_minus_y_plus_z[:,2] = box_minus_x_minus_y_plus_z[:,2] + box[2]
    
    box_minus_x_plus_y_minus_z = coord.copy()
    box_minus_x_plus_y_minus_z[:,0] = box_minus_x_plus_y_minus_z[:,0] - box[0]
    box_minus_x_plus_y_minus_z[:,1] = box_minus_x_plus_y_minus_z[:,1] + box[1]
    box_minus_x_plus_y_minus_z[:,2] = box_minus_x_plus_y_minus_z[:,2] - box[2]
    
    box_minus_x_minus_y_minus_z = coord.copy()
    box_minus_x_minus_y_minus_z[:,0] = box_minus_x_minus_y_minus_z[:,0] - box[0]
    box_minus_x_minus_y_minus_z[:,1] = box_minus_x_minus_y_minus_z[:,1] - box[1]
    box_minus_x_minus_y_minus_z[:,2] = box_minus_x_minus_y_minus_z[:,2] - box[2]     
    
    box_plus_y = coord.copy()
    box_plus_y[:,1] = box_plus_y[:,1] + box[1]
    
    box_plus_y_plus_z = coord.copy()
    box_plus_y_plus_z[:,1] = box_plus_y_plus_z[:,1] + box[1]
    box_plus_y_plus_z[:,2] = box_plus_y_plus_z[:,2] + box[2]
    
    box_plus_y_minus_z = coord.copy()
    box_plus_y_minus_z[:,1] = box_plus_y_minus_z[:,1] + box[1]
    box_plus_y_minus_z[:,2] = box_plus_y_minus_z[:,2] - box[2]
    
    box_minus_y = coord.copy()
    box_minus_y[:,1] = box_minus_y[:,1] - box[1]
    
    box_minus_y_plus_z = coord.copy()
    box_minus_y_plus_z[:,1] = box_minus_y_plus_z[:,1] - box[1]
    box_minus_y_plus_z[:,2] = box_minus_y_plus_z[:,2] + box[2]
    
    box_minus_y_minus_z = coord.copy()
    box_minus_y_minus_z[:,1] = box_minus_y_minus_z[:,1] - box[1]
    box_minus_y_minus_z[:,2] = box_minus_y_minus_z[:,2] - box[2]   
    
    box_plus_z = coord.copy()
    box_plus_z[:,2] = box_plus_z[:,2] + box[2]
    
    box_minus_z = coord.copy()
    box_minus_z[:,2] = box_minus_z[:,2] - box[2]
    
    coord_with_boundary = np.concatenate((coord, box_plus_x,  box_plus_x_plus_y, box_plus_x_minus_y,  box_plus_x_plus_z, box_plus_x_minus_z, \
            box_plus_x_plus_y_plus_z, box_plus_x_plus_y_minus_z, box_plus_x_minus_y_plus_z, box_plus_x_minus_y_minus_z,  box_minus_x,  box_minus_x_plus_y, box_minus_x_minus_y, box_minus_x_plus_z, \
            box_minus_x_minus_z, box_minus_x_plus_y_plus_z, box_minus_x_plus_y_minus_z, box_minus_x_minus_y_plus_z,  box_minus_x_minus_y_minus_z,  box_plus_y, box_plus_y_plus_z, box_plus_y_minus_z,\
            box_minus_y, box_minus_y_plus_z, box_minus_y_minus_z, box_plus_z, box_minus_z ), axis=0)
    
    return coord_with_boundary
    
    
def impose_boundary_condition_27_box(name, path, file):

    '''
    Extracts the coordinates of 'name_atom' and imposes the periodic boundary conditions
    '''
    
    data = extract_data(name, path, file)
    coord = data['coord']
    box = data['box']
    
    coord_with_boundary = impose_boundary_condition_27_box_given_cood(coord, box)
    
    return coord_with_boundary
    
    
def get_adjacency_matrix_pbc(dist, coord, box):

    '''
    Returns the adjacency matrix with periodic boundary conditions
    '''
    
    coord = np.squeeze(np.array(coord))
    box = np.squeeze(np.array(box))
    
    N = np.shape(coord)[0]
    dim = 3
    
    ###############################################################################################################################################
    
    # this part of the code is from Yang, Y. & Fergus, M. (2020) Pair distances with PBC, url = https://yangyushi.github.io/science/2020/11/02/pbc_py.html
    dist_nd_sq = np.zeros(N * (N - 1) // 2)  # to match the result of pdist
    for d in range(dim):
        pos_1d = coord[:, d][:, np.newaxis]  # shape (N, 1)
        dist_1d = pdist(pos_1d)  # shape (N * (N - 1) // 2, )
        dist_1d[ np.where( dist_1d > box[d] * 0.5)] = dist_1d[ np.where( dist_1d > box[d] * 0.5)] - box[d]
        dist_nd_sq += dist_1d ** 2  # d^2 = dx^2 + dy^2 + dz^2
    dist_nd = np.sqrt(dist_nd_sq)
    
    ###############################################################################################################################################
    
    M = dist_nd <= dist       # boolean vector whose entry {i,j} is True if the distance between node v_i and node v_j is less or equal than dist
    M1 = M*1
    
    row = []
    col = []
    for i in range(N-1):
        for j in range(i+1, N):
            row.append(i)
            col.append(j)    

    D = csr_matrix((M1, (row, col)), shape=(N, N))
    A = D + D.T
    
    return csr_matrix(A)
    
def get_adjacency_matrix(dist, coord):

    '''
    Returns the adjacency matrix without periodic boundary conditions  
    '''

    coord = np.asmatrix(coord)    
    D = distance_matrix(coord, coord)
    D1 = D <= dist #boolean matrix whose entry {i,j} is True if the distance between node v_i and node v_j is less or equal than dist
    
    A = D1*1
    A = A - scipy.sparse.diags(A.diagonal())
    
    return csr_matrix(A)
    
def get_distance_matrix(dist, coord):

    '''
    Returns the adjacency matrix without periodic boundary conditions  
    '''

    coord = np.asmatrix(coord)    
    D = distance_matrix(coord, coord)
    D1 = D <= dist #boolean matrix whose entry {i,j} is True if the distance between node v_i and node v_j is less or equal than dist
    
    return csr_matrix(D1)

    


def angle_between(H, O1, O2):

    '''
    Returns the angle in radians between the points H O1 O2
    '''
    
    O1H =  H - O1
    O1O2 = O2 - O1
    V1 = O1H/np.linalg.norm(O1H)
    V2 = O1O2/np.linalg.norm(O1O2)
    res = np.clip(np.dot(V1, V2), -1.0, 1.0)
    angle = np.arccos(res)
    
    return angle
   
   
def oriented_adjacency_matrix(dist, coord_O, coord_H1, coord_H2, n_nodes_begin, boundary):

    '''
    Returns the sparse directed adjacency matrix
    '''

    D = get_distance_matrix(dist, coord_O)
    
    A_dist0 = D*1
    A_dist1 = A_dist0 - scipy.sparse.diags(A_dist0.diagonal())
    A_dist = csr_matrix(A_dist1)
    n, _ = np.shape(A_dist)
    
    lower = ((-30)*np.pi)/180
    upper = ((30)*np.pi)/180
    
    for i in range(n):
        vect = A_dist[i,:]
        _,J,_ = find(vect)
        oxyg_1 = np.squeeze(np.array(coord_O[i,:]))
        
        # I select the coordinates of H1 among all the periodic coordinates. I chose the one with the shortest distance from O.
        
        if int(boundary) == 1:
            hyd_1 = np.squeeze(np.array(coord_H1[i,:]))
            hyd_2 = np.squeeze(np.array(coord_H2[i,:]))
            
        elif int(boundary) == 2:
            i_box = int(i/n_nodes_begin)
            rest_box =27-i_box                        
            list1 = np.arange(0,rest_box)
            list2 = np.arange(-i_box,0)
            list_box = np.sort(list(list1)+list(list2))
            
            distance_vect_1 = []
            for k in list_box:
                hyd_1 = np.squeeze(np.array(coord_H1[i+k*n_nodes_begin,:]))
                distance_OH1 = np.linalg.norm(oxyg_1-hyd_1)
                distance_vect_1.append(distance_OH1)
                
            kk = distance_vect_1.index(min(distance_vect_1))
            kkk = list_box[kk]
            
            hyd_1 = np.squeeze(np.array(coord_H1[i+kkk*n_nodes_begin,:]))
            
            # I select the coordinates of H2 among all the periodic coordinates. I chose the one with the shortest distance from O.
            distance_vect_2 = []
            for k in list_box:
                hyd_2 = np.squeeze(np.array(coord_H2[i+k*n_nodes_begin,:]))
                distance_OH2 = np.linalg.norm(oxyg_1-hyd_2)
                distance_vect_2.append(distance_OH2)
                
            mm = distance_vect_2.index(min(distance_vect_2))
            mmm = list_box[mm]
            
            hyd_2 = np.squeeze(np.array(coord_H2[i+mmm*n_nodes_begin,:]))
            
        else:
            print('Value boundary not valid')
            sys.exit(1)
            
        for j in J:
            oxyg_2 = np.squeeze(np.array(coord_O[j,:]))
            deg_1 = angle_between(hyd_1, oxyg_1, oxyg_2)
            deg_2 = angle_between(hyd_2, oxyg_1, oxyg_2) 
            
            if lower <= deg_1 <= upper or lower <= deg_2 <= upper:
                pass
            else:
                A_dist[i,j] = 0
            #endif
    
    A_dist.eliminate_zeros()

    return A_dist
  
  
def parameters_graph_directed():

    '''
    Returns the parameters that allow you to build the directed graph: 
    1) name_atom is the name of the atom acting as a node, 
    2) dist is the threshold on the distance between two nodes,
    3) boundary indicates whether or not there are periodic boundary conditions.
    '''
    
    # We choose the distance
    fmt = 'Threshold to define the distance: '
    dist = input(fmt)
    dist = float(dist)
    
    while dist < 0.0:
        print('ERROR: the distance must be positive')
        dist = input(fmt)
        dist = float(dist)
        
    print('Two nodes are connected if their distance is less than {} and if the angle is between + - 30 degrees\n'.format(dist))
    
    # We coose the atom name for the nodes
    fmt = 'Atom name which defines the nodes. Press ENTER to default name (OW1) '
    
    name = input(fmt)
    if not name:
        name = 'OW1'
        
    # impose periodic boundary condition
    fmt = 'Do you want to impose periodic boundary conditions: (1) no, (2) yes, 27 boxes: '
    
    boundary = input(fmt)
    
    return {'dist':dist, 'name':name, 'boundary':boundary}
        
        
def parameters_graph_undirected():

    '''
    Returns the parameters that allow you to build the undirected graph: 
    1) name_atom is the name of the atom acting as a node, 
    2) dist is the threshold on the distance between two nodes,
    3) boundary indicates whether or not there are periodic boundary conditions.
    '''

    # We choose the distance
    fmt = 'Threshold to define the edge: '
    dist = input(fmt)
    dist = float(dist)
    
    while dist < 0.0:
        print('ERROR: the distance must be positive')
        dist = input(fmt)
        dist = float(dist)
        
    # We coose the atom name for the nodes
    fmt = 'Atom name which defines the nodes. Press ENTER to default name (OW1) '
    
    name = input(fmt)
    if not name:
        name = 'OW1'
        
    # impose periodic boundary condition
    fmt = 'Do you want to impose periodic boundary conditions: (1) no, (2) yes, 27 boxes: '
    
    boundary = input(fmt)
    
    return {'dist':dist, 'name':name, 'boundary':boundary}


def get_box(file, dist, name, boundary, n_nodes_begin):

    '''
    Extracts the values of the box
    '''
        
    path, filename = os.path.split(file)

    data = extract_data(name, path, file)
    box = data['box']
        
    return box


def get_graph(file, dist, name, boundary, n_nodes_begin):

    '''
    Construts the adjacency matrix and the graph (undirected case). Returns:
    - A = adjacency matrix, 
    - G = the graph, 
    - pos = phisical position of the nodes in the central box,
    - coord = coordinates of the nodes with PBC,
    - n_edges = number of edges,
    - box = values of the box.
    '''
        
    path, filename = os.path.split(file)
    data = extract_data(name, path, file)
    coord = data['coord']
    box = data['box']
    
    if int(boundary) == 1:
        A = get_adjacency_matrix(dist, coord)
    elif int(boundary) == 2:
        A = get_adjacency_matrix_pbc(dist, coord, box)
    else:
        print('Value not valid')
        sys.exit(1) 
    
    G = nx.from_scipy_sparse_matrix(A)
    
    n_edges = G.number_of_edges()
    
    pos={}
    for i in range(n_nodes_begin):
        position = np.squeeze(np.array(coord[i, 0:2]))
        pos[i] = position
        
    return A, G, pos, coord, n_edges, box
    
    
def get_graph_G(file, dist, name, boundary, n_nodes_begin):

    '''
    Construts the adjacency matrix, but  return only the graph and the number of edges
    '''
        
    path, filename = os.path.split(file)
    
    data = extract_data(name, path, file)
    coord = data['coord']
    box = data['box']
    
    if int(boundary) == 1:
        A = get_adjacency_matrix(dist, coord)
    elif int(boundary) == 2:
        A = get_adjacency_matrix_pbc(dist, coord, box)
    else:
        print('Value not valid')
        sys.exit(1)
    
    G = nx.from_scipy_sparse_matrix(A)
    
    n_edges = G.number_of_edges()
        
    return G, n_edges
    
        
        
def get_graph_A(file, dist, name, boundary):

    '''
    Returns only the adjacency matrix
    '''
        
    path, filename = os.path.split(file)
    data = extract_data(name, path, file)
    coord = data['coord']
    box = data['box']
    
    if int(boundary) == 1:
        A = get_adjacency_matrix(dist, coord)
    elif int(boundary) == 2:
        A = get_adjacency_matrix_pbc(dist, coord, box)
    else:
        print('Value not valid')
        sys.exit(1)
        
    return A
     
     
def get_graph_directed(file, dist, name, boundary, n_nodes_begin):

    '''
    Construts the adjacency matrix and the graph
    '''
        
    path, filename = os.path.split(file)

    data_O = extract_data(name, path, file)
    coord_O = data_O['coord']
    box = data_O['box']
    
    data_H1 = extract_data('HW2', path, file)
    coord_H1 = data_H1['coord']
    
    data_H2 = extract_data('HW3', path, file)
    coord_H2 = data_H2['coord']
    
    if int(boundary) == 1:
        coord2_O = np.asmatrix(coord_O)
        coord2_H1 = np.asmatrix(coord_H1)
        coord2_H2 = np.asmatrix(coord_H2)
    elif int(boundary) == 2:
        coord2_O =  impose_boundary_condition_27_box_given_cood(coord_O, box)
        coord2_H1 =  impose_boundary_condition_27_box_given_cood(coord_H1, box)
        coord2_H2 =  impose_boundary_condition_27_box_given_cood(coord_H2, box)
    else:
        print('Value boundary not valid')
        sys.exit(1)

    A = oriented_adjacency_matrix(dist, coord2_O, coord2_H1, coord2_H2, n_nodes_begin, boundary)
    G = nx.from_scipy_sparse_matrix(A, create_using=nx.DiGraph)
    
    n_edges = G.number_of_edges()
    
    pos={}
    for i in range(n_nodes_begin):
        position = np.squeeze(np.array(coord2_O[i, 0:2]))
        pos[i] = position
        
    
    return A, G, pos, coord2_O, n_edges, box


def compute_NTC_directed(file, n_nodes_begin, path2, dist, name, boundary, normalization, beta):

    '''
    Returns bc and rc for NTC
    '''
    
    _, G, _, _, n_edges, _ = get_graph_directed(file, dist, name, boundary,n_nodes_begin)      
    
    # TOTAL DIRECTED COMMUNICABILITY        
    bc, rc = total_directed_communicability(G, t=beta)
    
    if normalization == 1:      
        sub_bc = {}
        sub_rc = {}
        for i in range(n_nodes_begin):
            sub_bc[i] = bc[i]
            sub_rc[i] = rc[i]
            
    elif normalization == 2:
        sub_bc = {}
        sub_rc = {}
        for i in range(n_nodes_begin):
            sub_bc[i] = bc[i]/n_edges
            sub_rc[i] = rc[i]/n_edges
    
    else:
        print('Value not valid')
        sys.exit(1)
            
    with open(os.path.join(path2, 'directed_NTC_beta_'+ str(beta)+'.txt'), 'a+') as fobj:
        for ii in range(n_nodes_begin):
            fobj.write('{:=4d}      {:=18.6f}      {:=18.6f}\n'.format(ii, sub_bc[ii],sub_rc[ii]))
        fobj.write('\n')


def compute_NTC(file, beta, n_nodes_begin, path2, dist, name, boundary, normalization):

    '''
    Returns NTC and the execution time
    '''    
    
    if normalization == 1:     
        G, _ = get_graph_G(file, dist, name, boundary,n_nodes_begin)         
        
        # TOTAL COMMUNICABILITY    
        start = time.process_time()
        sub_tot_comm = total_communicability(G,  t = beta)  
        end = time.process_time()
        execution_time = end-start
        
            
    elif normalization == 2:
        G, n_edges = get_graph_G(file, dist, name, boundary,n_nodes_begin)  
        
        # TOTAL COMMUNICABILITY  
        start = time.process_time()
        tot_comm = total_communicability(G, t = beta)          
        sub_tot_comm = {}
        for i in range(n_nodes_begin):
            sub_tot_comm[i] = tot_comm[i]/n_edges
        end = time.process_time()
        execution_time = end-start
            
            
    elif normalization == 3:
        G, _ = get_graph_G(file, dist, name, boundary,n_nodes_begin)  
        
        # TOTAL COMMUNICABILITY 
        start = time.process_time()
        tot_comm = total_communicability(G, t = beta)  
        tnc = total_network_communicability(G)
        sub_tot_comm = {}
        for i in range(n_nodes_begin):
            sub_tot_comm[i] = tot_comm[i]/tnc
        end = time.process_time()
        execution_time = end-start
            
    elif normalization == 4:
        A, G, _, _, _, _ = get_graph(file, dist, name, boundary, n_nodes_begin)
        
        # TOTAL COMMUNICABILITY    
        start = time.process_time()
        tot_comm = total_communicability(G, t = beta) 
        
        A_sparse = A.asfptype()
        vals, vecs = scipy.sparse.linalg.eigs(A_sparse, k=6, which = 'LR')
        spectral_radius = np.real(vals[0])
        sub_tot_comm = {}
        for i in range(n_nodes_begin):
            sub_tot_comm[i] = tot_comm[i]/math.exp(spectral_radius * beta)
        with open(os.path.join(path2, 'spectral_radius.txt'), 'a+') as fobj:
            fobj.write('\n')
            fobj.write('{:=18.6f} \n'.format(spectral_radius))
        end = time.process_time()
        execution_time = end-start
    
    else:
        print('Value not valid')
        sys.exit(1)
        
    
    with open(os.path.join(path2, 'NTC_beta_'+ str(beta)+ '_boundary_' + str(boundary)+'_normalization_'+str(normalization) + '.txt'), 'a+') as fobj:
        for ii in range(n_nodes_begin):
            fobj.write('{:=4d}      {:=18.6f}\n'.format(ii, sub_tot_comm[ii]))
        fobj.write('\n')
        
    
    
    return execution_time, sub_tot_comm
   
   
def katz_direct_solver_sparse_matrix(G, alpha):

    '''
    Returns Katz centrality using a direct solver for sparse matrices
    '''

    n = G.number_of_nodes()
    nodelist = list(np.arange(n))
    A = nx.adjacency_matrix(G, nodelist = nodelist)
    
    b = np.ones((n, 1))    
    
    centrality = spsolve(scipy.sparse.identity(n) - (alpha * A), b)
    kz = dict(zip(nodelist, map(float, centrality)))
    
    return kz


def compute_degree(file, n_nodes_begin, path2, dist, name, boundary):

    '''
    Returns degree centrality
    '''

    G, _= get_graph_G(file, dist, name, boundary,n_nodes_begin)     
    n = G.number_of_nodes()
    
    # DEGREE      
    degree = nx.degree_centrality(G)
    
    sub_deg = {}
    for i in range(n_nodes_begin):
        sub_deg[i] = degree[i]*(n-1)
    
    with open(os.path.join(path2, 'DEGREE'+ '_boundary_' + str(boundary) + '.txt'), 'a+') as fobj:
        fobj.write('\n')
        for ii in range(n_nodes_begin):
            fobj.write('{:=4d}      {:=18.1f}\n'.format(ii, sub_deg[ii]))
            
def compute_eigenvector(file, n_nodes_begin, path2, dist, name, boundary):

    '''
    Returns eigenvector centrality
    '''
    
    G, _= get_graph_G(file, dist, name, boundary,n_nodes_begin)     

    # EIGENVECTOR CENTRALITY
    eig = nx.eigenvector_centrality_numpy(G)
    
    sub_eig = {}
    for i in range(n_nodes_begin):
        sub_eig[i] = eig[i]
    
    with open(os.path.join(path2, 'EIGENVECTOR'+ '_boundary_' + str(boundary) + '.txt'), 'a+') as fobj:
        fobj.write('\n')
        for ii in range(n_nodes_begin):
            fobj.write('{:=4d}      {:=18.6e}\n'.format(ii, sub_eig[ii]))
                
                
def compute_directed_degree(file, n_nodes_begin, path2, dist, name, boundary):

    '''
    Returns bc and rc for degree centrality
    '''
    _, G, _, _, _, _ = get_graph_directed(file, dist, name, boundary,n_nodes_begin)      
    n = G.number_of_nodes()
    
    # DEGREE      
    bc = nx.out_degree_centrality(G) #np.squeeze(np.array(A@np.ones((n,1))))
    rc = nx.in_degree_centrality(G) #np.squeeze(np.array(np.ones((1,n))@A))
    
    sub_deg_bc = {}
    for i in range(n_nodes_begin):
        sub_deg_bc[i] = bc[i]*(n-1)
        
    sub_deg_rc = {}
    for i in range(n_nodes_begin):
        sub_deg_rc[i] = rc[i]*(n-1)
    
    with open(os.path.join(path2, 'directed_DEGREE.txt'), 'a+') as fobj:
        for ii in range(n_nodes_begin):
            fobj.write('{:=4d}      {:=18.1f}      {:=18.1f}\n'.format(ii, sub_deg_bc[ii], sub_deg_rc[ii]))
        fobj.write('\n')



def compute_directed_NTC_degree(file, n_nodes_begin, path2, dist, name, boundary, normalization, beta):

    '''
    Returns bc and rc for NTC
    '''

    _, G, _, _, n_edges, _ = get_graph_directed(file, dist, name, boundary,n_nodes_begin)      
    
    # TOTAL DIRECTED COMMUNICABILITY        
    bc, rc = total_directed_communicability(G, t=beta)
    
    if normalization == 1:      
        sub_bc = {}
        sub_rc = {}
        for i in range(n_nodes_begin):
            sub_bc[i] = bc[i]
            sub_rc[i] = rc[i]
            
    elif normalization == 2:
        sub_bc = {}
        sub_rc = {}
        for i in range(n_nodes_begin):
            sub_bc[i] = bc[i]/n_edges
            sub_rc[i] = rc[i]/n_edges
    
    else:
        print('Value not valid')
        sys.exit(1)
    
    with open(os.path.join(path2, 'directed_NTC_beta_'+ str(beta)+'.txt'), 'a+') as fobj:
        for ii in range(n_nodes_begin):
            fobj.write('{:=4d}      {:=18.6f}      {:=18.6f}\n'.format(ii, sub_bc[ii],sub_rc[ii]))
        fobj.write('\n')
        
    n = G.number_of_nodes()
    
    # DEGREE      
    bc = nx.out_degree_centrality(G) #np.squeeze(np.array(A@np.ones((n,1))))
    rc = nx.in_degree_centrality(G) #np.squeeze(np.array(np.ones((1,n))@A))
    
    sub_deg_bc = {}
    for i in range(n_nodes_begin):
        sub_deg_bc[i] = bc[i]*(n-1)
        
    sub_deg_rc = {}
    for i in range(n_nodes_begin):
        sub_deg_rc[i] = rc[i]*(n-1)
    
    with open(os.path.join(path2, 'directed_DEGREE.txt'), 'a+') as fobj:
        for ii in range(n_nodes_begin):
            fobj.write('{:=4d}      {:=18.1f}      {:=18.1f}\n'.format(ii, sub_deg_bc[ii], sub_deg_rc[ii]))
        fobj.write('\n')


def compute_subgraph(file, beta, n_nodes_begin, path2, dist, name, boundary):

    '''
    Compute subgraph centrality
    '''

    G, _ = get_graph_G(file, dist, name, boundary, n_nodes_begin)     
    
    # SUBGRAPH     
    start = time.process_time()
    subgraph = subgraph_centrality_Bucci(G, t=beta)
    end = time.process_time()
       
    
    sub_subgraph = {}
    for i in range(n_nodes_begin):
        sub_subgraph[i] = subgraph[i]
    
    with open(os.path.join(path2, 'SUBGRAPH_beta_'+ str(beta) + '_boundary_' + str(boundary) +'.txt'), 'a+') as fobj:
        fobj.write('\n')
        for ii in range(n_nodes_begin):
            fobj.write('{:=4d}      {:=18.6f}\n'.format(ii, sub_subgraph[ii]))
            
    execution_time = end-start
    
    return execution_time, subgraph
    
    
def compute_CLOSENESS(file, n_nodes_begin, path2, dist, name, boundary):

    '''
    Returns the execution time of the closeness centrality
    '''

    G, _ = get_graph_G(file, dist, name, boundary, n_nodes_begin)    
    start = time.process_time()
    cl = nx.closeness_centrality(G, u=None, distance=None,wf_improved=True)
    end = time.process_time()
    
    with open(os.path.join(path2, 'CLOSENESS_boundary_' + str(boundary) +'.txt'), 'a+') as fobj:
        fobj.write('\n')
        for ii in range(n_nodes_begin):
            fobj.write('{:=4d}      {:=18.6f}\n'.format(ii, cl[ii]))
            
    execution_time = end-start
    
    return execution_time
    
    
def compute_BETWEENNESS(file, n_nodes_begin, path2, dist, name, boundary):

    '''
    Returns the execution time of the betweenness centrality 
    '''

    G, _ = get_graph_G(file, dist, name, boundary, n_nodes_begin)    
    start = time.process_time()
    bw = nx.betweenness_centrality(G, normalized=False)
    end = time.process_time()
    
            
        
    with open(os.path.join(path2, 'BETWEENNESS_boundary_' + str(boundary) +'.txt'), 'a+') as fobj:
        fobj.write('\n')
        for ii in range(n_nodes_begin):
            fobj.write('{:=4d}      {:=18.6f}\n'.format(ii, bw[ii]))
            
    execution_time = end-start
    
    return execution_time
  
  
def compute_KATZ(file, n_nodes_begin, path2, dist, name, boundary, alpha):

    '''
    Returns the execution time of the Katz centrality using a direct solver for sparse matrices
    '''

    # direct method with dense matrix
    #G, _ = get_graph_G(file, dist, name, boundary, n_nodes_begin) 
    #kz = nx.katz_centrality_numpy(G, alpha = alpha, beta = 1.0, normalized=False) 
    
    # powers method
    #G, _ = get_graph_G(file, dist, name, boundary, n_nodes_begin) 
    #kz = nx.katz_centrality(G, alpha = alpha, beta = 1.0, normalized=False) 
    
    # direct method with sparse matrix
    G, _ = get_graph_G(file, dist, name, boundary, n_nodes_begin) 
    
    start = time.process_time()
    
    kz = katz_direct_solver_sparse_matrix(G, alpha)
    
    end = time.process_time()

    with open(os.path.join(path2, 'KATZ_'+str(alpha)+ '_boundary_' + str(boundary) +'.txt'), 'a+') as fobj:
        fobj.write('\n')
        for ii in range(n_nodes_begin):
            fobj.write('{:=4d}      {:=18.6f}\n'.format(ii, kz[ii]))
            
    execution_time = end-start
    
    return execution_time
    

def plot_network_3d(coord, dist, path, name_img, color = None, vmin = None, vmax = None): 

    '''
    Plot 3d of a frame
    '''

    n, _ = np.shape(coord)
    path2 = os.path.split(path)[0]
    
    if color == None:
        color = np.ones(n)
  
    if vmin == None: 
        vmin = min(color)
        vmax = max(color)

    fig = plt.figure(0, figsize=(10,5))
    ax1 = fig.add_subplot(1, 1, 1, projection='3d')
    plot1 = ax1.scatter(coord[:,0], coord[:,1], coord[:,2], c=color, cmap=plt.cm.jet, vmin=vmin, vmax=vmax, s=30  )
            
    for i in range(n):
        for j in range(i+1, n):
            matrix_i_j = np.linalg.norm(coord[i]-coord[j])
            if matrix_i_j <= dist:
                ax1.plot([coord[i,0], coord[j,0]], [coord[i,1], coord[j,1]],[coord[i,2], coord[j,2]], 'k-', linewidth=0.3)
                
    
    if any(color != np.ones(n)):
        fig.colorbar(plot1, ax = ax1)
    ax1.set_rasterized(True)
    fig.savefig(os.path.join(path2, 'fig_' + str(name_img) + '.png'))
    fig.savefig(os.path.join(path2, 'fig_' + str(name_img) + '.eps'))
    plt.show()
    plt.close()
    

def plot_network_2d(pos, G, range_nodes, path, name_img, color = None, vmin = None, vmax = None): 

    '''
    Plot 2d of a frame
    '''

    path2 = os.path.split(path)[0]
    
    n = len(range_nodes)
    
    if color == None:
        color = np.ones(n)
  
    if vmin == None: 
        vmin = min(color)
        vmax = max(color)

    fig, axes = plt.subplots(nrows=1, ncols=1,  figsize=(10, 5)) 
    
    ec = nx.draw_networkx_edges(G, pos, alpha=0.2, ax = axes)
    nc = nx.draw_networkx_nodes(G, pos, nodelist=range_nodes, node_color=color, cmap=plt.cm.jet, vmin=vmin, vmax=vmax, ax = axes)
    lc = nx.draw_networkx_labels(G, pos, font_size=6, ax = axes)

    if any(color != np.ones(n)):
        fig.colorbar(nc, ax = axes)
    
    axes.set_rasterized(True)
    fig.savefig(os.path.join(path2, 'fig_' + str(name_img) + '.png'))
    fig.savefig(os.path.join(path2, 'fig_' + str(name_img) + '.eps'))
    plt.show()
    plt.close()
    

def count_cycles(file, n_nodes_begin, path2, dist, name, boundary):

    '''
    Counts number of edges and cycles of length 3,4 and 5
    '''

    A = get_graph_A(file, dist, name, boundary)   
    G = nx.from_scipy_sparse_matrix(A)
    
    n_edges = G.number_of_edges()    
    n,_  = np.shape(A)
    
    with open(os.path.join(path2, 'number_edges.txt'), 'a+') as fobj:
        fobj.write('Number edges = {:d} \n'.format(n_edges))
        
    A2 = A@A
    
    n_nodes, _ = np.shape(A)
    range_nodes = np.arange(0, n_nodes)
    degree = nx.degree_centrality(G)
    degree = [degree.get(node)*(n_nodes-1) for node in range_nodes]
    f1 = [degree[i]*(degree[i]-1) for i in range(n_nodes)]
    F1 = sum(f1)/2   
    
    A3 = A2@A
    a3 = A3.diagonal().sum()
    F2 = a3/6
    
    A4 = A3@A
    m = n_edges
    a4 = A4.diagonal().sum()
    F5 = (a4 - 4*F1 - 2*m)/8
    
    f6 = [A3[i,i]*(degree[i]-2)/2 for i in range(n_nodes) if degree[i] > 2]
    F6 = sum(f6)
    
    A5 = A4@A
    a5 = A5.diagonal().sum()
    F8 = (1/10)*(a5 -30*F2 -10*F6)
    
    with open(os.path.join(path2, 'cycles_3.txt'), 'a+') as fobj:
        fobj.write('Number cycles 3 = {:f} \n'.format( F2))
        
    with open(os.path.join(path2, 'cycles_4.txt'), 'a+') as fobj:
        fobj.write('Number cycles 4 = {:f} \n'.format(F5))
        
    with open(os.path.join(path2, 'cycles_5.txt'), 'a+') as fobj:
        fobj.write('Number cycles 5 = {:f} \n'.format(F8))
        
    return n_edges, F2, F5, F8
            
    

def estrada_index(file, n_nodes_begin, path2, dist, name, boundary):

    '''
    Returns the Estrada index
    '''

    G, _= get_graph_G(file, dist, name, boundary, n_nodes_begin)   
    estrada_index_values = nx.estrada_index(G)
    
    with open(os.path.join(path2, 'estrada_index.txt'), 'a+') as fobj:
        fobj.write('Estrada index = {:f} \n'.format(estrada_index_values))
    
    return estrada_index_values
    
    
def watts_strogatz(file, n_nodes_begin, path2, dist, name, boundary):

    '''
    Returns the Wattz-Strogatz coefficient clustering
    '''

    G, _= get_graph_G(file, dist, name, boundary, n_nodes_begin)   
    watts_strogatz_values = nx.average_clustering(G)
    
    with open(os.path.join(path2, 'watts_strogatz.txt'), 'a+') as fobj:
        fobj.write('Watts-Strogatz coeff. clustering = {:f} \n'.format(watts_strogatz_values))
        
    return watts_strogatz_values
  
  
def transitivity_index(file, n_nodes_begin, path2, dist, name, boundary):

    '''
    Returns the trasitivity index
    '''

    G, _= get_graph_G(file, dist, name, boundary, n_nodes_begin) 

    number_of_triangles = sum(nx.triangles(G).values())/3
    
    range_nodes = np.arange(0, n_nodes_begin)
    
    degree_norm = nx.degree_centrality(G)
    degree = [degree_norm.get(node)*(n_nodes_begin-1) for node in range_nodes]
    P2 = sum([x*(x-1)/2 for x in degree])
    
    transitivity_index_values = 3*number_of_triangles/P2
    
    with open(os.path.join(path2, 'transitivity_index.txt'), 'a+') as fobj:
        fobj.write('Transitivity index = {:f} \n'.format(transitivity_index_values))
        
    return transitivity_index_values
            

def assortativity(file, n_nodes_begin, path2, dist, name, boundary): 

    '''
    Returns the assortativity coefficient of the graph G
    '''

    A, G, _, _, n_edges, _ = get_graph(file, dist, name, boundary, n_nodes_begin) 
    
    assortativity_values = nx.degree_assortativity_coefficient(G)

    A2 = A@A
    
    n_nodes, _ = np.shape(A)
    range_nodes = np.arange(0, n_nodes)
    degree = nx.degree_centrality(G)
    degree = [degree.get(node)*(n_nodes-1) for node in range_nodes]
    degree_minusone = np.asmatrix([i-1 for i in degree])
    f1 = [degree[i]*(degree[i]-1) for i in range(n_nodes)]
    F1 = sum(f1)/2   
    
    A3 = A2@A
    a3 = A3.diagonal().sum()
    F2 = a3/6
    
    m = n_edges    
    
    path2_path1 = F1/m       
    
    S2 = (degree_minusone@A@(degree_minusone.T))/2 - 3*F2
    
    path3_path2 = S2/F1
    
    with open(os.path.join(path2, 'assortativity.txt'), 'a+') as fobj:
        fobj.write('Assortativity = {:f} \n'.format(assortativity_values))
        fobj.write('P2/P1 = {:f} \n'.format(path2_path1))
        fobj.write('P3/P2 = {:f} \n'.format(path3_path2[0,0]))
        fobj.write('\n')
        
    return assortativity_values, path2_path1, path3_path2[0,0]
    
def entropy_subgraph(file, n_nodes_begin, path2, dist, name, boundary):

    '''
    Returns the entropy of a graph based on the subgraph centrality: S(G, \beta) = - \sum_{i=1}^N p_i *ln(p_i), where p_i = \frac{exp(\beta A)_{i i} }{Tr(exp(\beta A))}
    '''

    G, _= get_graph_G(file, dist, name, boundary, n_nodes_begin) 
    largest_cc = max(nx.connected_components(G), key=len)
    G1 = G.subgraph(largest_cc)    
    node_list = list(G1.nodes())
    
    sub_values = subgraph_centrality_Bucci(G1, t=1.0)

    EE = sum(sub_values.values())
    
    p_vect = [sub_values[ii]/EE for ii in node_list]
    
    H = [p_vect[ii]*log(p_vect[ii]) for ii in range(n) ]
    
    entropy_values = -sum(H)
    
    with open(os.path.join(path2, 'entropy_subgraph.txt'), 'a+') as fobj:
        fobj.write('Entropy = {:f} \n'.format(entropy_values))

    return entropy_values
    
    
def entropy_TC(file, n_nodes_begin, path2, dist, name, boundary):

    '''
    Returns the entropy of a graph based on the total communicability : TC(G, \beta) = - \sum_{i=1}^N p_i *ln(p_i), where p_i = \frac{(exp(\beta A) 1)_{i} }{1^T (exp(\beta A)) 1}
    '''

    G, _= get_graph_G(file, dist, name, boundary, n_nodes_begin) 
    largest_cc = max(nx.connected_components(G), key=len)
    G1 = G.subgraph(largest_cc)    
    node_list = list(G1.nodes())
    n = len(node_list)
    
    TC_values = total_communicability(G1, t = 1.0)  

    NTC = sum(TC_values.values())
    
    p_vect = [TC_values[ii]/NTC for ii in node_list]
    
    H = [p_vect[ii]*log(p_vect[ii]) for ii in range(n) ]
    
    entropy_values = -sum(H)
    
    with open(os.path.join(path2, 'entropy_TC.txt'), 'a+') as fobj:
        fobj.write('Entropy = {:f} \n'.format(entropy_values))

    return entropy_values
    
    
def entropy_von_neumann(file, n_nodes_begin, path2, dist, name, boundary):

    '''
    Returns the Von-Neumann entropy: E(G) = - Tr( p * ln(p)), where p is the density matrix, p = L/Tr(L), and L is the Laplacian graph
    '''
    
    G, _= get_graph_G(file, dist, name, boundary, n_nodes_begin) 
    largest_cc = max(nx.connected_components(G), key=len)
    G1 = G.subgraph(largest_cc)    
    A = nx.adjacency_matrix(G1)
     
    n = A.shape[0]
    e1 = np.asmatrix(np.ones((n,1)))
    degree = A@e1
    deg = np.squeeze(np.array(degree))
    
    L = scipy.sparse.diags(deg) - A
    
    trace_L = sum(deg)
    
    Den = L/trace_L
    
    w,_ = scipy.sparse.linalg.eigsh(Den, k=n-1, which='LM')
    
    H = [lambda_i * log(lambda_i) for lambda_i in w ]
    
    entropy_values = -sum(H)
    
    with open(os.path.join(path2, 'entropy_Von_Neumann.txt'), 'a+') as fobj:
        fobj.write('Entropy = {:f} \n'.format(entropy_values))

    return entropy_values
    
    
def entropy_CM_for_files(n_nodes_begin, path2, boundary, list_path_files, CM):

    '''
    Returns the entropy value given the subgraph values or the total communicability values
    
    TC(G, \beta) = - \sum_{i=1}^N p_i *ln(p_i), where p_i = \frac{(exp(\beta A) 1)_{i} }{1^T (exp(\beta A)) 1}
    
    S(G, \beta) = - \sum_{i=1}^N p_i *ln(p_i), where p_i = \frac{exp(\beta A)_{i i} }{Tr(exp(\beta A))}
    '''
    
    entropy_vect = []
    
    print('Select the file that contains the centrality measure values \n')
    
    root = Tk()
    root.withdraw()
    file_txt_centrality = filedialog.askopenfilename()
    root.destroy()
    
    entropy_vect = []
    
    with open(file_txt_centrality, 'r') as data:
        data.seek(0)
        for ff in range(2):
            data.readline()
            
        for frame in range(len(list_path_files)):
            sub = []
            for i in range(n_nodes_begin):
                cent = data.readline()
                cent = float(cent[9:30])
                sub.append(cent)
            data.readline()
            
            EE = sum(sub)
            p_vect = [sub[ii]/EE for ii in range(n_nodes_begin)]
    
            H = [p_vect[ii]*log(p_vect[ii]) for ii in range(n_nodes_begin) ]
            
            entropy_values = -sum(H)
            
            with open(os.path.join(path2, 'entropy_'+str(CM) +'_from_file.txt'), 'a+') as fobj:
                fobj.write('Entropy = {:f} \n'.format(entropy_values))
                
            entropy_vect.append(entropy_values)
            
    print('Average entropy = {} \n'.format(np.mean(entropy_vect)))
    
    return
            
def max_min_eigenvalues(n_nodes_begin, path2, boundary): 

    '''
    Returns the max and min eigenvalue
    '''
    
    print('Select the folder that contains the eigenvalues \n')
    
    root = Tk()
    root.withdraw()
    path = filedialog.askdirectory()
    root.destroy()
    
    n = n_nodes_begin
        
    min_eig = []
    max_eig = []
    diff = []
    sorted_list_files = [os.path.join(path, str(file)) for file in os.listdir(path)]
    sorted_list_files.sort(key=os.path.getctime) #sorted by creation time

    for file in sorted_list_files:    # the order in which the files are opened does not matter to me
        eig = []
        with open(file, 'r') as data:
            data.seek(0)
            for i in range(n):
                d = data.readline()
                eig.append(float(d))
        
        mi_e = min(eig)
        ma_e = max(eig)
        
        min_eig.append(mi_e)
        max_eig.append(ma_e)
        diff.append(ma_e + mi_e)
        
    print('Mean value of the maximum eigenvalue = {} \n'.format(np.mean(max_eig)))
    print('Mean value of the minimum eigenvalue = {} \n'.format(np.mean(min_eig)))
    print('Mean value of the differnce = {} \n'.format(np.mean(diff)))
    print('Maximum eigenvalue = {} \n'.format(max(max_eig)))
    
    return 
    
    
    
    
def bipartivity_for_files(n_nodes_begin, path2, boundary): 

    '''
    Returns the bipartivity value given the eigenvalues. The eigenvalues are in a folder. In the folder, there is a file.txt, with the eigenvalues, for each frame of the trajectory
    '''
    
    bipartivity_vect = []
    
    print('Select the folder that contains the eigenvalues \n')
    
    root = Tk()
    root.withdraw()
    path = filedialog.askdirectory()
    root.destroy()
    
    n = n_nodes_begin
    sorted_list_files = [os.path.join(path, str(file)) for file in os.listdir(path)]
    sorted_list_files.sort(key=os.path.getctime) #sorted by creation time

    for file in sorted_list_files:    
        eig = []
        
        with open(file, 'r') as data:
            data.seek(0)
            for i in range(n):
                d = data.readline()
                eig.append(float(d))
        eig.sort(reverse=True)
        num = sum(np.cosh(eig))
        den = sum(np.exp(eig))
        bipartivity_values = num/den
        
        with open(os.path.join(path2, 'bipartivity.txt'), 'a+') as fobj:
            fobj.write('Bipartivity = {:f} \n'.format(bipartivity_values))
            
        bipartivity_vect.append(bipartivity_values)
                        
    print('Average number of bipartivity value = {} \n'.format(np.mean(bipartivity_vect)))
        
    return 
    
def bipartivity(file, n_nodes_begin, path2, path3, dist, name, boundary,ii): 

    '''
    Returns the bipartivity value of the graph G. Firstly it computes the eigenvalues and it saves them in a folder, then it computes the bipartivity measure
    '''
   
    eig = eigenvalues(file, n_nodes_begin, path3, dist, name, boundary, ii)
    eig.sort(reverse=True)
    num = sum(np.cosh(eig))
    den = sum(np.exp(eig))
    bipartivity_values = num/den    
    
    with open(os.path.join(path2, 'bipartivity.txt'), 'a+') as fobj:
        fobj.write('Bipartivity = {:f} \n'.format(bipartivity_values))
        
    return bipartivity_values
            
def sparsity_A(file, n_nodes_begin, path2, dist, name, boundary, frame):

    '''
    Sparsity plot of th adjacency matrix 
    '''
 
    A = get_graph_A(file, dist, name, boundary) 
    n = np.shape(A)[0]
    
    if n == n_nodes_begin*27:
    
        fig = plt.figure(0)
        plt.spy(A, markersize=0.3, color = 'k')
        plt.title('27 boxes')
        fig.savefig(os.path.join(path2, 'sparsity A_frame' + str(frame) + '_27_boxes.png'))
        plt.show()
        plt.close()
        
        fig2 = plt.figure(0)
        plt.spy(A[0:n_nodes_begin, 0:n_nodes_begin], markersize=0.3, color = 'k')
        plt.title('central box')
        fig2.savefig(os.path.join(path2, 'sparsity A_frame' + str(frame) + '_central_box.png'))
        plt.show()
        plt.close()
    else:
    
        fig2 = plt.figure(0)
        plt.spy(A[0:n_nodes_begin, 0:n_nodes_begin], markersize=0.3, color = 'k')
        plt.title('central box')
        fig2.savefig(os.path.join(path2, 'sparsity A_frame' + str(frame) + '_central_box.png'))
        plt.show()
        plt.close()
   

    
def density(file, n_nodes_begin, path2, dist, name, boundary):

    '''
    Returns the density of the central box, defining as n_nodes/volume
    '''

    box = get_box(file, dist, name, boundary,n_nodes_begin) 

    volume = box[0]*box[1]*box[2]
    dens = n_nodes_begin/volume
    
    with open(os.path.join(path2, 'density.txt'), 'a+') as fobj:
        fobj.write('density = {:f} \n'.format(dens))
    
    return dens
    
    
def density_graph(file, n_nodes_begin, path2, dist, name, boundary):

    '''
    Returns the density of the graph G
    '''

    G, _ = get_graph_G(file, dist, name, boundary,n_nodes_begin) 

    dens = nx.density(G)
    
    with open(os.path.join(path2, 'graph_density.txt'), 'a+') as fobj:
        fobj.write('graph density = {:f} \n'.format(dens))
    
    return dens


def eigenvalues(file, n_nodes_begin, path3, dist, name, boundary, ii):

    '''
    Computes the eigenvalues of the adjacency matrix. It can take a long time
    '''

    A = get_graph_A(file, dist, name, boundary) 
    n = np.shape(A)[0]
    
    w = scipy.linalg.eigh(A.todense(), eigvals_only=True)
     
    
    with open(os.path.join(path3, 'eigenvalues_'+str(ii)+'.txt'), 'a+') as fobj:
        for ii in range(n):
            fobj.write('{:f} \n'.format(w[ii]))
        fobj.write('\n')
        
    return list(w)
    
    
def save_matrix_matlab_format(file, n_nodes_begin, new_path, dist, name, boundary, kk):

    '''
    Saves a matrix in matlab format A.mat
    '''

    A = get_graph_A(file, dist, name, boundary) 
    
    scipy.io.savemat(os.path.join(new_path,'A_' + str(kk) +'.mat'), {'A':A})   
    
    return
    
def energy(file, n_nodes_begin, path2, dist, name, boundary):

    '''
    Returns the energy of the graph E(G) = \sum_{i=1}^N |\lambda_i|
    '''

    print('Select the folder that contains the eigenvalues \n')
    
    root = Tk()
    root.withdraw()
    path = filedialog.askdirectory()
    root.destroy()
    
    n = n_nodes_begin
        
    energy_vect = []
    sorted_list_files = [os.path.join(path, str(file)) for file in os.listdir(path)]
    sorted_list_files.sort(key=os.path.getctime) #sorted by creation time
    
    for file in sorted_list_files:    # the order in which the files are opened does not matter to me
        eig = []
        with open(file, 'r') as data:
            data.seek(0)
            for i in range(n):
                d = data.readline()
                eig.append(float(d))
                
        energy = np.sum(np.abs(eig))
        
        with open(os.path.join(path2, 'energy_graph.txt'), 'a+') as fobj:
            fobj.write('{:f} \n'.format(energy))
            fobj.write('\n')
        
        energy_vect.append(energy)
        
        
    print('Mean value of the energy = {} \n'.format(np.mean(energy_vect)))
    
    return
    
def ASPL_diameter_isolated(file, n_nodes_begin, path2, dist, name, boundary):

    '''
    Returns the average shortest path length, the diameter of the graph and the number of isolated points. We consider the largest connected component
    '''
    
    G, _ = get_graph_G(file, dist, name, boundary, n_nodes_begin) 
    largest_cc = max(nx.connected_components(G), key=len)
    G_s = G.subgraph(largest_cc)
    ASPL = nx.average_shortest_path_length(G)
    diam = nx.diameter(G)
    n_G_s = G_s.number_of_nodes()
    n_isolated = n_nodes_begin - n_G_s
    
    
    return ASPL, diam, n_isolated
    
def algebraic_connectivity(file, n_nodes_begin, path2, dist, name, boundary):

    '''
    Returns the algebraic connectivity of the graph
    '''
    
    G, _ = get_graph_G(file, dist, name, boundary, n_nodes_begin) 
    AG = nx.algebraic_connectivity(G)
    
    
    return AG
