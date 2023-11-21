# WaterNetworks
This is the second version of the code.  The code is implemented in Python 3.7 and utilizes the MDAnalysis library [2,3] to extract information from the molecular dynamics simulations  to construct the graphs.    Additionally, it  heavily relies on the NetworkX 2.6.3 software package  [1] for the creation, manipulation, and analysis of graphs, the Python packages NetworkSNS [4] and DyNetX [5], and open-source packages including NumPy [6], SciPy [7] and Matplotlib [8]. 

In this section,  the instructions for installing the code will be explained, as well as the capabilities it offers for calculations. 


## INSTALLATION: 
In order to use the code provided in the Github repository, 
 (https://github.com/ChiaraFaccio/WaterNetworks), it is necessary to have certain packages and libraries installed on your computer, namely NumPy, SciPy, Matplotlib, MDAnalysis, and NetworkX. It is necessary that these dependencies are installed before proceeding. Once installed, or their presence has been verified, one downloads the Python scripts **md_functions_water_network.py**, **md_link_to_main.py** and **md_main_water.py** from the Github repository and saves them in the same folder. The first and the second script contain all the functions that are called in the main script. To run the code, the user is requires to open the command prompt and navigate to the folder where the scripts are located. The main script is executed by entering \texttt{python md\_main\_water.py} or \texttt{python3 md\_main\_water.py} in the terminal, depending on the Python version used. 

 ## USER GUIDE
 The code is designed with a user-friendly interface, ensuring ease of use. It interacts with the user in different ways: for example by presenting questions to which the user responds by typing the answer, or by selecting it from a list of options. In other situations, the user is prompted to select files or folders. In this case, the code automatically opens a dialog box, allowing the user to choose the desired file or folder conveniently. This approach improves usability by making user interaction with the code simple and intuitive.

The code combines two Python libraries: MDAnalysis for the part concerning the manipulation of the molecular dynamics trajectories, and NetworkX for creating the graphs and calculating most of the centrality measures and global metrics.

The process begins by uploading the data through MDAnalysis. The code requests the user to upload the topology file, which contains a list of atoms, residues, their connectivity, etc. Accepted formats for the topology file are PSF, PDB, CRD, and GRO. The user has the possibility to include a trajectory file, which consists of a sequence of coordinates corresponding to the topology. If no trajectory file is provided, only the structure described in the topology file will be used. In contrast, if multiple trajectory files are provided, they will be concatenated. Single frame can be in PDB, CRD, or GRO formats, while time series data can be in DCD, XTC, TRR, or XYZ formats.

At this point, the user has the option to select whether he wants to construct an undirected or directed graph, and he defines the threshold to determine when two nodes are connected, measured in nanometers (nm). Additionally, the user can choose whether to apply periodic boundary conditions during the graph construction. In this version of the code,  for undirected networks, the imposition of the periodic boundary conditions in the graphs works also in the case of a triclinic box. However, for directed graphs, the implementation is limited to boxes with orthogonal faces.

The next step involves choosing the atoms to be used as nodes. The code allows the user to select atoms from different molecule types. For example, given a trajectory involving water, TFA, and HYD molecules (according to GROMACS notation), the user can consider the oxygen atoms OW1 from water and the atoms O05 and O06 from TFA as graph vertices. Furthermore, the user can decide  to save only the values of CMs for some atoms/nodes. In the aforementioned example, for instance, it is possible to select and save only the CMs values of the oxygen atoms in water.
 

Once all these choices are made, the code constructs a graph for each frame of the trajectory, enabling analysis using the tools from graph theory. The following schematic representation demonstrates the  computations that can be performed:

- **UNDIRECTED GRAPHS**
    - Centrality measure computation
       - total communicability
       - degree
       - subgraph centrality
       - closeness centrality
       - betweenness centrality
       - Katz centrality
       - eigenvector centrality
       - non-backtracking total communicability
   - Computation of other metrics
       - to count edges and cycles of length 3, 4 and 5  
       - to compute Estrada index  
       - to compute Watts-Strogatz clustering coefficient  
       - to compute transitivity index  
       - to compute bipartivity measures  
       - to compute the degree assortativity coefficient  
       - to compute eigenvalues  
       - to plot the sparsity pattern of a frame  
       - to compute entropy of the graph with subgraph  
       - to compute entropy of the graph with TC  
       - to compute entropy of the graph with graph Laplacian (von Neumann entropy)  
       - to compute the mean value of min/max eigenvalues
       - to compute the mean values of the density of the graphs  
       - to compute the mean values of the density of the boxes (number of molecules/volume box)  
       - to compute the energy of the graph 
       - to compute ASPL and diameter  
       - to compute the algebraic connectivity 
   - to plot a single graph in 3D 
   - to plot a single graph in 2D
   - to save the adjacency matrices in MATLAB format
   - Dynamic graph metrics
       - to compute dynamic communicability
       - to compute aggregated degree
   - Weighted graph metrics
       - to compute degree e TC 
       - to compute Estrada index
       - to compute the entropy of the graph with TC 
       - to compute the entropy of the graph with graph Laplacian (von Neumann entropy)
- **DIRECTED GRAPHS**
   - Centrality measures computation
       - in- and out-degree 
       - left and right eigenvector centrality of $A$ 
       - HITS 
       - gTC 

In the case of undirected graphs, the code allows the computation of various centrality measures and global metrics, and allows one to plot the graphs in 2D or 3D, considering a sequence of graphs as dynamic graphs, to assign a positive weight to each edge of the network. The user can also choose to save the adjacency matrices in MATLAB format if they want to perform additional analyses with MATLAB. For undirected graphs, only the computation of some centrality measures for digraphs is implemented. However, it is possible to modify the code, as it is released under the BSD-3-Clause license. Making changes is straightforward: it is sufficient to add a calculation option in the main file **python md_main_water.py** and the related functions in the other two scripts, using the previous functions as templates.

During the execution of the code, each result will be written in a .txt file, named with an intuitive name. These files will be saved in the same folder containing the topology/trajectory files; for this reason, we recommend having each trajectory in its own folder.


## REFERENCES:  
[1] Hagberg, A., Chult, D. S. & Swart, P. (2008) Exploring Network Structure, Dynamics, and Function using NetworkX. In Varoquaux, G., Vaught, T. & Millman, J., editors, Proceedings of the 7th Python in Science Conference, pages 11 – 15, Pasadena, CA USA, Networkx

[2] R. J. Gowers, M. Linke, J. Barnoud, T. J. Reddy, M. N. Melo, S. L. Seyler, J. Domanski,
D. L. Dotson, S. Buchoux, I. M. Kenney, et al. Mdanalysis: a python package for the
rapid analysis of molecular dynamics simulations. In Proceedings of the 15th python in
science conference, volume 98, page 105. SciPy Austin, TX, 2016.

[3] N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein. Mdanalysis: a
toolkit for the analysis of molecular dynamics simulations. Journal of computational
chemistry, 32(10):2319–2327, 2011.

[4] A. Bucci. NetworkSNS (v3.0.0). Zenodo. https://doi.org/10.5281/zenodo.8073704,
June 2023.

[5] Rossetti, G. (2020) DyNetx: dynamic network analysis library, v0.2.1, Zenodo, doi: 10.5281/zenodo.3953119

[6] Harris, C. R., Millman, K. J., van der Walt, S. J., Gommers, R., Virtanen, P., Cournapeau, D., Wieser, E., Taylor, J., Berg, S., Smith, N. J., Kern, R., Picus, M., Hoyer, S., van Kerkwijk, M. H., Brett, M., Haldane, A., del R´ıo, J. F., Wiebe, M., Peterson, P., Gerard-Marchant, P., Sheppard, K., Reddy, T., Weckesser, W., Abbasi, ´ H., Gohlke, C. & Oliphant, T. E. (2020) Array programming with NumPy. Nature, 585(7825), 357–362 NumPy

[7] Virtanen, P., Gommers, R., Oliphant, T. E., Haberland, M., Reddy, T., Cournapeau, D., Burovski, E., Peterson, P., Weckesser, W., Bright, J., van der Walt, S. J., Brett, M., Wilson, J., Millman, K. J., Mayorov, N., Nelson, A. R. J., Jones, E., Kern, R., Larson, E., Carey, C. J., Polat, ˙I., Feng, Y., Moore, E. W., VanderPlas, J., Laxalde, D., Perktold, J., Cimrman, R., Henriksen, I., Quintero, E. A., Harris, C. R., Archibald, A. M., Ribeiro, A. H., Pedregosa, F., van Mulbregt, P. & SciPy 1.0 Contributors (2020) SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. Nature Methods, 17, 261–272. Scipy

[8] Hunter, J. D. (2007) Matplotlib: A 2D graphics environment. Computing in Science & Engineering, 9(3),90–95. Matplotlib
