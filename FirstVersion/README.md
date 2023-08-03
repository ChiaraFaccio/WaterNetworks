# WaterNetworks
This is the companion software of the article " Structural analysis of water networks" of  Benzi, M., Daidone, I., Faccio, C.,  Zanetti-Polzi, L.

## INSTALLATION:  
List of Python modules that you need to install on the pc:
- NumPy
- SciPy
- Matplotlib
- NetworkX

To install the code, you must download the two Python scripts and save them in the same folder. Then run the main_water.py script from the terminal. All the obtained results will be saved in the folder where the two scripts are located.

## USE:  
The code proposes a series of questions to be answered by the user. In the case of questions that involve choosing a file/folder, the code will automatically open a window and the user can select the file or folder. At first, the user is asked to import the .gro files, which contain the coordinates of the atoms. These can be in two formats:

1) a single .gro file with all the frames one after the other. In this case, the code generates a folder with the same file name, and within it, the algorithm creates a subfolder called "frames." Within the latter, there are individual frames, renamed with a number. If the file has already been split, or if a folder with this name already exists, the code returns an error.
2) if you have a folder with individual frames, the name of each frame should be an increasing number. Be careful, select the folder that directly contains the .gro files.

The code then asks:

- if you have undirected graphs (answer 'no' for directed graphs, version not fully implemented);
- in the case of undirected graphs, the distance threshold (in nm) to define two connected molecules;
- which type of atom to consider as a node (by default, oxygen);
- whether periodic boundary conditions are to be imposed (answer 2 for yes).

At this point, you can decide whether to calculate centrality measures, global metrics, whether you want a 2D or 3D graph of a frame, whether you need to write adjacency matrices in Matlab format, or whether you want to exit.

At the end of each computation, the results are copied into a .txt file and saved in the initial folder. The code asks again what you want to calculate.

## REFERENCES:  
The code makes heavy use of the NetworkX 2.6.3 module [1] in Python 3.7, the Python package NetworkSNS [2] which is distributed under the BSD 2-Clause License, the code [3] for imposing the periodic boundary conditions, and the open-source packages NumPy [4], SciPy [5], and Matplotlib [6].

[1] Hagberg, A., Chult, D. S. & Swart, P. (2008) Exploring Network Structure, Dynamics, and Function using NetworkX. In Varoquaux, G., Vaught, T. & Millman, J., editors, Proceedings of the 7th Python in Science Conference, pages 11 – 15, Pasadena, CA USA, [Networkx](https://networkx.org/)

[2] Bucci, A. (2021) NetworkSNS, https://github.com/alb95/NetworkSNS.git.

[3] Yang, Y. & Fergus, M. (2020) Pair distances with PBC, url = https://yangyushi.github.io/science/2020/11/02/pbc_py.html.

[4] Harris, C. R., Millman, K. J., van der Walt, S. J., Gommers, R., Virtanen, P., Cournapeau, D., Wieser, E.,
Taylor, J., Berg, S., Smith, N. J., Kern, R., Picus, M., Hoyer, S., van Kerkwijk, M. H., Brett, M., Haldane, A.,
del R´ıo, J. F., Wiebe, M., Peterson, P., Gerard-Marchant, P., Sheppard, K., Reddy, T., Weckesser, W., Abbasi, ´
H., Gohlke, C. & Oliphant, T. E. (2020) Array programming with NumPy. Nature, 585(7825), 357–362 [NumPy](https://numpy.org/)

[5] Virtanen, P., Gommers, R., Oliphant, T. E., Haberland, M., Reddy, T., Cournapeau, D., Burovski, E., Peterson,
P., Weckesser, W., Bright, J., van der Walt, S. J., Brett, M., Wilson, J., Millman, K. J., Mayorov, N., Nelson,
A. R. J., Jones, E., Kern, R., Larson, E., Carey, C. J., Polat, ˙I., Feng, Y., Moore, E. W., VanderPlas, J., Laxalde,
D., Perktold, J., Cimrman, R., Henriksen, I., Quintero, E. A., Harris, C. R., Archibald, A. M., Ribeiro, A. H.,
Pedregosa, F., van Mulbregt, P. & SciPy 1.0 Contributors (2020) SciPy 1.0: Fundamental Algorithms for
Scientific Computing in Python. Nature Methods, 17, 261–272. [Scipy](https://scipy.org/)

[6] Hunter, J. D. (2007) Matplotlib: A 2D graphics environment. Computing in Science & Engineering, 9(3),90–95. [Matplotlib](https://matplotlib.org/)
