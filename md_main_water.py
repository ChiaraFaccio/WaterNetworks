#!/usr/bin/env python3

# Chiara Faccio, Scuola Normale Superiore, Pisa
# chiara.faccio@sns.it
# March 2023

from md_link_to_main import *

 
if __name__ == '__main__':  
    
    # We choose the file.gro
    print('Hello, this code allows you to apply graph theory to boxes of water molecules. \nThe code uses MDAnalysis tools to access data in molecular dynamics trajectories. \nTo begin with, the topology and trajectory files need to be selected. The topology file can be in PSF, PDB, CRD, or GRO format. \nIf no trajectory files are provided, only a structure is loaded. \nThe supported trajectory file formats include single frames like \
 PDB, CRD, and GRO, and time-series data like DCD, XTC, TRR, and XYZ.\n ')
    
    path, topol = select_topology_file()
    
    answ = input('Do you want to upload a trajectory file? (yes/no) ')
    
    if answ == 'yes':
        traj = select_trajectory_file()
        U = mda.Universe(topol, traj)
    else:
        U = mda.Universe(topol)
        
    fmt_und = 'Do you want undirected graph? (yes/no) (directed graphs are implemented only in the case of orthogonal boxes) '
    undirected = input(fmt_und)
        
    if undirected.lower() == 'yes':
    
        par = parameters_graph_undirected(topol)
        dist = par['dist']                       # threshold for edges
        boundary = par['boundary']               # (1 = no boundary, 2 = pbc)
        list_molecules = par['list_molecules']   # which molecules consider
        list_nodes = par['list_nodes']           # list_nodes = for each molecule which atom to consider as a node
        nodes_to_save = par['nodes_to_save']     # nodes_to_save = the name of the molecule for which you would like to save the values of centrality measures of its nodes
        
        n_frames = len(U.trajectory)
        
        begin_loop = True
        
        while begin_loop:
        
            fmt = "What do you want to do?\n \
            1) Centrality measure computation \n \
            2) Computation of other metrics \n \
            3) to plot a single box in 3D \n \
            4) to plot a single box in 2D \n \
            5) to save adjacency matrix in MATLAB format \n \
            6) Dynamic graph metrics \n \
            7) Weighted graph metrics \n \
            8) exit \n "
            
            
            res = int(input(fmt))
            
            if res == 1:
                fmt2 = 'Which centrality measures do you want to compute? \n \
                1) node total communicability \n \
                2) degree \n \
                3) subgraph centrality \n \
                4) closeness centrality \n \
                5) betweenness centrality \n \
                6) Katz centrality \n \
                7) eigenvector centrality \n \
                8) non-backtracking total communicability \n '
                
                res2 = input(fmt2)
                if int(res2) == 1:
                
                    res1_1(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path)
                                       
                elif int(res2) == 2:
                    
                    res_1_2(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path)
                    
                elif int(res2) == 3:
                    
                    res_1_3(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path)                 
                                             
                elif int(res2) == 4:
                    
                    res_1_4(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) 
                                       
                elif int(res2) == 5:
                
                    res_1_5(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) 
                    
                elif int(res2) == 6:
                    
                    res_1_6(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) 
                    
                elif int(res2) == 7:
                    
                    res_1_7(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) 
                    
                elif int(res2) == 8:
                    
                    res_1_8(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) 
                    
            elif res == 2:
                fmt = 'What do you want to do? \n \
                1) to count edges and cycles of length 3, 4 and 5 \n \
                2) to compute Estrada index \n \
                3) to compute Watts-Strogatz clustering coefficient \n \
                4) to compute transitivity index \n \
                5) to compute bipartivity measures \n \
                6) to compute the degree assortativity coefficient \n \
                7) to compute eigenvalues \n \
                8) sparsity pattern of a frame \n \
                9) entropy of the graph with subgraph \n \
                10) entropy of the graph with TC \n \
                11) entropy of the graph with Laplacian graph (von Neumann entropy) \n \
                12) mean value of min/max max_min_eigenvalues (eigenvalues already calculated) \n \
                13) mean values of the density of the graph \n \
                14) mean values of the density of the boxes ({number of molecules}/{volume box}) \n \
                15) energy of the graph (eigenvalues already calculated) \n \
                16) ASPL and diameter \n \
                17) algebraic connectivity \n'

                res_metric = int(input(fmt))
                
                if res_metric == 1:
                    
                    res_2_1(U, list_molecules, list_nodes, boundary, dist, n_frames, path) 
                                          
                elif res_metric == 2:
                    
                    res_2_2(U, list_molecules, list_nodes, boundary, dist, n_frames, path) 
                    
                elif res_metric == 3:
                
                    res_2_3(U, list_molecules, list_nodes, boundary, dist, n_frames, path) 
                    
                elif res_metric == 4:
                
                    res_2_4(U, list_molecules, list_nodes, boundary, dist, n_frames, path) 
                                        
                elif res_metric == 5:    
                
                    res_2_5(U, list_molecules, list_nodes, boundary, dist, n_frames, path) 
                    
                elif res_metric == 6:
                
                    res_2_6(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) 

                elif res_metric == 7:

                    res_2_7(U, list_molecules, list_nodes, boundary, dist, n_frames, path) 
                        
                elif res_metric == 8:
                
                    res_2_8(U, list_molecules, list_nodes, boundary, dist, n_frames, path) 
                  
                elif res_metric == 9:
                
                    res_2_9(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) 
                    
                elif res_metric == 10:
                
                    res_2_10(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) 
                    
                elif res_metric == 11:
                
                    res_2_11(U, list_molecules, list_nodes, boundary, dist, n_frames, path) 
                    
                elif res_metric == 12:
                
                    res_2_12(U, list_molecules, list_nodes, boundary, dist, n_frames, path) 
                        
                elif res_metric == 13:
                
                    res_2_13(U, list_molecules, list_nodes, boundary, dist, n_frames, path) 
                    
                elif res_metric == 14: 
                
                    res_2_14(U, list_molecules, list_nodes, boundary, dist, n_frames, path) 
                                   
                elif res_metric == 15:
                
                    res_2_15(U, list_molecules, list_nodes, boundary, dist, n_frames, path) 
                    
                elif res_metric == 16:
                
                    res_2_16(U, list_molecules, list_nodes, boundary, dist, n_frames, path) 
                    
                elif res_metric == 17:
                
                    res_2_17(U, list_molecules, list_nodes, boundary, dist, n_frames, path)
                    
                else:
                    print('Value not valid')
                    continue
                
                #endif              
                    
            elif res == 3:
            
                res_3(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) 
                
            elif res == 4:
            
                res_4(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) 
                    
            elif res == 5:
            
                res_5(U, list_molecules, list_nodes, boundary, dist, n_frames, path) 

            elif res == 6:
                
                fmt_1 = 'What do you want to do?\n \
                1) to compute dynamic communicability \n \
                2) to compute aggregated degree \n'
                
                resp = int(input(fmt_1))
                
                if resp == 1:
                
                    res_6_1(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) 
                    
                elif resp == 2:
                
                    res_6_2(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) 
                    
                elif resp == 3:
                    
                    res_6_3(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) 
                                     
                
            elif res == 7:    
            
                fmt_type = 'How are edge weights defined? \n \
                1) w(i,j) = 1/d(i,j)  \n \
                2) w(i,j) = 1/d^2(i,j)  \n \
                3) w(i,j) = 1/d^3(i,j)  \n \
                4) w(i,j) = 1/sqrt(d(i,j)) \n '
                
                weight_edges = int(input(fmt_type))  
                if weight_edges == 1:
                    str_weight = '1fracd'
                elif weight_edges ==2:
                    str_weight = '1fracd2'
                elif weight_edges ==3:
                    str_weight = '1fracd3'
                elif weight_edges ==4:
                    str_weight = '1fracsqrt(d)'
                
                
                fmt_7 = 'What do you want to do? \n \
                1) to compute degree e TC \n \
                2) to compute Estrada index \n \
                3) entropy of the graph with TC \n \
                4) entropy of the graph with graph Laplacian (von Neumann entropy) \n'
                
                res_7 = int(input(fmt_7))
              
                if res_7 == 1:
                
                    res_7_1(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path, str_weight, weight_edges) 
                        
                elif res_7 == 2:
                
                    res_7_2(U, list_molecules, list_nodes, boundary, dist, n_frames, path, str_weight, weight_edges) 
                    
                elif res_7 == 3:
                
                    res_7_3(U, list_molecules, list_nodes, boundary, dist, n_frames, path, str_weight, weight_edges) 
                    
                elif res_7 == 4:
                
                    res_7_4(U, list_molecules, list_nodes, boundary, dist, n_frames, path, str_weight, weight_edges) 
            
            
            elif res == 8:
                sys.exit(1)
                
            else:
                print('Value not valid')
                continue
                
    else:

        par = parameters_graph_directed(topol)
        dist = par['dist']                       # threshold for edges
        boundary = par['boundary']               # (1 = no boundary, 2 = pbc)
        list_molecules = par['list_molecules']   # which molecules consider
        list_nodes = par['list_nodes']           # list_nodes = for each molecule which atom to consider as a node
        nodes_to_save = par['nodes_to_save']     # nodes_to_save = the name of the molecule for which you would like to save the values of centrality measures of its nodes
        
        n_frames = len(U.trajectory)
        
        begin_loop = True
        
        while begin_loop:
        
            fmt = 'What do you want to do?\n \
            1) Centrality measures computation \n \
            2) exit \n'
            
            res = int(input(fmt))
            
            if res == 1:
                fmt2 = 'Which centrality measures do you want to compute? \n \
                1) in- and out-degree  \n \
                2) left and right eigenvector centrality of A \n \
                3) HITS \n \
                4) gTC \n '
                
                res2 = input(fmt2) 
                    
                if int(res2) == 1:
                
                    res_10_2(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) 
                    
                elif int(res2) == 2:
                
                    res_10_3(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) 
                    
                elif int(res2) == 3:
                
                    res_10_4(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) 
                    
                elif int(res2) == 4:
                
                    res_10_5(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) 
                    
                else:
                    print('Value not valid')
                    continue

            
            elif res == 2:
                sys.exit(1)
                
            else:
                print('Value not valid')
                continue
                    

    sys.exit(1)
    
    
    
