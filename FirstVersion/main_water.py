#!/usr/bin/env python3

# Chiara Faccio, Scuola Normale Superiore, Pisa
# chiara.faccio@sns.it
# May 2022

from functions_water_network import *
import itertools
import tracemalloc

 
if __name__ == '__main__':  
    
    # We choose the file.gro
    print('Hello, this code allows you to apply graph theory to a box of water molecules.\nFirst of all, we select the file.gro from which to extract the data.\nYou have two possibilities:\n \
       (1) a single file.gro of trajectories, \n \
       (2) a folder containing several file.gro, one for each frame of the trajectory. \n ')
    fmt = 'What do you choose? (1 or 2) '
    
    res = input(fmt)
    if not res:
        sys.exit(1)
        
    try:
        if int(res) == 1:
            path, name_file = divide_file()
            list_number = [int(os.path.splitext(file)[0]) for file in os.listdir(path)]
            list_number.sort()
            list_path_files = [os.path.join(path, str(file)+'.gro') for file in list_number]
            
        elif int(res) == 2:
            print("Choose the folder : ")
            root = Tk()
            root.withdraw()
            path = filedialog.askdirectory()
            root.destroy()
            name_file = os.path.splitext(os.path.basename(path))[0]
            
            print("You have chosen the folder {} ".format(name_file))
            
            try:
            
                list_number = [int(os.path.splitext(file)[0]) for file in os.listdir(path)]
                file = os.listdir(path)[0]
                list_number.sort()
                list_path_files = [os.path.join(path, str(file)+'.gro') for file in list_number]
                
            except NameError:
                print('ERROR: The files .gro inside the folder must be renamed with numbers in ascending order!')
                sys.exit(1)
        
    except ValueError:
        print('Unrecognized keyword: {}'.format(res))
        sys.exit(1)
        
    fmt_und = 'Do you want undirected graph? (yes/no) '
    undirected = input(fmt_und)
        
    if undirected.lower() == 'yes':
    
        par = parameters_graph_undirected()
        dist = par['dist']             #threshold for edges
        name = par['name']             #name of nodes
        boundary = par['boundary']     #int values (1 = no boundary, 2 = 27 boxes)
        
        n_nodes_begin = extract_number_nodes(list_path_files[0], name)    # extract number of nodes in the box
        
        begin_loop = True
        
        while begin_loop:
        
            fmt = "What do you want to do?\n \
            1) to compute centrality measures \n \
            2) to compute other metrics \n \
            3) to plot a single box in 3D \n \
            4) to plot a single box in 2D \n \
            5) to save adjacency matrix in matlab format \n \
            6) exit \n "
            
            
            res = int(input(fmt))
            
            if res == 1:
                fmt2 = 'Which centrality measures do you want to compute? \n \
                1) node total communicability \n \
                2) degree \n \
                3) subgraph centrality \n \
                4) closeness centrality \n \
                5) betweenness centrality \n \
                6) Katz centrality \n \
                7) eigenvector centrality \n'
                
                res2 = input(fmt2)
                if int(res2) == 1:
                    fmt = 'Do you want normalized NTC? \n \
                    1) no \n \
                    2) yes, w.r.t. number of edges \n \
                    3) yes, w.r.t. total network communicability \n \
                    4) yes, w.r.t. exp( beta * lambda_1) \n'
                    
                    normalization = int(input(fmt))
                    
                    fmt = 'Value of beta = '
                    beta = float(input(fmt))
                    
                    while beta <= 0:
                        print('ERROR: beta is a positive parameter')
                        
                        fmt = 'Value of beta = '
                        beta = float(input(fmt))
                    
                    path2 = os.path.split(path)[0]
                    
                    with open(os.path.join(path2, 'NTC_beta_'+ str(beta)+ '_boundary_'+str(boundary)+'_normalization_'+str(normalization)+'.txt'), 'a+') as fobj:
                        fobj.write('Type of normalization = {} \n'.format(normalization))
                        fobj.write('\n')
                    
                    ex_time=[]
                    
                    for kk, file in enumerate(list_path_files):  
                    
                        t,_ = compute_NTC(file, beta, n_nodes_begin, path2, dist, name, boundary, normalization)
                        ex_time.append(t)
                    print('mean value execution times = ', np.mean(ex_time))                             
                                       
                elif int(res2) == 2:
                
                    path2 = os.path.split(path)[0]
                
                    with open(os.path.join(path2, 'DEGREE'+ '_boundary_' + str(boundary) +'.txt'), 'a+') as fobj:
                        fobj.write('\n')
                    
                    for kk, file in enumerate(list_path_files):  
                    
                        compute_degree(file, n_nodes_begin, path2, dist, name, boundary)
                    
                elif int(res2) == 3:
                    
                    fmt = 'Value of beta = '
                    beta = float(input(fmt))
                    
                    while beta <= 0:
                        print('ERROR: beta is a positive parameter')
                        
                        fmt = 'Value of beta = '
                        beta = float(input(fmt))

                    path2 = os.path.split(path)[0]
                    
                    with open(os.path.join(path2, 'SUBGRAPH_beta_'+ str(beta) + '_boundary_' + str(boundary) +'.txt'), 'a+') as fobj:
                        fobj.write('\n')
                        
                    ex_time=[]
                    
                    for kk, file in enumerate(list_path_files):  
                    
                        t, _ = compute_subgraph(file, beta, n_nodes_begin, path2, dist, name, boundary)
                        ex_time.append(t)
                    print('mean value execution times = ', np.mean(ex_time))                    
                                             
                elif int(res2) == 4:
                
                    path2 = os.path.split(path)[0]
                
                    with open(os.path.join(path2, 'CLOSENESS_boundary_' + str(boundary)+'.txt'), 'a+') as fobj:
                        fobj.write('\n')
                        
                    ex_time=[]
                    
                    for kk, file in enumerate(list_path_files):  
                    
                        t = compute_CLOSENESS(file, n_nodes_begin, path2, dist, name, boundary)
                        ex_time.append(t)
                    print('mean value execution times = ', np.mean(ex_time))
                                       
                elif int(res2) == 5:
                
                    path2 = os.path.split(path)[0]
                
                    with open(os.path.join(path2, 'BETWEENNESS_boundary_' + str(boundary)+'.txt'), 'a+') as fobj:
                        fobj.write('\n')
                        
                    ex_time=[]
                    
                    for kk, file in enumerate(list_path_files):  
                    
                        t = compute_BETWEENNESS(file, n_nodes_begin, path2, dist, name, boundary)
                        ex_time.append(t)
                    print('mean value execution times = ', np.mean(ex_time))
                    
                elif int(res2) == 6:
                    
                    fmt = 'Value of alpha = '
                    alpha = float(input(fmt))
                
                    path2 = os.path.split(path)[0]
                
                    with open(os.path.join(path2, 'KATZ_'+str(alpha)+ '_boundary_' + str(boundary)+'.txt'), 'a+') as fobj:
                        fobj.write('\n')
                        
                    ex_time=[]
                    
                    for kk, file in enumerate(list_path_files):  
                    
                        t = compute_KATZ(file, n_nodes_begin, path2, dist, name, boundary, alpha)
                        ex_time.append(t)
                    print('mean value execution times = ', np.mean(ex_time))
                    
                elif int(res2) == 7:
                
                    path2 = os.path.split(path)[0]
                
                    with open(os.path.join(path2, 'EIGENVECTOR'+ '_boundary_' + str(boundary) +'.txt'), 'a+') as fobj:
                        fobj.write('\n')
                    
                    for kk, file in enumerate(list_path_files):  
                    
                        compute_eigenvector(file, n_nodes_begin, path2, dist, name, boundary)
                    
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
                11) entropy of the graph with Laplacian (Von Neumann entropy) \n \
                12) mean value of min/max max_min_eigenvalues(n_nodes_begin, path2, boundary) \n \
                13) mean values of the density of the graph \n \
                14) mean values of the density of the boxes ({number of molecules}/{volume box}) \n \
                15) energy of the graph \n \
                16) ASPL and diameter \n \
                17) algebraic connectivity \n'

                res_metric = int(input(fmt))
                
                if res_metric == 1:
                
                    path2 = os.path.split(path)[0]
    
                    edges = []
                    cycles3 = []
                    cycles4 = []
                    cycles5 = []

                    for ii, file in enumerate(list_path_files):  
                        n_edges, F2, F5, F8 = count_cycles(file, n_nodes_begin, path2, dist, name, boundary) # F2 = number cycles of length 3; F5 = number cycles of length 4; F8 = number cycles of length 5
                            
                        edges.append(n_edges)
                        cycles3.append(F2)
                        cycles4.append(F5)
                        cycles5.append(F8)
                        
                    print('Average number of edges = {} \n'.format(np.mean(edges)))
                    print('Average number of cycles of length 3 = {} \n'.format(np.mean(cycles3)))
                    print('Average number of cycles of length 4 = {} \n'.format(np.mean(cycles4)))
                    print('Average number of cycles of length 5 = {} \n'.format(np.mean(cycles5)))
                                          
                elif res_metric == 2:
                
                    path2 = os.path.split(path)[0]
                    estrada_index_vect = []
                    
                    for ii, file in enumerate(list_path_files):  
                    
                        estrada_index_values = estrada_index(file, n_nodes_begin, path2, dist, name, boundary)
                        estrada_index_vect.append(estrada_index_values)
                        
                    print('Average number of estrada index = {} \n'.format(np.mean(estrada_index_vect)))
                    
                elif res_metric == 3:
                
                    path2 = os.path.split(path)[0]
                    watts_strogatz_vect = []
                    
                    for ii, file in enumerate(list_path_files):  
                        watts_strogatz_values = watts_strogatz(file, n_nodes_begin, path2, dist, name, boundary)
                        watts_strogatz_vect.append(watts_strogatz_values)
                        
                    print('Average number of  Watts-Strogatz coeff. clustering = {} \n'.format(np.mean(watts_strogatz_vect)))
                    
                elif res_metric == 4:
                
                    path2 = os.path.split(path)[0]
                    transitivity_index_vect = []
                    
                    for ii, file in enumerate(list_path_files):  
                        transitivity_index_values = transitivity_index(file, n_nodes_begin, path2, dist, name, boundary)
                        transitivity_index_vect.append(transitivity_index_values)
                        
                    print('Average number of transitivity index = {} \n'.format(np.mean(transitivity_index_vect)))
                    
                elif res_metric == 5:    
                    
                    path2 = os.path.split(path)[0]
                    
                    fmt_compute_eig = 'Have you already calculated the eigenvalues? (yes/no) (For each graph, the eigenvalues are written in a .txt file. All .txt files are inside a folder)  \n'
                
                    compute_eig = input(fmt_compute_eig)
                    
                    if compute_eig == 'yes':
                        bipartivity_for_files(n_nodes_begin, path2, boundary)
                        
                    elif compute_eig == 'no':
                    
                        bipartivity_vect = []
                        path3 = os.path.join(path2, 'eigenvalues' )
                        os.mkdir(path3)
                    
                        for ii, file in enumerate(list_path_files):  
                            bipartivity_values = bipartivity(file, n_nodes_begin, path2, path3,  dist, name, boundary,ii)
                            bipartivity_vect.append(bipartivity_values)
                        
                        print('Average number of bipartivity valuex = {} \n'.format(np.mean(bipartivity_vect)))
                            
                    else: 
                        print('Invalid answer')
                        continue                  
                    
                elif res_metric == 6:
                
                    path2 = os.path.split(path)[0]
                    assortativity_vect = []
                    path2_path1_vect = []
                    path3_path2_vect = []
                    
                    for ii, file in enumerate(list_path_files): 
                        assortativity_values, path2_path1, path3_path2  = assortativity(file, n_nodes_begin, path2, dist, name, boundary)
                        assortativity_vect.append(assortativity_values)
                        path2_path1_vect.append(path2_path1)
                        path3_path2_vect.append(path3_path2)
                        
                    print('Average number of degree assortativity coefficient = {} \n'.format(np.mean(assortativity_vect)))
                    print('Average number of P2/P1 = {} \n'.format(np.mean(path2_path1_vect)))
                    print('Average number of P3/P2 = {} \n'.format(np.mean(path3_path2_vect)))
                    
                elif res_metric == 7:

                    path2 = os.path.split(path)[0]
                    path3 = os.path.join(path2, 'eigenvalues' )
                    os.mkdir(path3)

                    for kk, file in enumerate(list_path_files): 
                        eigenvalues(file, n_nodes_begin, path3, dist, name, boundary, kk)
                        
                elif res_metric == 8:
                    fmt_frame = 'Which frame do you want to see? (Python starts counting from 0) \n'
                
                    frame = int(input(fmt_frame))
                    path2 = os.path.split(path)[0]
                    file  = list_path_files[frame]
                    sparsity_A(file, n_nodes_begin, path2, dist, name, boundary, frame)
                  
                elif res_metric == 9:
                    path2 = os.path.split(path)[0]
                    entropy_vect = []
                    
                    for ii, file in enumerate(list_path_files):  
                    
                        entropy_values = entropy_subgraph(file, n_nodes_begin, path2, dist, name, boundary)
                        entropy_vect.append(entropy_values)
                        
                    print('Average number of entropy = {} \n'.format(np.mean(entropy_vect)))
                    
                elif res_metric == 10:
                    path2 = os.path.split(path)[0]
                    entropy_vect = []
                    
                    for ii, file in enumerate(list_path_files):  
                    
                        entropy_values = entropy_TC(file, n_nodes_begin, path2, dist, name, boundary)
                        entropy_vect.append(entropy_values)
                        
                    print('Average number of entropy = {} \n'.format(np.mean(entropy_vect)))
                    
                elif res_metric == 11:
                    path2 = os.path.split(path)[0]
                    entropy_vect = []
                    
                    for ii, file in enumerate(list_path_files):  
                    
                        entropy_values = entropy_von_neumann(file, n_nodes_begin, path2, dist, name, boundary)
                        entropy_vect.append(entropy_values)
                        
                    print('Average number of entropy = {} \n'.format(np.mean(entropy_vect)))
                    
                elif res_metric == 12:
                
                    path2 = os.path.split(path)[0]
                    max_min_eigenvalues(n_nodes_begin, path2, boundary)
                    
                        
                elif res_metric == 13:
            
                    density_v = []
                    path2 = os.path.split(path)[0]

                    for kk, file in enumerate(list_path_files):  
                        density_ = density_graph(file, n_nodes_begin, path2, dist, name, boundary)
                        density_v.append(density_)
                        
                    print('mean density ', np.mean(density_v))
                    
                elif res_metric == 14: 
                    density_v = []
                    path2 = os.path.split(path)[0]

                    for kk, file in enumerate(list_path_files):  
                        density_ = density(file, n_nodes_begin, path2, dist, name, boundary)
                        density_v.append(density_)
                        
                    fig, axes = plt.subplots(nrows=1, ncols=1,  figsize=(20, 10)) 
                    plt.plot(density_v)
                    plt.xlabel('frames')
                    plt.ylabel('density')
                    plt.show()
                    fig.savefig(os.path.join(path2, 'density.png'))
                    plt.close()
                        
                    print('mean density ', np.mean(density_v))
                    
                    with open(os.path.join(path2, 'density.txt'), 'a+') as fobj:
                        for ii in range(len(density_v)):
                            fobj.write('{:f} \n'.format(density_v[ii]))
                        fobj.write('\n')
            
                        
                elif res_metric == 15:
                    path2 = os.path.split(path)[0]
                    
                    energy(file, n_nodes_begin, path2, dist, name, boundary)
                    
                elif res_metric == 16:
                
                    path2 = os.path.split(path)[0]
                    ASPL_vect = []
                    diam_vect = []
                    isolated_vect = []
                    
                    for ii, file in enumerate(list_path_files):  
                    
                        ASPL, diam, n_isolated = ASPL_diameter_isolated(file, n_nodes_begin, path2, dist, name, boundary)
                        ASPL_vect.append(ASPL)
                        diam_vect.append(diam)
                        isolated_vect.append(n_isolated)
                        
                        with open(os.path.join(path2, 'ASPL_diameter_n_isolated.txt'), 'a+') as fobj:
                            fobj.write('ASPL = {:f} \n'.format(ASPL))
                            fobj.write('diameter = {:f} \n'.format(diam))
                            fobj.write('n. isolated points = {:f} \n'.format(n_isolated))
                            fobj.write('\n')
                        
                    print('Average number of ASPL = {} \n'.format(np.mean(ASPL_vect)))
                    print('Average number of diam = {} \n'.format(np.mean(diam_vect)))
                    print('Average number of isolated points = {} \n'.format(np.mean(isolated_vect)))
                    
                elif res_metric == 17:
                
                    path2 = os.path.split(path)[0]
                    AG_vect = []
                    for ii, file in enumerate(list_path_files):  
                    
                        AG = algebraic_connectivity(file, n_nodes_begin, path2, dist, name, boundary)
                        AG_vect.append(AG)
                        
                        with open(os.path.join(path2, 'algebraic_connectivity.txt'), 'a+') as fobj:
                            fobj.write('AG = {:f} \n'.format(AG))
                            fobj.write('\n')
                        
                    print('Average number of algebraic connectivity = {} \n'.format(np.mean(AG)))     
                    
                    
                    
                    
                else:
                    print('Value not valid')
                    continue
                
                #endif              
                    
            elif res == 3:
            
                fmt_frame = 'Which frame do you want to plot? (Python starts counting from 0) \n'
                
                frame = int(input(fmt_frame))
                
                fmt = 'Do you want to color the nodes according to some measure of centrality already calculated? (yes / no) '
                
                res = input(fmt)
                
                if res.lower() == 'no':
                    file = list_path_files[frame]
                    _, _, _, coord, _ , _= get_graph(file, dist, name, 1,n_nodes_begin)   
                    plot_network_3d(coord[0:n_nodes_begin, 0:n_nodes_begin], dist, path, 'plot 3D')
                    
                elif res.lower() == 'yes':
                    print('Choose the file.txt with the values: \n')
                    path2 = os.path.split(path)[0]
                    root = Tk()
                    root.withdraw()
                    file_txt_centrality = filedialog.askopenfilename()
                    root.destroy()
                    
                    fmt_color = 'Chromatic scale or two colors to identify LDL and HDL : (chro/bicolor)   '
                    
                    res_color = input(fmt_color)
                    if res_color == 'bicolor':
                        fmt_bicolor = 'Threshold to identify LDL/HDL phase : '
                        try:
                            threshold_color = float(input(fmt_bicolor))
                            print('If the value of the node is <= threshold, then the node is in LDL phase (BLUE), otherwise it is in HDL phase (RED)')
                        except ValueError:
                            print('ERROR: value not valid')
                            continue
                    
                    color = []
                    values = []
                    
                    with open(file_txt_centrality, 'r') as data:
                        data.seek(0)
                        for ff in range(((n_nodes_begin)+1)*frame +2):
                            data.readline()

                        for i in range(n_nodes_begin):
                            cent = data.readline()
                            cent = float(cent[9:30])
                            values.append(cent)
                            if res_color == 'chro':
                                color.append(cent)                    
                            elif res_color == 'bicolor':
                                if cent <= threshold_color:
                                    color.append(0)
                                else:
                                    color.append(1)
                    print(max(values))
                    fmt = 'Do you want vmin and vmax by default? (yes/no) '
                    res = input(fmt)
                    if res.lower() == 'yes':
                        vmin = None
                        vmax = None
                    elif res.lower() == 'no':
                        fmt = 'vmin = '
                        vmin = float(input(fmt))
                        fmt = 'vmax = '
                        vmax = float(input(fmt))
                    else: 
                        print('Value not valid')
                        continue
                    
                    fmt = 'What name do you want to save the plot with?  '
                    name_img = input(fmt)                
                    
                    file = list_path_files[frame]
                    _, _, _, coord, _, _= get_graph(file, dist, name, 1,n_nodes_begin)   
                    plot_network_3d(coord[0:n_nodes_begin, 0:n_nodes_begin], dist, path, name_img, color = color, vmin = vmin, vmax =vmax)
                    
                else:
                    print('Value not valid')
                    continue
                    
            elif res == 4:
                fmt_frame = 'Which frame do you want to plot? (Python starts counting from 0) \n'
                
                frame = int(input(fmt_frame))
                
                fmt = 'Do you want to color the nodes according to some measure of centrality already calculated? (yes / no) '
                
                res = input(fmt)
                
                if res.lower() == 'no':
                    file = list_path_files[frame]
                    _, G, pos, _, _, _ = get_graph(file, dist, name, 1,n_nodes_begin)  
                    
                    range_nodes = range(0,n_nodes_begin)
                    
                    plot_network_2d(pos, G, range_nodes, path, 'plot 2D')
                    
                elif res.lower() == 'yes':
                    print('Choose the file.txt with the values: \n')
                    path2 = os.path.split(path)[0]
                    root = Tk()
                    root.withdraw()
                    file_txt_centrality = filedialog.askopenfilename()
                    root.destroy()
                    
                    fmt_color = 'Chromatic scale or two colors to identify LDL and HDL : (chro/bicolor)   '
                    
                    res_color = input(fmt_color)
                    if res_color == 'bicolor':
                        fmt_bicolor = 'Threshold to identify LDL/HDL phase : '
                        try:
                            threshold_color = float(input(fmt_bicolor))
                            print('If the value of the node is <= threshold, then the node is in LDL phase (BLUE), otherwise it is in HDL phase (RED)')
                        except ValueError:
                            print('ERROR: value not valid')
                            continue
                    
                    color = []
                    
                    with open(file_txt_centrality, 'r') as data:
                        data.seek(0)
                        for ff in range(((n_nodes_begin)+1)*frame +2):
                            data.readline()

                        for i in range(n_nodes_begin):
                            cent = data.readline()
                            cent = float(cent[7:])
                            if res_color == 'chro':
                                color.append(cent)                    
                            elif res_color == 'bicolor':
                                if cent <= threshold_color:
                                    color.append(0)
                                else:
                                    color.append(1)     
                    
                    fmt = 'Do you want vmin and vmax by default? (yes/no) '
                    res = input(fmt)
                    if res.lower() == 'yes':
                        vmin = None
                        vmax = None
                    elif res.lower() == 'no':
                        fmt = 'vmin = '
                        vmin = float(input(fmt))
                        fmt = 'vmax = '
                        vmax = float(input(fmt))
                    else: 
                        print('Value not valid')
                        continue
                    
                    fmt = 'What name do you want to save the plot with?  '
                    name_img = input(fmt)                
                    
                    file = list_path_files[frame]
                    _, G, pos, _, _ , _= get_graph(file, dist, name, 1, n_nodes_begin)  
                    range_nodes = range(0,n_nodes_begin)
                    
                    plot_network_2d(pos, G, range_nodes, path, name_img, color = color, vmin = vmin, vmax =vmax)
                
                
                    
                else:
                    print('Value not valid')
                    continue
                    
            elif res == 5:
                    
                path2 = os.path.split(path)[0]
                new_path = os.path.join(path2, 'matrices_matlab_format' )
                os.mkdir(new_path)
                
                for kk, file in enumerate(list_path_files):  
                    save_matrix_matlab_format(file, n_nodes_begin, new_path, dist, name, boundary, kk)

            
            elif res == 6:
                sys.exit(1)
                
            else:
                print('Value not valid')
                continue
                
    else:


        par = parameters_graph_directed()
        dist = par['dist']             #threshold for edges
        name = par['name']             #name of nodes
        boundary = par['boundary']     #int values (1 = no boundary, 2 = 27 boxes)
        
        n_nodes_begin = extract_number_nodes(list_path_files[0], name)    # extract number of nodes in the box
        
        begin_loop = True
        
        while begin_loop:
        
            fmt = 'What do you want to do?\n \
            1) to compute centrality measures \n \
            2) exit \n'
            
            res = int(input(fmt))
            
            if res == 1:
                fmt2 = 'Which centrality measures do you want to compute? \n \
                1) node total communicability \n \
                2) degree \n \
                3) node total communicability and degree \n '
                
                res2 = input(fmt2)
                if int(res2) == 1:
                    fmt = 'Do you want normalized NTC? \n \
                    1) no \n \
                    2) yes, w.r.t. number of edges \n'
                    
                    normalization = int(input(fmt))
                    
                    fmt = 'Value of beta = '
                    beta = float(input(fmt))
                    
                    while beta <= 0:
                        print('ERROR: beta is a positive parameter')
                        
                        fmt = 'Value of beta = '
                        beta = float(input(fmt))
                    
                    path2 = os.path.split(path)[0]
                    
                    with open(os.path.join(path2, 'directed_NTC_beta_'+ str(beta)+'.txt'), 'a+') as fobj:
                        fobj.write('Type of normalization = {} \n'.format(normalization))
                        fobj.write('# node                bc                      rc    \n')
                        


                    for kk, file in enumerate(list_path_files): 
                    
                        compute_NTC_directed(file, n_nodes_begin, path2, dist, name, boundary, normalization, beta)
                    
                elif int(res2) == 2:
                
                    path2 = os.path.split(path)[0]

                    with open(os.path.join(path2, 'directed_DEGREE.txt'), 'a+') as fobj:
                        fobj.write('# node                  bc                        rc    \n')
                        fobj.write('\n')

                    for kk, file in enumerate(list_path_files):
                        compute_directed_degree(file, n_nodes_begin, path2, dist, name, boundary)
                    
                elif int(res2) == 3:
                    fmt = 'Do you want normalized NTC? \n \
                    1) no \n \
                    2) yes, w.r.t. number of edges \n'
                    
                    normalization = int(input(fmt))
                    
                    fmt = 'Value of beta = '
                    beta = float(input(fmt))
                    
                    while beta <= 0:
                        print('ERROR: beta is a positive parameter')
                        
                        fmt = 'Value of beta = '
                        beta = float(input(fmt))
                    
                    path2 = os.path.split(path)[0]
                    
                    with open(os.path.join(path2, 'directed_NTC_beta_'+ str(beta)+'.txt'), 'a+') as fobj:
                        fobj.write('Type of normalization = {} \n'.format(normalization))
                        fobj.write('# node                bc                      rc    \n')
                        
                        
                    with open(os.path.join(path2, 'directed_DEGREE.txt'), 'a+') as fobj:
                        fobj.write('# node                  bc                        rc    \n')
                        fobj.write('\n')
                        

                    for kk, file in enumerate(list_path_files):
                    
                        compute_directed_NTC_degree(file, n_nodes_begin, path2, dist, name, boundary, normalization, beta)
                    
                
                else:
                    print('Value not valid')
                    continue

            
            elif res == 2:
                sys.exit(1)
                
            else:
                print('Value not valid')
                continue
                    

    sys.exit(1)
    
    
    
