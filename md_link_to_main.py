# This script contains the responses to the various "if"s of the "md_main_water.py" script

# Chiara Faccio, Scuola Normale Superiore, Pisa
# chiara.faccio@sns.it
# March 2023

from md_functions_water_network import *

def res1_1(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path):
    
    fmt = 'Do you want normalized NTC? \n \
    1) no \n \
    2) yes, w.r.t. number of edges \n \
    3) yes, w.r.t. total network communicability (TNC) \n \
    4) yes, w.r.t. exp( beta * lambda_1) \n \
    5) yes, w.r.t. number of nodes in the graph \n \
    6) yes, w.r.t. the TNC of a complete graph \n \
    7) yes, w.r.t. the average TNC  of 100 random graphs with the same number of nodes and edges as the original graph  \n\
    8) yes, w.r.t the average degree \n'
    
    normalization = int(input(fmt))
    
    fmt = 'Value of beta (beta > 0) = '
    beta = float(input(fmt))
    
    while beta <= 0:
        print('ERROR: beta is a positive parameter')
        
        fmt = 'Value of beta (beta > 0) = '
        beta = float(input(fmt))
    
    
    with open(os.path.join(path, 'NTC_beta_'+ str(beta)+ '_boundary_' + str(boundary)+'_normalization_'+str(normalization) + '_dist_'+ str(dist) +'.txt'), 'a+') as fobj:
        fobj.write('Type of normalization = {},  molecules to consider = {}, for each molecule which atom to consider as a node = {}, nodes_to_save = {}\n'.format(normalization, list_molecules, list_nodes, nodes_to_save))
        fobj.write('\n')
    
    Tc_mean = []
    
    for aaa in range(n_frames):
    
        tc = compute_NTC(U, aaa, beta, list_molecules, list_nodes, nodes_to_save, path, boundary, normalization, dist)
        Tc_mean.append(np.mean(tc))
                       
    print('mean value TC = ', np.mean(Tc_mean))  
    return
    
    
def res_1_2(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path):

    deg_m = []

    with open(os.path.join(path, 'DEGREE'+ '_boundary_' + str(boundary) +'dist_'+str(dist)+'.txt'), 'a+') as fobj:
        fobj.write('Molecules to consider = {}, for each molecule which atom to consider as a node = {}, nodes_to_save = {}\n'.format(list_molecules, list_nodes, nodes_to_save))
    
    for bbb in range(n_frames):
    
        deg = compute_degree(U, bbb, list_molecules, list_nodes, nodes_to_save, path, boundary, dist)
        
        deg_m.append(np.mean(deg))
        
    print('mean degree = ', np.mean(deg_m))
    
    return
                    

def res_1_3(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path):

    fmt = 'Value of beta (beta > 0) = '
    beta = float(input(fmt))
    
    while beta <= 0:
        print('ERROR: beta is a positive parameter')
        
        fmt = 'Value of beta (beta > 0) = '
        beta = float(input(fmt))
    
    with open(os.path.join(path, 'SUBGRAPH_beta_'+ str(beta) + '_boundary_' + str(boundary) +'.txt'), 'a+') as fobj:
        fobj.write('Molecules to consider = {}, for each molecule which atom to consider as a node = {}, nodes_to_save = {}\n'.format(list_molecules, list_nodes, nodes_to_save))
    
    for ccc in range(n_frames):
    
        _ = compute_subgraph(U, ccc, beta, list_molecules, list_nodes, nodes_to_save, path, boundary, dist)
    
    return
                    
 
def res_1_4(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) :
    with open(os.path.join(path, 'CLOSENESS_boundary_' + str(boundary)+'.txt'), 'a+') as fobj:
        fobj.write('Molecules to consider = {}, for each molecule which atom to consider as a node = {}, nodes_to_save = {}\n'.format(list_molecules, list_nodes, nodes_to_save))
    
    for ddd in range(n_frames):
    
        _ = compute_CLOSENESS(U, ddd, list_molecules, list_nodes, nodes_to_save, path, boundary, dist)
    
    return
    

def res_1_5(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) :
    with open(os.path.join(path, 'BETWEENNESS_boundary_' + str(boundary)+'.txt'), 'a+') as fobj:
        fobj.write('Molecules to consider = {}, for each molecule which atom to consider as a node = {}, nodes_to_save = {}\n'.format(list_molecules, list_nodes, nodes_to_save))
    
    for eee in range(n_frames):
    
        _ = compute_BETWEENNESS(U, eee, list_molecules, list_nodes, nodes_to_save, path, boundary, dist)
    
    return
    
    
def res_1_6(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) :
    fmt = 'Value of alpha = '
    alpha = float(input(fmt))

    with open(os.path.join(path, 'KATZ_'+str(alpha)+ '_boundary_' + str(boundary)+'.txt'), 'a+') as fobj:
        fobj.write('Molecules to consider = {}, for each molecule which atom to consider as a node = {}, nodes_to_save = {}\n'.format(list_molecules, list_nodes, nodes_to_save))
    
    for fff in range(n_frames):
    
        _ = compute_KATZ(U, fff, alpha, list_molecules, list_nodes, nodes_to_save, path, boundary, dist)
    
    return
    
    
def res_1_7(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) :
    with open(os.path.join(path, 'EIGENVECTOR'+ '_boundary_' + str(boundary) +'.txt'), 'a+') as fobj:
        fobj.write('Molecules to consider = {}, for each molecule which atom to consider as a node = {}, nodes_to_save = {}\n'.format(list_molecules, list_nodes, nodes_to_save))
    
    for ggg in range(n_frames):
    
        compute_eigenvector(U, ggg, list_molecules, list_nodes, nodes_to_save, path, boundary, dist)
        
    return
    
 
def res_1_8(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) :
    fmt = 'Value of beta (beta > 0) = '
    beta = float(input(fmt))
    
    while beta <= 0:
        print('ERROR: beta is a positive parameter')
        
        fmt = 'Value of beta (beta > 0) = '
        beta = float(input(fmt))
    
    
    with open(os.path.join(path, 'NBT_TC_beta_'+ str(beta)+ '_boundary_' + str(boundary)+'_dist_'+ str(dist) +'.txt'), 'a+') as fobj:
        fobj.write('Molecules to consider = {}, for each molecule which atom to consider as a node = {}, nodes_to_save = {}\n'.format(list_molecules, list_nodes, nodes_to_save))
        fobj.write('\n')
    
    for hhh in range(n_frames):
        compute_NBT_total_communicability(U, hhh, beta, list_molecules, list_nodes, nodes_to_save, path, boundary, dist)
        
    return
    
    
def res_2_1(U, list_molecules, list_nodes, boundary, dist, n_frames, path) :
    edges = []
    cycles3 = []
    cycles4 = []
    cycles5 = []

    for iii in range(n_frames):
        n_edges, F2, F5, F8 = count_cycles(U, iii, list_molecules, list_nodes, path, boundary, dist) # F2 = number cycles of length 3; F5 = number cycles of length 4; F8 = number cycles of length 5
            
        edges.append(n_edges)
        cycles3.append(F2)
        cycles4.append(F5)
        cycles5.append(F8)
        
    print('Average number of edges = {} \n'.format(np.mean(edges)))
    print('Average number of cycles of length 3 = {} \n'.format(np.mean(cycles3)))
    print('Average number of cycles of length 4 = {} \n'.format(np.mean(cycles4)))
    print('Average number of cycles of length 5 = {} \n'.format(np.mean(cycles5)))
    
    return
    

def res_2_2(U, list_molecules, list_nodes, boundary, dist, n_frames, path) :
    estrada_index_vect = []
    
    for jjj in range(n_frames):  
    
        estrada_index_values = estrada_index(U, jjj, list_molecules, list_nodes, path, boundary, dist)
        estrada_index_vect.append(estrada_index_values)
        
    print('Average number of estrada index = {} \n'.format(np.mean(estrada_index_vect)))
    
    return
    
    
def res_2_3(U, list_molecules, list_nodes, boundary, dist, n_frames, path) :

    watts_strogatz_vect = []
    
    for kkk in range(n_frames):   
        watts_strogatz_values = watts_strogatz(U, kkk, list_molecules, list_nodes, path, boundary, dist)
        watts_strogatz_vect.append(watts_strogatz_values)
        
    print('Average number of  Watts-Strogatz coeff. clustering = {} \n'.format(np.mean(watts_strogatz_vect)))
    
    return
    
    
def res_2_4(U, list_molecules, list_nodes,  boundary, dist, n_frames, path) :
                
    transitivity_index_vect = []
    
    for lll in range(n_frames):   
        transitivity_index_values = transitivity_index(U, lll, list_molecules, list_nodes, path, boundary, dist)
        transitivity_index_vect.append(transitivity_index_values)
        
    print('Average number of transitivity index = {} \n'.format(np.mean(transitivity_index_vect)))
    
    return

    
def res_2_5(U, list_molecules, list_nodes, boundary, dist, n_frames, path) :
                    
    fmt_compute_eig = 'Have you already calculated the eigenvalues? (yes/no) (For each graph, the eigenvalues are written in a .txt file. All .txt files are inside a folder)  \n'

    compute_eig = input(fmt_compute_eig)
    
    if compute_eig == 'yes':
        
        G, _ = get_graph_G(U, 0, list_molecules, list_nodes, boundary, dist)
        n = len(G)
        bipartivity_for_files(path, n)
        
    elif compute_eig == 'no':
    
        bipartivity_vect = []
        path3 = os.path.join(path, 'eigenvalues' )
        try:
            os.mkdir(path3)
        except:
            print('Unable to create folder because it already exists')
            return
           
    
        for mmm in range(n_frames):   
            bipartivity_values = bipartivity(U, mmm, list_molecules, list_nodes,  path, path3, boundary, dist)
            bipartivity_vect.append(bipartivity_values)
        
        print('Average number of bipartivity valuex = {} \n'.format(np.mean(bipartivity_vect)))
            
    else: 
        print('Invalid answer')
        
    return

    
def res_2_6(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) :
                
    assortativity_vect = []
    path2_path1_vect = []
    path3_path2_vect = []
    
    for nnn in range(n_frames):  
        assortativity_values, path2_path1, path3_path2  = assortativity(U, nnn, list_molecules, list_nodes, nodes_to_save, path, boundary, dist)
        assortativity_vect.append(assortativity_values)
        path2_path1_vect.append(path2_path1)
        path3_path2_vect.append(path3_path2)
        
    print('Average number of degree assortativity coefficient = {} \n'.format(np.mean(assortativity_vect)))
    print('Average number of P2/P1 = {} \n'.format(np.mean(path2_path1_vect)))
    print('Average number of P3/P2 = {} \n'.format(np.mean(path3_path2_vect)))
    
    return
    
    
def res_2_7(U, list_molecules, list_nodes, boundary, dist, n_frames, path) :

    path3 = os.path.join(path, 'eigenvalues' )
    try:
        os.mkdir(path3)
    except:
        print('Unable to create folder because it already exists')
        return

    for ooo in range(n_frames):   
        eigenvalues(U, ooo, list_molecules, list_nodes, path3, boundary, dist)
        
    return
    
    
def res_2_8(U, list_molecules, list_nodes, boundary, dist, n_frames, path) :
    fmt_frame = 'Which frame do you want to see? (Python starts counting from 0) \n'

    frame = int(input(fmt_frame))

    sparsity_A(U, frame, list_molecules, list_nodes, boundary, dist, path)
    return
    
    
def res_2_9(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) :
    

    entropy_vect = []
    fmt = 'Value of beta = '
    beta = float(input(fmt))
    
    for ppp in range(n_frames):   
    
        entropy_values = entropy_subgraph(U, ppp, beta, list_molecules, list_nodes, nodes_to_save, path, boundary, dist)
        entropy_vect.append(entropy_values)
        
    print('Average number of entropy (using subgraph) = {} \n'.format(np.mean(entropy_vect)))
        
    return
    
    
def res_2_10(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) :
                    
    entropy_vect = []
    fmt = 'Value of beta = '
    beta = float(input(fmt))
    
    for qqq in range(n_frames): 
    
        entropy_values = entropy_TC(U, qqq, beta, list_molecules, list_nodes, nodes_to_save, path, boundary, dist)
        entropy_vect.append(entropy_values)
        
    print('Average number of entropy (using TC) = {} \n'.format(np.mean(entropy_vect)))
    
    return
    
    
def res_2_11(U, list_molecules, list_nodes, boundary, dist, n_frames, path) :

    entropy_vect = []
    
    for rrr in range(n_frames):
    
        entropy_values = entropy_von_neumann(U, rrr, list_molecules, list_nodes, path, boundary, dist)
        entropy_vect.append(entropy_values)
        
    print('Average number of entropy (Von Neumann entropy) = {} \n'.format(np.mean(entropy_vect)))
    return

    
def res_2_12(U, list_molecules, list_nodes, boundary, dist, n_frames, path) :

    G, _ = get_graph_G(U, 0, list_molecules, list_nodes, boundary, dist)
    n = len(G)
    max_min_eigenvalues(n)  
    return


def res_2_13(U, list_molecules, list_nodes, boundary, dist, n_frames, path) :
            
    density_g = []

    for sss in range(n_frames):
        density_ = density_graph(U, sss, list_molecules, list_nodes, boundary, dist, path)
        density_g.append(density_)
        
    print('mean graph density ', np.mean(density_g))
    return


def res_2_14(U, list_molecules, list_nodes,  boundary, dist, n_frames, path) :
    density_v = []

    for ttt in range(n_frames):
        G, _ = get_graph_G(U, ttt, list_molecules, list_nodes, boundary, dist)
        n = len(G)
        density__ = density(U, ttt, list_molecules, list_nodes, n)
        density_v.append(density__)
        
    print('mean box density ', np.mean(density_v))
    
    with open(os.path.join(path, 'density.txt'), 'a+') as fobj:
        for ii in range(len(density_v)):
            fobj.write('{:f} \n'.format(density_v[ii]))
        fobj.write('\n')
        
    return


def res_2_15(U, list_molecules, list_nodes, boundary, dist, n_frames, path) :

    G, _ = get_graph_G(U, 0, list_molecules, list_nodes, boundary, dist)
    n = len(G)
    energy(path, n)
    
    return
                    
    
def res_2_16(U, list_molecules, list_nodes, boundary, dist, n_frames, path) :
                
    ASPL_vect = []
    diam_vect = []
    isolated_vect = []
    
    for uuu in range(n_frames):
    
        ASPL, diam, n_isolated = ASPL_diameter_isolated(U, uuu, list_molecules, list_nodes, boundary, dist)
        ASPL_vect.append(ASPL)
        diam_vect.append(diam)
        isolated_vect.append(n_isolated)
        
        with open(os.path.join(path, 'ASPL_diameter_n_isolated.txt'), 'a+') as fobj:
            fobj.write('ASPL = {:f} \n'.format(ASPL))
            fobj.write('diameter = {:f} \n'.format(diam))
            fobj.write('n. isolated points = {:f} \n'.format(n_isolated))
            fobj.write('\n')
        
    print('Average number of ASPL = {} \n'.format(np.mean(ASPL_vect)))
    print('Average number of diam = {} \n'.format(np.mean(diam_vect)))
    print('Average number of isolated points = {} \n'.format(np.mean(isolated_vect)))
    
    return
    

def res_2_17(U, list_molecules, list_nodes, boundary, dist, n_frames, path):
                
    AG_vect = []
    for kkkk in range(n_frames):
    
        AG = algebraic_connectivity(U, kkkk, list_molecules, list_nodes, boundary, dist)
        AG_vect.append(AG)
        
        with open(os.path.join(path, 'algebraic_connectivity.txt'), 'a+') as fobj:
            fobj.write('AG = {:f} \n'.format(AG))
            fobj.write('\n')
        
    print('Average number of algebraic connectivity = {} \n'.format(np.mean(AG)))     
    
    return
    
    
def res_3(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) :
            
    dict_nodes = range_nodes_to_save(U, 0, list_molecules, list_nodes, nodes_to_save)
    range_nodes = list(dict_nodes.keys())

    fmt_frame = 'Which frame do you want to plot? (Python starts counting from 0) \n'
    
    frame = int(input(fmt_frame))
    
    fmt = 'Do you want to color the nodes according to some measure of centrality already calculated? (yes / no) '
    
    res_c = input(fmt)
    
    if res_c.lower() == 'no':
        _, _, _, coord, _ , _= get_graph(U, frame, list_molecules, list_nodes, 1, dist, nodes_to_save)
        coord2 = coord[range_nodes,:]
        plot_network_3d(coord2, dist, path, 'plot 3D')
        
    elif res_c.lower() == 'yes':
        print('Choose the file.txt with the values: \n')
        
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
                return
        
        color = []
        values = []
        
        n_nodes = len(range_nodes)
        
        
        with open(file_txt_centrality, 'r') as data:
            data.seek(0)
            for ff in range(((n_nodes)+1)*frame +2):
                data.readline()

            for i in range(n_nodes):
                cent = data.readline()
                 
                nn_for_read = 0
                while cent[nn_for_read] != '>':
                    nn_for_read +=1

                cent = float(cent[nn_for_read+1 : nn_for_read+25])
                values.append(cent)
                if res_color == 'chro':
                    color.append(cent)                    
                elif res_color == 'bicolor':
                    if cent <= threshold_color:
                        color.append(0)
                    else:
                        color.append(1)
        print(min(values))
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
            return
        
        fmt = 'What name do you want to save the plot with?  '
        name_img = input(fmt)                
        
        _, _, _, coord, _, _= get_graph(U, frame, list_molecules, list_nodes, 1, dist, nodes_to_save) 
        coord2 = coord[range_nodes,:]           
        plot_network_3d(coord2, dist, path, name_img, color = color, vmin = vmin, vmax =vmax)
        
    else:
        print('Value not valid')
        
    return


def res_4(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) :
    fmt_frame = 'Which frame do you want to plot? (Python starts counting from 0) \n'
    
    frame = int(input(fmt_frame))
    
    fmt = 'Do you want to color the nodes according to some measure of centrality already calculated? (yes / no) '
    
    res_c2 = input(fmt)
    
    dict_nodes = range_nodes_to_save(U, 0, list_molecules, list_nodes, nodes_to_save)
    range_nodes = list(dict_nodes.keys())
    n_nodes = len(range_nodes)
    
    if res_c2.lower() == 'no':
    
        _, G, pos, _, _, _ = get_graph(U, frame, list_molecules, list_nodes, 1, dist, nodes_to_save)
        
        plot_network_2d(pos, G, range_nodes, path, 'plot 2D')
        
    elif res_c2.lower() == 'yes':
        print('Choose the file.txt with the values: \n')

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
                return
        
        color = []
        
        
        with open(file_txt_centrality, 'r') as data:
            
            data.seek(0)
            for ff in range(((n_nodes)+1)*frame +2):
                data.readline()

            for i in range(n_nodes):
                cent = data.readline()
                
                nn_for_read = 0
                while cent[nn_for_read] != '>':
                    nn_for_read +=1

                cent = float(cent[nn_for_read+1 : nn_for_read+25])

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
            return
        
        fmt = 'What name do you want to save the plot with?  '
        name_img = input(fmt)                
        
        _, G, pos, _, _ , _= get_graph(U, frame, list_molecules, list_nodes, 1, dist, nodes_to_save)
        
        plot_network_2d(pos, G, range_nodes, path, name_img, color = color, vmin = vmin, vmax =vmax)
    
    
        
    else:
        print('Value not valid')
        
    return
      
    
def res_5(U, list_molecules, list_nodes, boundary, dist, n_frames, path) :
                    
    new_path = os.path.join(path, 'matrices_matlab_format' )
    try:
        os.mkdir(new_path)
    except:
        print('Unable to create folder because it already exists')
        return
    
    for vvv in range(n_frames):  
        save_matrix_matlab_format(U, vvv, list_molecules, list_nodes, boundary, dist, new_path)
        
    return


def res_6_1(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) :
                
    fmt = 'Value of alpha (0 < alpha < 1/sigma, where sigma = max{rho(A^t), t = 0,1,...,T}) = '
    alpha = float(input(fmt))
    
    while alpha <= 0:
        print('ERROR: alpha is a positive parameter')
        
        fmt = 'Value of alpha = '
        alpha = float(input(fmt))
    
    katz_dynamic_graph(U, list_molecules, list_nodes, boundary, dist, n_frames, path, alpha, nodes_to_save)
    
    return
    
    
def res_6_2(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) :
                
    aggregated_degree(U, list_molecules, list_nodes, boundary, dist, n_frames, path, nodes_to_save)
    return


def res_7_1(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path, str_weight, weight_edges) :
            
    fmt_b = 'Value of beta for TC (beta > 0) = '
    beta = float(input(fmt_b))
    
    while beta <= 0:
        print('ERROR: beta is a positive parameter')
        
        fmt_b = 'Value of beta (beta > 0) = '
        beta = float(input(fmt_b))
    
    with open(os.path.join(path, 'weighted_NTC_beta_'+ str(beta)+'_undirected_boundary_'+str(boundary)+'_dist_'+str(dist)+'_weight_'+str(str_weight)+'.txt'), 'a+') as fobj:
        fobj.write('\n')
        
    with open(os.path.join(path, 'weighted_DEGREE_undirected_boundary_'+str(boundary)+'_dist_'+str(dist)+'_weight_'+str(str_weight)+'.txt'), 'a+') as fobj:
        fobj.write('\n')
        
    
    for cccc in range(n_frames):
        compute_weighted_NTC_degree(U, cccc, beta, list_molecules, list_nodes, nodes_to_save, path, boundary, dist, weight_edges= weight_edges)
        
    return
    
    
def res_7_2(U, list_molecules, list_nodes, boundary, dist, n_frames, path, str_weight, weight_edges) :
                
    fmt_b = 'Value of beta = '
    beta = float(input(fmt_b))
    
    fmt_method = 'What method do you use to compute the Estrada index? (1 = eigenvalues, 2 = subgraph centrality values) '
    method = int(input(fmt_method))
    
    EE_vect = []
    
    for dddd in range(n_frames):  
        EE = compute_weighted_Estrada_index(U, dddd, beta, list_molecules, list_nodes, path, boundary, dist, method, weight_edges = weight_edges)
        EE_vect.append(EE)
        
    print('Average weighted EE = ', np.mean(EE_vect))
    
    return
    
    
def res_7_3(U, list_molecules, list_nodes, boundary, dist, n_frames, path, str_weight, weight_edges) :
                    
    entropy_vect = []
    
    for eeee in range(n_frames): 
    
        entropy_values = weighted_entropy_TC(U, eeee, list_molecules, list_nodes, dist, path, boundary, weight_edges = weight_edges)
        entropy_vect.append(entropy_values)
        
    print('Average number of entropy (using TC) = {} \n'.format(np.mean(entropy_vect)))
        
    return
    
    
def res_7_4(U, list_molecules, list_nodes, boundary, dist, n_frames, path, str_weight, weight_edges) :
    entropy_vect = []
    
    for ffff in range(n_frames):  
    
        entropy_values = weighted_entropy_von_neumann(U, ffff, list_molecules, list_nodes, dist, path, boundary, weight_edges = weight_edges)
        entropy_vect.append(entropy_values)
        
    print('Average number of entropy (Von Neumann entropy) = {} \n'.format(np.mean(entropy_vect)))
    
    return
    
    
def res_10_1(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path):

    fmt = 'Do you want normalized TC? \n \
    1) no \n \
    2) yes, w.r.t. number of edges \n'
    
    normalization = int(input(fmt))
    
    fmt = 'Value of beta (beta > 0) = '
    beta = float(input(fmt))
    
    while beta <= 0:
        print('ERROR: beta is a positive parameter')
        
        fmt = 'Value of beta (beta > 0) = '
        beta = float(input(fmt))
    
    
    with open(os.path.join(path, 'directed_TC_beta_'+ str(beta)+ '_boundary_' + str(boundary)+'_normalization_'+str(normalization) + '_dist_'+ str(dist) +'.txt'), 'a+') as fobj:
        fobj.write('Type of normalization = {},  molecules to consider = {}, for each molecule which atom to consider as a node = {}, nodes_to_save = {}\n'.format(normalization, list_molecules, list_nodes, nodes_to_save))
        fobj.write('\n')
    
    for gggg in range(n_frames):
    
        compute_TC_directed(U, gggg, beta, list_molecules, list_nodes, nodes_to_save, path, boundary, normalization, dist)


    return
    
    
def res_10_2(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) :
                
    with open(os.path.join(path, 'directed_DEGREE'+ '_boundary_' + str(boundary) +'dist_'+str(dist)+'.txt'), 'a+') as fobj:
        fobj.write('Molecules to consider = {}, for each molecule which atom to consider as a node = {}, nodes_to_save = {}\n'.format(list_molecules, list_nodes, nodes_to_save))
        fobj.write('\n')
    
    for hhhh in range(n_frames):
    
        compute_directed_degree(U, hhhh, list_molecules, list_nodes, nodes_to_save, path, boundary, dist)
        
    return
    
    
def res_10_3(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path):
                
    with open(os.path.join(path, 'directed_EIG'+ '_boundary_' + str(boundary) +'dist_'+str(dist)+'.txt'), 'a+') as fobj:
        fobj.write('Molecules to consider = {}, for each molecule which atom to consider as a node = {}, nodes_to_save = {}\n'.format(list_molecules, list_nodes, nodes_to_save))
        fobj.write('\n')
    
    for hhhh in range(n_frames):
    
        compute_directed_eigenvector(U, hhhh, list_molecules, list_nodes, nodes_to_save, path, boundary, dist)
        
    return

    
def res_10_4(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path):
    
    with open(os.path.join(path, 'directed_HITS'+ '_boundary_' + str(boundary) +'dist_'+str(dist)+'.txt'), 'a+') as fobj:
        fobj.write('Molecules to consider = {}, for each molecule which atom to consider as a node = {}, nodes_to_save = {}\n'.format(list_molecules, list_nodes, nodes_to_save))
        fobj.write('\n')
    
    for hh in range(n_frames):
    
        compute_directed_HITS(U, hh, list_molecules, list_nodes, nodes_to_save, path, boundary, dist)
    

    return

    
def res_10_5(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path):
    
    with open(os.path.join(path, 'directed_gTC'+ '_boundary_' + str(boundary) +'dist_'+str(dist)+'.txt'), 'a+') as fobj:
        fobj.write('Molecules to consider = {}, for each molecule which atom to consider as a node = {}, nodes_to_save = {}\n'.format(list_molecules, list_nodes, nodes_to_save))
        fobj.write('\n')
    
    for hh in range(n_frames):
    
        compute_directed_gTC(U, hh, list_molecules, list_nodes, nodes_to_save, path, boundary, dist)
    

    return
    
    
def res_10_6(U, list_molecules, list_nodes, nodes_to_save, boundary, dist, n_frames, path) :
    fmt = 'Do you want normalized NTC? \n \
    1) no \n \
    2) yes, w.r.t. number of edges \n'
    
    normalization = int(input(fmt))
    
    fmt = 'Value of beta (beta > 0) = '
    beta = float(input(fmt))
    
    while beta <= 0:
        print('ERROR: beta is a positive parameter')
        
        fmt = 'Value of beta (beta > 0) = '
        beta = float(input(fmt))
    
    
    with open(os.path.join(path, 'directed_bTC_beta_'+ str(beta)+ '_boundary_' + str(boundary)+'_normalization_'+str(normalization) + '_dist_'+ str(dist) +'.txt'), 'a+') as fobj:
        fobj.write('Type of normalization = {},  molecules to consider = {}, for each molecule which atom to consider as a node = {}, nodes_to_save = {}\n'.format(normalization, list_molecules, list_nodes, nodes_to_save))
        fobj.write('\n')
    
    for gggg in range(n_frames):
    
        compute_bTC_directed(U, gggg, beta, list_molecules, list_nodes, nodes_to_save, path, boundary, normalization, dist)
        
    return
    
    