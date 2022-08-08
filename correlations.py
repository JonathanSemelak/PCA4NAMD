import numpy as np
import MDAnalysis as md
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.kernel_ridge import KernelRidge
from scipy.signal import savgol_filter
from functions import*
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV

def covar_coord_arrh(rfile,topfile,e0file,e1file,temperature,alpha,framespertraj,filter,nframes_animation,amplitude):
    kB=0.000003173 #kB in hartree
    kt=kB*temperature
    print('Reading coordinates')
    natoms, nframes = get_natoms_nframes(rfile)
    print("natoms:", natoms)
    print("nframes:", nframes)
    names, resids = get_names_and_resids(topfile,natoms)
    atomic_numbers = get_atomic_numbers(names)
    print('Loading trajectory')
    r = get_r_md(topfile,rfile,natoms,nframes)
    # r = get_r(rfile,natoms,nframes) #this is another way to load the trajectory but without MDAnalysis lib
    print('Loading Energies evolution')
    e0 = np.loadtxt(e0file)
    e1 = np.loadtxt(e1file)
    energy = np.column_stack((e0,e1))
    print('Calculating distance matrix')
    dmat = get_distance_matrix(r)
    print('Processing')
    energy = flip(energy,framespertraj)
    energy_diff_raw = energy[:,1]-energy[:,0]
    energy_diff_smoothened = savgol_filter(energy[:,1]-energy[:,0], filter, 3)
    energy_diff_derivative = get_derivative(energy_diff_smoothened)
    np.savetxt("energy_diff_raw.dat",energy_diff_raw,fmt='%10.5f')
    np.savetxt("energy_diff_smoothened.dat",energy_diff_smoothened,fmt='%10.5f')
    np.savetxt("energy_diff_derivative.dat",energy_diff_derivative,fmt='%10.5f')


    print('Calculating averages')
    r_av = np.array([np.mean(r,axis=0)])
    # r = np.append(r,r_av,axis=0) #average coordinates could be appended as last frame

    print('Evaluating displacement vector (r - <r>)')
    disp = r[:]- r_av

    print('Calculating the relaxation pathway (c)')
    # The relaxation (c) pathway is a 3*natom-dimensional vector
    # The i-th component is ci=<dispi*arrh_termi>/sqrt(<dispi**2><arrh_termi**2>)
    # Where:
    # dispi = ri-<ri>
    # arrh_termi = sign(-enerdiff)exp(-|enerdiff|/alphakt) - <sign(-enerdiff)exp(-|enerdiff|/alphakt)>
    # i is an atomic coordinate, not a frame

    avgenerdiff=np.average(energy_diff_smoothened) #not used

    #calculates arrh_term and <arrh_term> (this is the same for every coordinate of c)
    arrh_term=np.zeros(nframes)
    for i in range(0,nframes):
        temp=energy_diff_smoothened[i]
        arrh_term[i]=np.sign(-temp)*np.exp(-abs((temp)/(alpha*kt)))
    arrh_term_av=np.average(arrh_term)
    np.savetxt("arrh_term.dat",arrh_term,fmt='%10.5f')

    #writes P(|enerdiff|)=|enerdiff|*exp(-|enerdiff|/alphakt)/Qenerdiff
    #for windows of size window_size
    #where Q_enerdiff=sum[exp(-|enerdiff|/alphakt)] (sum over the frames in the window)
    #this is only used to check the parameter alpha adecuacy
    P_enerdiff, P_enerdiff_by_windows = get_P_ener_diffs(energy_diff_smoothened,alpha,kt,nframes,50)
    np.savetxt('P_enerdiff.dat',P_enerdiff)
    np.savetxt('P_enerdiff_by_windows.dat',P_enerdiff_by_windows,fmt='%10.5f')

    #calculates c, first in an array shape
    c_array=np.zeros((nframes,natoms,3))
    for i in range(0,nframes):
        c_array[i]=disp[i]*(arrh_term[i]-arrh_term_av)
    c_array=np.mean(c_array,axis=0)

    #calculates the normalization factor (i.e. sqrt(<dispi**2><arrh_termi**2>))
    arrh_term2=np.square(arrh_term)
    arrh_term2_av=np.average(arrh_term2)
    disp2=np.square(disp)
    disp2_av=np.mean(disp2,axis=0)
    norm_factor=np.sqrt(disp2_av*arrh_term2_av)

    #normalizes c and reshape it in a vector-like shape
    c_array=c_array/norm_factor #carray esta bien
    c_vector=np.reshape(c_array, 3*natoms)
    c_vector_norm=np.sqrt(np.sum(np.square(c_vector)))
    c_vector=c_vector/c_vector_norm
    np.savetxt('coordinate.dat',c_vector,fmt='%10.5f')

    #computes centrality
    centrality=get_centrality(c_vector,natoms)
    np.savetxt('centrality.dat',centrality,fmt='%10.5f')

    #writes reordered-mode
    write_nmd('input-writer.nmd',names,natoms,resids,r,c_vector)

    #writes energy-diff mode animation
    write_animation('animode-energy-diff.xyz',atomic_numbers,r_av,c_vector,natoms,nframes_animation,amplitude)

    #writes md-energy-diff projection
    wt=np.zeros(nframes)
    for i in range(0,nframes):
        disp_vector=np.reshape(disp[i], 3*natoms)
        wt[i]=np.dot(c_vector,disp_vector)
    write_projection('projection-md-energy-diff.xyz',atomic_numbers,r_av,disp,c_vector,0,natoms,nframes,wt)


    print('-------------------------------------------')
    print('The following files have been written:')
    print('energy_diff_raw.dat')
    print('energy_diff_smoothened.dat')
    print('energy_diff_derivative.dat')
    print('arrh_term.dat')
    print('coordinate.dat')
    print('centrality.dat')
    print('input-writer.nmd')
    print('animode-energy-diff.xyz')
    print('projection-md-energy-diff.xyz')
    print('-------------------------------------------')


def covar_mulliken_arrh(mfile,topfile,e0file,e1file,temperature,alpha,framespertraj,filter):
    kB=0.000003173 #kB in hartree
    kt=kB*temperature
    print('Reading coordinates')
    natoms, nframes = get_natoms_nframes_mulliken(mfile)
    names, resids = get_names_and_resids(topfile,natoms)
    print("natoms:", natoms)
    print("nframes:", nframes)
    print('Loading Mulliken charges evolution')
    q = get_mull(mfile,natoms,nframes)
    print('Loading Energies evolution')
    e0 = np.loadtxt(e0file)
    e1 = np.loadtxt(e1file)
    energy = np.column_stack((e0,e1))
    print('Processing')
    energy = flip(energy,framespertraj)
    energy_diff_raw = energy[:,1]-energy[:,0]
    energy_diff_smoothened = savgol_filter(energy[:,1]-energy[:,0], filter, 3)
    energy_diff_derivative = get_derivative(energy_diff_smoothened)
    np.savetxt("energy_diff_raw.dat",energy_diff_raw,fmt='%10.5f')
    np.savetxt("energy_diff_smoothened.dat",energy_diff_smoothened,fmt='%10.5f')
    np.savetxt("energy_diff_derivative.dat",energy_diff_derivative,fmt='%10.5f')

    print('Calculating averages')
    q_av = np.array([np.mean(q,axis=0)])
    # r = np.append(r,r_av,axis=0) #average coordinates could be appended as last frame

    print('Evaluating displacement vector (q - <q>)')
    disp = q[:]- q_av

    print('Calculating the relaxation pathway (c)')
    # The relaxation (c) pathway is a natom-dimensional vector
    # The i-th component is ci=<dispi*arrh_termi>/sqrt(<dispi**2><arrh_termi**2>)
    # Where:
    # dispi = qi-<qi>
    # arrh_termi = sign(-enerdiff)exp(-|enerdiff|/alphakt) - <sign(-enerdiff)exp(-|enerdiff|/alphakt)>
    # i is an mulliken (mulliken) coordinate, not a frame
    avgenerdiff=np.average(energy_diff_smoothened) #not used

    #calculates arrh_term and <arrh_term> (this is the same for every coordinate of c)
    arrh_term=np.zeros(nframes)
    for i in range(0,nframes):
        temp=energy_diff_smoothened[i]
        arrh_term[i]=np.sign(-temp)*np.exp(-abs((temp)/(alpha*kt)))
    arrh_term_av=np.average(arrh_term)
    np.savetxt("arrh_term.dat",arrh_term,fmt='%10.5f')

    #writes P(|enerdiff|)=|enerdiff|*exp(-|enerdiff|/alphakt)/Qenerdiff
    #for windows of size window_size
    #where Q_enerdiff=sum[exp(-|enerdiff|/alphakt)] (sum over the frames in the window)
    #this is only used to check the parameter alpha adecuacy
    P_enerdiff, P_enerdiff_by_windows = get_P_ener_diffs(energy_diff_raw,alpha,kt,nframes,50)
    np.savetxt('P_enerdiff.dat',P_enerdiff)
    np.savetxt('P_enerdiff_by_windows.dat',P_enerdiff_by_windows,fmt='%10.5f')

    #calculates c, first in an array shape
    c_array=np.zeros((nframes,natoms))
    for i in range(0,nframes):
        c_array[i]=disp[i]*(arrh_term[i]-arrh_term_av)
    c_array=np.mean(c_array,axis=0)

    #calculates the normalization factor (i.e. sqrt(<dispi**2><arrh_termi**2>))
    arrh_term2=np.square(arrh_term)
    arrh_term2_av=np.average(arrh_term2)
    disp2=np.square(disp)
    disp2_av=np.mean(disp2,axis=0)
    norm_factor=np.sqrt(disp2_av*arrh_term2_av)

    #normalizes c and reshape it in a vector-like shape
    c_array=c_array/norm_factor #carray esta bien
    c_vector=np.reshape(c_array, natoms)
    c_vector_norm=np.sqrt(np.sum(np.square(c_vector)))
    c_vector=c_vector/c_vector_norm
    np.savetxt('coordinate_mull.dat',c_vector,fmt='%10.5f')

    write_pdb_with_charges('test.pdb',topfile,natoms,names,resids,c_vector)

    print('-------------------------------------------')
    print('The following files have been written:')
    print('energy_diff_raw.dat')
    print('energy_diff_smoothened.dat')
    print('energy_diff_derivative.dat')
    print('arrh_term.dat')
    print('P_enerdiff.dat')
    print('coordinate_mull.dat')
    print('-------------------------------------------')


def KRR_coord(rfile,topfile,e0file,e1file,filter,framespertraj,\
    temperature,alpha,minweight,dofiltered,testsize,kernelalpha,kernelgama,search_best_params):
    if (dofiltered):
        print('Reading coordinates')
        natoms, nframes = get_natoms_nframes(rfile)
        print("natoms:", natoms)
        print("nframes:", nframes)
        print("Filtering...")
        r_filtered,energy_diff_filtered,remainingframes= \
        read_and_filter(rfile,topfile,e0file,e1file,filter,framespertraj,temperature,alpha,minweight,natoms,nframes)
        efectiveframes = len(remainingframes)
        print("efectiveframes:", efectiveframes)
    else:
        print("Reading filtered files")
        natoms, efectiveframes = get_natoms_nframes('qm_filtered.xyz')
        print("natoms:", natoms)
        print("efectiveframes:", efectiveframes)
        energy_diff_filtered=np.loadtxt('energy_diff_filtered.dat',usecols=[1])
        remainingframes=np.loadtxt('energy_diff_filtered.dat',usecols=[0])
        remainingframes=remainingframes.astype(int)
        r_filtered = get_r_md(topfile,'qm_filtered.xyz',natoms,efectiveframes)

    print('Calculating averages')
    r_av = np.array([np.mean(r_filtered,axis=0)])
    # r = np.append(r,r_av,axis=0) #average coordinates could be appended as last frame

    print('Evaluating displacement vector (r - <r>)')
    disp = r_filtered[:]- r_av
    disp_vector=np.zeros((efectiveframes,3*natoms))

    for i in range(0,efectiveframes):
        jj=0
        for j in range(0,natoms):
            disp_vector[i][jj]=disp[i][j][0]
            disp_vector[i][jj+1]=disp[i][j][1]
            disp_vector[i][jj+2]=disp[i][j][2]
            jj=jj+3


    # here disp_vector is X and energy_diff_raw is Y
    print('Fitting KRR (rbf)')
    print('Test size is:', testsize)
    np.savetxt("disp_vector.dat",disp_vector,fmt='%10.5f')
    X_train, X_test, y_train, y_test = train_test_split(disp_vector, \
    energy_diff_filtered, test_size=testsize, random_state=42)
    if(search_best_params):
        print('Searching best parameters...')
        krr = GridSearchCV(
        KernelRidge(kernel="rbf", gamma=0.1),
        param_grid={"alpha": [1e0, 0.1, 1e-2, 1e-3], "gamma": np.logspace(-2, 2, 5)},
        )
        krr.fit(X_train, y_train)
        print(f"Best KRR with params: {krr.best_params_} and R2 score: {krr.best_score_:.3f}")
    else:
        print('Using specified parameters alpha=',kernelalpha,'and gamma=', kernelgama)
        krr = KernelRidge(alpha=kernelalpha,kernel='rbf',gamma=kernelgama)
        krr.fit(X_train, y_train)
    y_hat=krr.predict(X_test)
    prediction = np.column_stack((y_test,y_hat))
    np.savetxt("prediction.dat",prediction,fmt='%10.5f')
