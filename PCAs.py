import numpy as np
import MDAnalysis as md
from sklearn.metrics.pairwise import euclidean_distances
from scipy.signal import savgol_filter
from functions import*

def correlation_coord(rfile,topfile,e0file,e1file,temperature,alpha,framespertraj,filter,nframes_animation,amplitude):
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
    energy_diff_smoothened = savgol_filter(energy[:,1]-energy[:,0], filter, 3)
    energy_diff_derivative = get_derivative(energy_diff_smoothened)

    np.savetxt("energy_diff_raw.dat",energy[:,1]-energy[:,0],fmt='%10.5f')
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


def correlation_mulliken(mfile,e0file,e1file,temperature,alpha,framespertraj,filter):
    kB=0.000003173 #kB in hartree
    kt=kB*temperature
    print('Reading coordinates')
    natoms, nframes = get_natoms_nframes_mulliken(mfile)
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
    energy_diff_smoothened = savgol_filter(energy[:,1]-energy[:,0], filter, 3)
    energy_diff_derivative = get_derivative(energy_diff_smoothened)
    np.savetxt("energy_diff_raw.dat",energy[:,1]-energy[:,0],fmt='%10.5f')
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
    arrh_term_nosign=np.exp(-abs((energy_diff_smoothened)/(alpha*kt)))
    window_size=50
    nwindows=int(nframes/window_size)
    P_enerdiff=np.zeros((nwindows,2))
    Q_enerdiff=np.zeros(nwindows)
    ediff_times_arr_term_nosing=energy_diff_smoothened*arrh_term
    for i in range(0,nwindows):
        start=i*window_size
        end=start+window_size
        Q_enerdiff[i]=np.sum(arrh_term_nosign[start:end])
        P_enerdiff[i][1]=np.average(ediff_times_arr_term_nosing[start:end])/Q_enerdiff[i]
        P_enerdiff[i][0]=(start+end-1)/2
    np.savetxt('P_enerdiff.dat',P_enerdiff,fmt='%10.5f')

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

    print('-------------------------------------------')
    print('The following files have been written:')
    print('energy_diff_raw.dat')
    print('energy_diff_smoothened.dat')
    print('energy_diff_derivative.dat')
    print('arrh_term.dat')
    print('P_enerdiff.dat')
    print('coordinate_mull.dat')
    print('-------------------------------------------')
