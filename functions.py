import numpy as np
import MDAnalysis as md
from sklearn.metrics.pairwise import euclidean_distances
from MDAnalysis.coordinates.memory import MemoryReader

def get_atomic_numbers(names):
    dicZ=[("H",1),("He",2),("Li",3),("Be",4),("B",5),("C",6),("N",7),("O",8),("F",9),("Ne",10),("P",15),("S",16),("Fe",26)]
    atomic_numbers=[]
    for i in range(0,len(names)):
      for j in range(0,len(dicZ)):
        if (dicZ[j][0]==names[i]):
          atomic_numbers.append(dicZ[j][1])
    return atomic_numbers

def get_natoms_nframes(rfile):
    temp=open(rfile)
    temp=temp.readlines()
    nlines=len(temp)
    natoms=int(temp[0])
    nframes=int(nlines/(natoms+2))
    return natoms, nframes

def get_r_md(topfile,rfile,natoms,nframes): #This is probably not the best way to use this library
    u = md.Universe(topfile,rfile,dt=0.0005) #arbitrary dt value for us
    r=np.zeros((nframes,natoms,3))
    for i in range(0,nframes):
      for j in range(0,natoms):
          r[i][j][0]=u.trajectory[i][j][0]
          r[i][j][1]=u.trajectory[i][j][1]
          r[i][j][2]=u.trajectory[i][j][2]
    return r

def get_r(rfile,natoms,nframes):
    temp=open(rfile)
    temp=temp.readlines()
    r=np.zeros((nframes,natoms,3))
    for i in range(0,nframes):
      for j in range(0,natoms):
          x,y,z=temp[i*(natoms+2)+j+2].split()[1:4]
          r[i][j][0]=float(x)
          r[i][j][1]=float(y)
          r[i][j][2]=float(z)
    return r

def get_names_and_resids(topfile,natoms):
    u = md.Universe(topfile,dt=0.0005) #arbitrary dt value for us
    names=[]
    resids=[]
    for i in range(0,natoms):
      name=u.atoms[i].name
      resid=u.atoms[i].resid
      names.append(name)
      resids.append(resid)
    return names, resids

def get_distance_matrix(r):
  dmat=[]
  for i in (range(0,len(r))):
    m=euclidean_distances(r[i],r[i])
    dmat.append(m)
  return np.array(dmat)

def flip(energy,framespertraj):
  fliped_energy=energy
  doflip=False
  for i in range(0,len(energy)):
    temp=energy[i][1]-energy[i][0]
    if (temp < 0.000):
      doflip=True
      temp=-temp
    if((i % framespertraj == 0) and (i != 0)):
      doflip=False
    if(doflip):
      temp0=energy[i][0]
      temp1=energy[i][1]
      fliped_energy[i][0]=temp1
      fliped_energy[i][1]=temp0
  return fliped_energy

def get_derivative(energydiffsmooth):
    enerdifderiv=np.zeros(len(energydiffsmooth))
    for i in range(0,len(energydiffsmooth)-1):
      temp=energydiffsmooth[i]-energydiffsmooth[i+1]
      if(i>0):
        temp1=energydiffsmooth[i-1]-energydiffsmooth[i]
      if(i==0):
        temp1=energydiffsmooth[i+1]-energydiffsmooth[i+2]
      temp2=(temp-temp1)
      if((abs(temp2)>0.002)):
        enerdifderiv[i] = 0
      else:
        enerdifderiv[i]=temp
      if(i==(len(energydiffsmooth)-2)):
        enerdifderiv[i+1]=enerdifderiv[i]
    return enerdifderiv

def get_centrality(c_vector,natoms):
    centrality=np.zeros((natoms,2))
    ii=0
    for i in range(0,natoms):
      centrality[i][0]=i+1
      centrality[i][1]=np.sqrt(c_vector[ii]**2+c_vector[ii+1]**2+c_vector[ii+2]**2)
      ii=ii+3
    return centrality

def write_nmd(fname,names,natoms,resids,r,c_vector):
    title='MOL'
    f = open(fname, "w")
    f.write(title)
    f.write("\n")
    f.write('names '+' '.join(names))
    f.write("\n")
    f.write('resids ')
    for i in range(0,natoms):
      f.write(str(resids[i])+' ')
    f.write("\n")
    f.write('coordinates')
    for i in range(0,natoms):
      f.write(' '+f'{r[0][i][0]:.3f}'+' '+f'{r[0][i][1]:.3f}'+' '+f'{r[0][i][2]:.3f}')
    f.write("\n")
    f.write('mode 1')
    for i in range(0,3*natoms):
      f.write(' '+f'{c_vector[i]:.3f}')
    f.close()

def write_animation(fname,atomic_numbers,rav,c_vector,natoms,nframes_animation,amplitude):
    f = open(fname, "w")
    rav_vector=np.reshape(rav, 3*natoms)
    oneovernsteps=1/nframes_animation
    for i in range(0,nframes_animation):
        f.write('           '+str(natoms))
        f.write('\n')
        f.write('\n')
        ii=0
        for j in range(0,natoms):
            z=atomic_numbers[j]
            re1=rav_vector[ii]+c_vector[ii]*amplitude*np.sin(i*2*np.pi*oneovernsteps)
            re2=rav_vector[ii+1]+c_vector[ii+1]*amplitude*np.sin(i*2*np.pi*oneovernsteps)
            re3=rav_vector[ii+2]+c_vector[ii+2]*amplitude*np.sin(i*2*np.pi*oneovernsteps)
            f.write(' '+str(z)+'   '+ f'{re1:.8f}' + '       ' + f'{re2:.8f}' + '       ' + f'{re3:.8f}')
            f.write('\n')
            ii=ii+3


def write_projection(fname,atomic_numbers,r_av,disp,c_vector,refindex,natoms,nframes,wt):
    f = open(fname, "w")
    r_av_vector=np.reshape(r_av, 3*natoms)
    dispref_vector=np.reshape(disp[refindex], 3*natoms)
    for i in range(0,nframes):
        f.write('           '+str(natoms))
        f.write('\n')
        f.write('\n')
        ii=0
        for j in range(0,natoms):
            z=atomic_numbers[j]
            re1=r_av_vector[ii]+dispref_vector[ii]+c_vector[ii]*wt[i]
            re2=r_av_vector[ii+1]+dispref_vector[ii+1]+c_vector[ii+1]*wt[i]
            re3=r_av_vector[ii+2]+dispref_vector[ii+2]+c_vector[ii+2]*wt[i]
            f.write(' '+str(z)+'   '+ f'{re1:.8f}' + '       ' + f'{re2:.8f}' + '       ' + f'{re3:.8f}')
            f.write('\n')
            ii=ii+3

def get_natoms_nframes_mulliken(mfile):
    temp=open(mfile)
    sum=False
    natoms=0
    for line in temp:
        if 'Total' in line:
            break
        if(sum):
            natoms=natoms+1
        if 'Atom' in line:
            sum=True
    temp.seek(0)
    temp=temp.readlines()
    nlines=len(temp)
    nframes=int(nlines/(natoms+5))
    return natoms, nframes

def get_mull(mfile,natoms,nframes):
    temp=open(mfile)
    temp=temp.readlines()
    q=np.zeros((nframes,natoms))
    for j in range(0,natoms):
        charge=temp[j+3].split()[2]
        q[0][j]=float(charge)
    start=3+natoms+2
    for i in range(0,nframes-1):
      for j in range(0,natoms):
          charge=temp[start+i*(natoms+5)+j+3].split()[2]
          q[i+1][j]=float(charge)
    return q

def write_pdb_with_charges(fname,topfile,natoms,names,resids,q):
    u = md.Universe(topfile,dt=0.0005) #arbitrary dt value for us
    r=u.trajectory[0]
    f = open(fname, "w")
    f.write('CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1           1')
    f.write('\n')
    for i in range(0,natoms):
        newline='ATOM      '+str(i)+'  '+str(names[i])+'       X   '+str(resids[i])+'       '
        newline=newline+f'{r[i][0]:.3f}' + '  ' + f'{r[i][1]:.3f}' + '  ' + f'{r[i][1]:.3f}'+'  '
        newline=newline+'  0.00  '+ f'{q[i]:.2f}'+'           '+str(names[i])
        f.write(newline)
        f.write('\n')

def get_P_ener_diffs(energy_diff_smoothened,alpha,kt,nframes,window_size):
    arrh_term_nosign=np.exp(-abs((energy_diff_smoothened)/(alpha*kt)))
    window_size=50
    nwindows=int(nframes/window_size)
    P_enerdiff_by_windows=np.zeros((nwindows,2))
    Q_enerdiff_by_windows=np.zeros(nwindows)
    ediff_times_arr_term_nosing=energy_diff_smoothened*arrh_term_nosign
    for i in range(0,nwindows):
        start=i*window_size
        end=start+window_size
        Q_enerdiff_by_windows[i]=np.sum(arrh_term_nosign[start:end])
        P_enerdiff_by_windows[i][1]=np.sum(ediff_times_arr_term_nosing[start:end])/Q_enerdiff_by_windows[i]
        P_enerdiff_by_windows[i][0]=(start+end-1)/2
    # Q_enerdiff=np.sum(arrh_term_nosign)
    Q_enerdiff=1
    P_enerdiff=arrh_term_nosign/Q_enerdiff

    return P_enerdiff, P_enerdiff_by_windows



def get_data_filtered_by_weight(r,energy_diff_raw,nframes,minweight,temperature,alpha):
    kB=0.000003173 #kB in hartree
    kt=kB*temperature
    PE, temp = get_P_ener_diffs(energy_diff_raw,alpha,kt,nframes,50)
    PE=PE/(np.max(PE))
    np.savetxt('PE.dat',PE)
    framestodelete=[]
    remainingframes=[]
    for j in range(0,nframes):
        if PE[j]<minweight:
            framestodelete.append(j)
        else:
            remainingframes.append(j)
    # efectiveframes=nframes-len(framestodelete)
    r_filtered=np.delete(r,framestodelete,0)
    energy_diff_filtered=np.delete(energy_diff_raw,framestodelete)
    return r_filtered, energy_diff_filtered, remainingframes


def write_xyz_of_selectedframes(topfile,r_filtered,natoms,efectiveframes):
    u = md.Universe(topfile,r_filtered,in_memory=True,dt=0.0005) #arbitrary dt value for us
    p = u.select_atoms("all")
    with md.Writer("qm_filtered.xyz", p.n_atoms) as W:
      for ts in u.trajectory:
        W.write(p)

def read_and_filter(rfile,topfile,e0file,e1file,filter,framespertraj,temperature,alpha,minweight,natoms,nframes):
    print('Loading trajectory')
    r = get_r_md(topfile,rfile,natoms,nframes)
    print('Loading energies evolution')
    e0 = np.loadtxt(e0file)
    e1 = np.loadtxt(e1file)
    print('Processing')
    energy = np.column_stack((e0,e1))
    energy = flip(energy,framespertraj)
    energy_diff_raw = energy[:,1]-energy[:,0]
    r_filtered,energy_diff_filtered,remainingframes= \
    get_data_filtered_by_weight(r,energy_diff_raw,nframes,minweight,temperature,alpha)
    print(np.shape(r_filtered), np.shape(energy_diff_filtered))
    efectiveframes=len(energy_diff_filtered)
    print("nframes used for analysis:", efectiveframes, " (",efectiveframes*100/nframes," %)")
    frames = [ j for j in range(nframes) ]
    frames = np.array(frames)
    np.savetxt("energy_diff_raw_withindex.dat",np.column_stack((frames,energy_diff_raw)),fmt='%10.5f')
    np.savetxt("energy_diff_filtered.dat",np.column_stack((remainingframes,energy_diff_filtered)),fmt='%10.5f')
    write_xyz_of_selectedframes(topfile,r_filtered,natoms,efectiveframes)
    return r_filtered,energy_diff_filtered,remainingframes
