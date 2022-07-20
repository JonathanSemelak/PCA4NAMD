import numpy as np
import MDAnalysis as md
from sklearn.metrics.pairwise import euclidean_distances

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
