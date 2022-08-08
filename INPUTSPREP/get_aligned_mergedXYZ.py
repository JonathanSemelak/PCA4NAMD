import numpy as np
import MDAnalysis as md
from MDAnalysis.analysis import align
import sys

file='qm.xyz'
filemull='mulliken'
fileout='SMA.out'
mergemulliken=True
start=1
end=249
exceptions=[]
minlen = 16000
minminlen = 100
natoms = 9
outputname = 'merged_'+file
outputnamemull = 'merged_'+filemull
maxlinestoread = minlen*(natoms+2)
maxlinestoreadmull = minlen*(natoms+5)
readall = True
doalign = True
usepdb = True
pdbname = 'qm.pdb'
dt=0.0005
minminlen = 100

for i in range(start,end+1):
    isexception=False
    l=0
    path=str(i)+'/'+fileout
    temp=open(path)
    print('checking convergence on traj ', i)
    for line in temp:
      if line.startswith(' poblacion1'):
        l=l+1
    temp.seek(0)
    if (l<minminlen):
        exceptions.append(i)
        isexception=True
    if (not isexception):
      for line in temp:
        if 'No convergence' in line:
          exceptions.append(i)
          break
      temp.seek(0)

for i in range(start,end+1):
    if (i not in exceptions):
      path=str(i)+'/'+file
      temp=open(path)
      foundE0=0
      foundE1=0
      foundp1=0
      foundp2=0
      count=0
      for line in temp:
        if(count<minlen):
          if 'Final energy' in line:
            foundE0=foundE0+1
          if line.startswith(' STATE 1      ENERGY='):
            foundE1=foundE1+1
          if line.startswith(' poblacion1'):
            foundp1=foundp1+1
          if line.startswith(' poblacion2'):
            foundp2=foundp2+1
            count=count+1 #because this appears last
          if (not ((foundE0 == foundE1) and (foundE1 == foundp1) and (foundp1 == foundp2))):
            print('Something went wrong when reading traj ',i)
            print('An energy value or population is missing.')
            print('Stopping program...')
            sys.exit()
      temp.seek(0)

print("Total exceptions:", len(exceptions))

if (len(exceptions)>0 and exceptions[0]<= start):
  for i in range(start,end+1):
    newstart=start
    if (i not in exceptions):
      newstart=i
      break
  start=newstart

#If readall == True, a trajectory is generated assembling all qm.xyz files
if (readall):
  output=open(outputname, 'w')
  for i in range(start,end+1):
      if (i not in exceptions):
        path=str(i)+'/'+file
        temp=open(path)
        count=0
        for line in temp:
          if(count<maxlinestoread):
            output.write(line)
            count=count+1
output.close()
#If doalign == True, the trajectory is  aligned to the first frame coordinates
if (doalign):
  outputnamealigned='aligned_'+outputname
  u = md.Universe(pdbname,outputname,dt=dt)
  ref = u.select_atoms('resid 1')
  alignment = align.AlignTraj(u, ref, filename=outputnamealigned)
  alignment.run()

if (mergemulliken):
  outputmull=open(outputnamemull, 'w')
  for i in range(start,end+1):
    if (i not in exceptions):
      path=str(i)+'/'+filemull
      temp=open(path)
      count=0
      for line in temp:
        if(count<maxlinestoreadmull):
          outputmull.write(line)
          count=count+1
  outputmull.close()
