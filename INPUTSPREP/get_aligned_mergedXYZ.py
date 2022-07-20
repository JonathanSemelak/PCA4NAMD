import numpy as np
import MDAnalysis as md
from MDAnalysis.analysis import align

file='qm.xyz'
filemull='mulliken'
fileout='SMA.out'
mergemulliken=True
start=1
end=100
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
pdbname = 'qmpdb.pdb'
dt=0.0005
minminlen = 100

for i in range(start,end+1):
    isexception=False
    l=0
    path=str(i)+'/'+fileout
    temp=open(path)
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
