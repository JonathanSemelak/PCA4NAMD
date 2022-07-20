import numpy as np
import matplotlib.pyplot as plt
dt = 0.0005
pstofs = 1000
dt = dt * pstofs
file='SMA.out'
start=1
end=100
# exceptions=[3,25,31,42]
exceptions=[]
reasons=[]
forceminlen = True
minlenforced = 16000
printall = False
printaverage = True
minminlen = 100
for i in range(start,end+1):
    isexception=False
    l=0
    path=str(i)+'/'+file
    temp=open(path)
    for line in temp:
      if line.startswith(' poblacion1'):
        l=l+1
    temp.seek(0)
    if (l<minminlen):
        exceptions.append(i)
        reasons.append('Trajectory is too short')
        isexception=True
    if (not isexception):
      for line in temp:
        if 'No convergence' in line:
          exceptions.append(i)
          reasons.append('No convergence')
          break
      temp.seek(0)

print("Total exceptions:", len(exceptions))
for i in range(0,len(exceptions)):
    print(exceptions[i],reasons[i])

if (len(exceptions)>0 and exceptions[0]<= start):
  for i in range(start,end+1):
    newstart=start
    if (i not in exceptions):
      newstart=i
      break
  start=newstart
for i in range(start,end+1):
    l=0
    if (i not in exceptions):
      populationS0, populationS1 = [], []
      path=str(i)+'/'+file
      temp=open(path)
      for line in temp:
        if line.startswith(' poblacion1'):
          l=l+1
      temp.seek(0)
      print("traj ", i, " len is: ", l)
      if(i==start):
        minlen=l
        minlenindex=i
      else:
        if (l < minlen):
          minlen=l
          minlenindex=i

if (forceminlen == True):
  if (minlen<minlenforced):
      print("forced minlen is smaller than real minlen")
      print("not forcing")
      print("min len is ", minlen, "(traj: ", minlenindex, ")")
  else:
      minlen=minlenforced
      print("min len is ", minlen, "(forced)")


ndec=0
E0, E1 = [], []
for i in range(start,end+1):
    if (i not in exceptions):
      populationS0, populationS1 = [], []
      path=str(i)+'/'+file
      temp=open(path)
      count=0
      for line in temp:
        if(count<minlen):
          if 'Final energy' in line:
            e0=float(line.split()[7])
            E0.append(e0)
          if line.startswith(' STATE 1      ENERGY='):
            e1=float(line.split()[3])
            E1.append(e0+e1)
          if line.startswith(' poblacion1'):
            value=float(line.split()[1])
            populationS0.append(value)
          if line.startswith(' poblacion2'):
            value=float(line.split()[1])
            populationS1.append(value)
            count=count+1 #because this appears last
      populationS0=np.array(populationS0)
      populationS1=np.array(populationS1)
      if(populationS1[minlen-1]<0.5):
          ndec=ndec+1
      if(i==start):
        populationS0_av=populationS0
        populationS1_av=populationS1
      else:
        populationS0_av=populationS0_av+populationS0
        populationS1_av=populationS1_av+populationS1
      if (printall):
        x = [ j for j in range(minlen) ]
        x = np.array(x)*dt
        plt.plot(x,populationS0, label=str(i))
        plt.plot(x,populationS1, label=str(i))
E0=np.array(E0)
E1=np.array(E1)
print("decayed: ", ndec*100/(end-len(exceptions)), "%")

if (printaverage):
  x = [ j for j in range(minlen) ]
  x = np.array(x)*dt
  populationS0_av=populationS0_av/(end-len(exceptions))
  populationS1_av=populationS1_av/(end-len(exceptions))
  plt.plot(x,populationS0_av, label='average S0')
  plt.plot(x,populationS1_av, label='average S1')

if (printaverage or printall):
  plt.legend()
  plt.xlabel('time (fs)')
  plt.ylabel('Population S1')
  plt.show()
  plt.savefig('S0S1.png')

dataname='test.dat'
np.savetxt("E0.dat",E0,fmt='%10.5f')
np.savetxt("E1.dat",E1,fmt='%10.5f')



# data = np.loadtxt(dataname, usecols=[0,1], skiprows =1)
