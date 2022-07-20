from functions import*
from PCAs import*

PCA_type='Coordinates'
rfile='qm.xyz'
topfile='qm.pdb'
e0file='E0.dat'
e1file='E1.dat'
temperature=100
alpha=100
framespertraj=16000
filter=99
nframes_animation=100
amplitude=5

print('PCA_type: ', PCA_type)
if (PCA_type=='Coordinates'):
    correlation_coord(rfile,topfile,e0file,e1file,temperature,alpha,framespertraj,filter,nframes_animation,amplitude)
elif (PCA_type=='Mulliken'):
    print('PCA_type wrong value')
elif (PCA_type=='Kernel_Coordinates'):
    print('PCA_type wrong value')
elif (PCA_type=='Kernel_Mulliken'):
    print('PCA_type wrong value')
else:
    print('PCA_type wrong value')
    print('Available options are: Coordinates, Mulliken, Kernel_Coordinates, Kernel_Mulliken')
