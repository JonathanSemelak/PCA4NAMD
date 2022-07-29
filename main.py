from functions import*
from correlations import*

Analysis_type='Coordinates-PCA-like'
mfile='mulliken'
rfile='qm.xyz'
topfile='qm.pdb'
e0file='E0.dat'
e1file='E1.dat'
temperature=300
alpha=100
framespertraj=16000
filter=99
nframes_animation=100
amplitude=5

print('Analysis_type: ', Analysis_type)
if (Analysis_type=='Coordinates-PCA-like'):
    covar_coord_arrh(rfile,topfile,e0file,e1file,temperature,alpha,framespertraj,filter,nframes_animation,amplitude)
elif (Analysis_type=='Mulliken-PCA-like'):
    covar_mulliken_arrh(mfile,topfile,e0file,e1file,temperature,alpha,framespertraj,filter)
elif (Analysis_type=='Coordinates-KRR'):
    KRR_coord(rfile,topfile,e0file,e1file,filter,framespertraj)
#     print('PCA_type wrong value')
# elif (Analysis_type=='Mulliken-KRR'):
#     print('PCA_type wrong value')
else:
    print('Analysis_type wrong value')
    print('Available options are: Coordinates-PCA-like, Mulliken-PCA-like, Coordinates-KRR, Mulliken-KRR')
