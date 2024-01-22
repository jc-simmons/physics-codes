import numpy as np
from scipy import linalg 
from scipy.sparse.linalg import eigs, LinearOperator
from scipy.special import hermite 
from itertools import combinations_with_replacement



def mvp(vector):
    

    result = np.zeros(len(vector))
    
    for i in range(d):
        
        interim = vector.copy()
        
        for j in range(d):
            
            interim = dim_matmul(interim, keo_mat[i,j,:,:], n, j)
           # interim = dim_matmul(interim, btest, n, j)
            
            
        result = result + interim


    #result = result + pot(vector)

    return result



def dim_matmul(vector, matrix, m, dim):
    
    u = np.reshape(vector, (int(m**(d-dim-1)), m, int(m**dim)) , order = 'F' )
    
    v = np.zeros(u.shape)
    
    
    for i in range(m):
        for j in range(m):
            
            v[:,i,:] = v[:,i,:] + u[:,j,:]*matrix[i,j]
            
    
    
    return np.ravel(v, order = 'F')
     


def pot(x):
    
    return x**2/2
            

def hob(n,x,d):

    result = 1/np.sqrt(2**n * np.factorial(n))*(1/np.pi)**(1/4)*np.exp(-x**2/2)*hermite(n)(x)
    
    return result




n=8
d=1


fcoefs = np.array([0.5, 0.5, 0.5])
fexp = np.array([2, 2, 2])


hob_positionOP = np.array([np.sqrt((i+1)/2) for i in range(n-1)])

dvr_pnts, dvr_trans = linalg.eigh_tridiagonal(np.zeros(n), hob_positionOP, eigvals_only=False)



keo_mat = np.zeros((d,d,n,n))  

for dim in range(d):
    for deriv in range(d):
        if dim==deriv:
            
            for i in range(n):
                for j in range(n):
                    
                    if i==j:
                        keo_mat[dim,deriv,i,j] = (i+1/2)*fcoefs[dim]
                    elif abs(i-j)==2:
                        keo_mat[dim,deriv,i,j] = -np.sqrt(max(i,j)*(max(i,j)-1))/2*fcoefs[dim]

        else:
            
            keo_mat[dim,deriv,:,:] = np.eye(n)
            
 
    
    
# replace derivative KEO matrices with DVR transformed counterpart            
#keo_mat[dim,dim,:,:] = [np.transpose(dvr_trans) @ keo_mat[dim,dim,:,:] @ dvr_trans for dim in range(d)]   

for i in range(d):
    for j in range(d):
        keo_mat[i,j,:,:] = np.transpose(dvr_trans) @ keo_mat[i,j,:,:] @ dvr_trans
        
        
ktest = keo_mat[0,0,:,:]


btest=np.zeros((n,n))
for i in range(n):
    btest[i,i] = 0.5 + i


A = LinearOperator( (n**d, n**d), matvec = mvp)
   
evals, evecs = eigs(A, k=4, which ='SM')