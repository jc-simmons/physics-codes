#!/usr/bin/env python3

import numpy as np
from scipy.sparse.linalg import eigsh
from numpy import kron, identity
import timeit
#from mpi4py import MPI

#========================================================

# Density Matrix Renormalization Group (DMRG) alorithms applied to a Heisenberg Spin Chain
# 
# Project for High-Performance Computational Physics course
#
# All necessary functions in single file

#----setup--------------------------------

#iflag is eigenvalue method used to find superblock groundstate
# 0 lanczos,
# 1 lanczos using MPI 
# 2 numpy.eigh
# 3 scipy.sparse.eigsh
iflag=0


chain_len=22     # target size of 1D chain system
trunc_limit=15   # truncation size, eigenvectors of reduced density operator kept
sys_size=1       # system length initialization
sys_dim=2        # system dimension initialization


#Sz and Sp spin operators for single site
sz=np.array([[0.5,0],[0,-0.5]], dtype='d')
sp=np.array([[0, 1], [0, 0]], dtype='d')

site={}
site["Sp"]=sp
site["Sz"]=sz

#define starting system and environment Hamiltonian
hamiltonian=np.zeros((2,2), dtype='d')

sys_block={"H" : hamiltonian}
sys_block["Sz"] = sz
sys_block["Sp"] = sp 

env_block={"H" : hamiltonian}
env_block["Sz"] = sz
env_block["Sp"] = sp 


#this is in case of using parallel Lanczos, requires same starting vector on each node
np.random.seed(0)


if iflag==1:
    
   comm=MPI.COMM_WORLD
   rank=comm.Get_rank()
   size=comm.Get_size()


#-----------------------------------------




#Hamtiltonian interaction between block1 and block2
def interaction(block1,block2):
    
    interaction= (1/2)*(kron(block1["Sp"],block2["Sp"].conj().T)+kron(block1["Sp"].conj().T,block2["Sp"]))\
    +kron(block1["Sz"],block2["Sz"])
       
    return interaction


#tacks on to current block
def enlarge(block):
    
    n=np.shape(block["H"])[0]
    
    enlarged_block={}
    enlarged_block["H"] = kron(block["H"],identity(2)) + interaction(block,site)
    enlarged_block["Sp"]= kron(identity(n),sp)
    enlarged_block["Sz"]= kron(identity(n),sz)

    return enlarged_block


# Attempt to create a more efficient Kronecker product using MPI
# communication overhead seems to hinder efficinecy at least for the matrix sizes tested
def my_kron(A,B):


        n_A=np.shape(A)[0]
        n_B=np.shape(B)[0]
        prod=np.zeros((n_A*n_B,n_A*n_B))

        final=np.zeros((n_A*n_B,n_A*n_B))
   
        irank=rank

        while irank<n_A*n_B:

            r_i = irank // n_A
            c_i = irank % n_A

            prod[r_i*n_B:r_i*n_B+n_B,c_i*n_B:c_i*n_B+n_B] = A[r_i,c_i]*B

            irank=irank+size


        comm.Barrier()

        comm.Reduce(prod,final,op=MPI.SUM,root=0)



        return final



# MPI and non-MPI versions of Lanzcos iteration
# Keeps iterating until groundstate energy is converged to within 'tol'
# communication seems to hinder efficiency, step is highlighted below
def lanczos(A):

   if iflag==0:

     tol=10**(-8)
     en_prev=np.zeros(1)
     en_next=np.ones(1)

     n=A.shape[1]

     V=np.zeros((n,n+1))

     v_start= np.random.rand(n)
     v_start=v_start/np.linalg.norm(v_start)
     V[:,0]=v_start.copy()

     V_i=V[:,0]*1


     z=np.zeros(n)
     T=np.zeros((n,n))
     alpha=np.zeros(n+1)
     beta=np.zeros(n+1)
     beta[0]=0
     super_groundstate=np.zeros(n)

     i=0
     en_low=np.zeros(1)

     while (abs(en_next-en_prev)>tol or abs(en_low-en_next)>tol):

        en_prev=en_next*1

        z=np.dot(A,V_i)

        alpha[i]=np.dot(V[:,i],z)

        z=z-alpha[i]*V[:,i]-beta[i]*V[:,i-1]

        beta[i+1]=np.linalg.norm(z)
        V[:,i+1]=z/beta[i+1]

        T[i,i]=alpha[i]
        T[i-1,i]=beta[i]
        T[i,i-1]=beta[i]

        V_i=V[:,i+1]*1

        vals,vecs =np.linalg.eig(T[:i+1,:i+1])
        en_next=vals[0]*1
        T_groundstate=vecs[:,0]*1
        super_groundstate=np.dot(V[:,:i+1],T_groundstate)


        if en_low>en_next:
          en_low=en_next*1

        i+=1




   elif iflag==1:

     tol=10**(-8)
     en_prev=np.zeros(1)
     en_next=np.ones(1)

     n=A.shape[1]

     segs=int(n/size)
     chunks=np.ones(size,dtype='int')*segs
     displ=np.zeros(size,dtype='int')

     for i in range(size):
      displ[i]=i*segs


     z=np.zeros(n)
     T=np.zeros((n,n))
     alpha=np.zeros(n+1)
     beta=np.zeros(n+1)
     beta[0]=0
     z0=np.zeros((segs))
     super_groundstate=np.zeros(n)


     V=np.zeros((n,n+1))
     v_start=np.random.rand(n)
     v_start=v_start/np.linalg.norm(v_start)
     V[:,0]=v_start.copy()


     V_i=V[:,0]*1

     i=0
     en_low=np.zeros(1)

     while (abs(en_next-en_prev)>tol or abs(en_low-en_next)>tol):

         en_prev=en_next*1
 
         z0=np.dot(A[displ[rank]:displ[rank]+segs,:],V_i)
         
         #need z array, this affects parallel efficiency
         comm.Gatherv(z0[:],[z,chunks,displ,MPI.DOUBLE],root=0)
         comm.Bcast(z,root=0)         

         alpha[i]=np.dot(V[:,i],z)
 
         z=z-alpha[i]*V[:,i]-beta[i]*V[:,i-1]

         beta[i+1]=np.linalg.norm(z)
         V[:,i+1]=z/beta[i+1]

         T[i,i]=alpha[i]
         T[i-1,i]=beta[i]
         T[i,i-1]=beta[i]

         V_i=V[:,i+1]*1

         vals,vecs =np.linalg.eig(T[:i+1,:i+1])
         en_next=vals[0]*1
         T_groundstate=vecs[:,0]*1
         super_groundstate=np.dot(V[:,:i+1],T_groundstate)

         if en_low>en_next:
           en_low=en_next*1

         i+=1
         comm.Barrier()




   return en_next, super_groundstate



#performs a single DMRG step, either Infinite or Finite algorithms 
def dmrg_step(sys_block,env_block,sys_size):

    
    sys_dim=np.shape((sys_block["H"]))[0]
    env_dim=np.shape((env_block["H"]))[0]
    
    #grow each by a site
    sys_block= enlarge(sys_block)
    env_block = enlarge(env_block)

    
    sys_size+=1

    
    sys_dim=sys_dim*2
    env_dim=env_dim*2

    # Hamiltonians of system/environment in superblock basis + interaciton between
    superH = kron(sys_block['H'],identity(env_dim))+kron(identity(sys_dim),env_block['H'])\
             + interaction(sys_block,env_block)
 



    # iflag --> eigenvalue method selection
    
    # Lanczos or Lanczos with parallel matrix-vector products
    if (iflag==0 or iflag==1):
      energy,super_groundstate = lanczos(superH)

    # basic Python eigenvalue routine
    if iflag==2:
    
      energies, super_states = np.linalg.eigh(superH)
      energy=energies[0]
      super_groundstate=super_states[:,0]    

    # very efficient due to ARPACK backend
    if iflag==3:
   
      energy, super_groundstate = eigsh(superH, k=1, which="SA")


    # contructs reduced density operator 
    groundstate_matrix =np.reshape(super_groundstate, (sys_dim,env_dim) ,order='C')

    reduced_density = np.dot(groundstate_matrix, groundstate_matrix.conj().T)
    
    evals, evecs = np.linalg.eigh(reduced_density)
    
    
    # take the eigenvectors of reduced density operator associated with largest trunc_limit eigenvectors    
    trunc_states= evecs[:,-trunc_limit:]
    
    trunc_states=np.flip(trunc_states,axis=1)
    
    transfer=trunc_states
    


    # transforms all enlarged system block operators to new basis, independent of the type of lattice
    for key, val in sys_block.items():
        
        sys_block[key]= transfer.conj().T @ val @ transfer

    
    # if the size of the englarged block is above the maximum, bring it back down
    # allows m to grow to trunc_limit where it is maintained
    if sys_dim > trunc_limit:
        sys_dim = trunc_limit
    
        
    return  sys_block , energy, sys_size, sys_dim



# this is to calculate a full matrix in the desired basis (up/down spin for Heisenberg chain of length L)
# repeatedly enlarges by a site 
def exact_ground(L):
    
    
    Hfull=dict(sys_block)
    
    for i in range(L-1):
    
        Hfull=enlarge(Hfull)
 
    vals, vecs = eigsh(Hfull["H"], k=1 , which="SA")
    
    print('Groundstate energy: ' , vals[0]/L)
    
    return 



# grows an input block up to size chain_len
# sys_dim grows to trunc_limit
# sys_size grows to size of the system block lattice (not the Hilbert space)
def infinite(block,sys_size):
      

    while 2*sys_size<chain_len:
        #environment is a mirror of the system
        block1=dict(block)
        block2=dict(block)
        
        [block,energy,sys_size, sys_dim] = dmrg_step(block1,block2,sys_size)
  
    
   
    
    if iflag==1:

      if rank==0:

         print('Groundstate energy: ' ,energy/chain_len)

    else:

      print('Groundstate energy: ' ,energy/chain_len)


    return #sys_size, sys_dim  

    


def finite(block,sys_size):
    
    block_list = {}
   
    block_list[0,sys_size]=block   #system blocks
    block_list[1,sys_size]=block   #environment blocks
 
    #growth according to infinite algorithm
    while 2*sys_size<chain_len:
        
        block1=dict(block)
        block2=dict(block)
        
        [block,energy,sys_size, sys_dim] = dmrg_step(block1,block2,sys_size)

        #record
        block_list[0,sys_size]=dict(block)
        block_list[1,sys_size]=dict(block)

        
        
    nsweep=1
    
    #use as system and environment index    
    i=0
    j=1
    
    
    tol=1.e-8 #
    epsilon_p=0.
    epsilon_q=energy/(sys_size*2)
    
    # run until converged to within tol
    while (abs(epsilon_q-epsilon_p)>tol):
        
        if chain_len-sys_size-2>1:

            epsilon_p=epsilon_q*1.
            
            enblock=dict(block_list[j,chain_len-sys_size-2])
            [block,energy,sys_size, sys_dim]=dmrg_step(block,enblock,sys_size)  
            
            #keep recording only the growing block
            block_list[i,sys_size]=dict(block)
            
            epsilon_q = energy/chain_len
        
        
        else:
            
            nsweep+=1
            print('---start sweep--- ', nsweep)
            sys_size=1
            # swap the system/environment blocks
            block=dict(block_list[j,sys_size])
            
            # swaps values of system and environment
            h=i*1.
            i=j*1.
            j=h*1.


    if iflag==1:
      # only print info once for 0th node if parallel
      if rank==0:
    
         print('Groundstate energy: ' ,energy/chain_len)
    
    else: 
     
      print('Groundstate energy: ' ,energy/chain_len)    


    return 




start=timeit.default_timer()

infinite(sys_block,sys_size)
#finite(sys_block,sys_size)

stop=timeit.default_timer()


# only print timing info once for 0th node if parallel
if iflag==1:

    if rank==0:

  
       print('time: ', stop-start)
else:

    print('time: ', stop-start)
