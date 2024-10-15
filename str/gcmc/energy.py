import numpy as np
import torch
from gpu_sel import *

def structure_eann(filename):
    # used for select a unoccupied GPU
    #gpu_sel()
    # gpu/cpu
    #device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    device = torch.device("cpu")
    # same as the atomtype in the file input_density
    atomtype=['Cu','Ce','O','C']
    #load the serilizable model
    pes=torch.jit.load("EANN_PES_DOUBLE.pt")
    # FLOAT: torch.float32; DOUBLE:torch.double for using float/double in inference
    pes.to(device).to(torch.double)
    # set the eval mode
    pes.eval()
    # save the lattic parameters
    cell=np.zeros((3,3),dtype=np.float)
    period_table=torch.tensor([1,1,1],dtype=torch.float,device=device)   # same as the pbc in the periodic boundary condition
    npoint=0
    rmse=0
    with open("./%s" %filename,'r') as f1:
        while True:
            string=f1.readline()
            if not string or npoint==1: break
            string=f1.readline()
            cell[0]=np.array(list(map(float,string.split())))
            string=f1.readline()
            cell[1]=np.array(list(map(float,string.split())))
            string=f1.readline()
            cell[2]=np.array(list(map(float,string.split())))
            string=f1.readline()
            species=[]
            cart=[]
            while True:
                string=f1.readline()
                if "abprop" in string: break
                tmp=string.split()
                tmp1=list(map(float,tmp[2:8]))
                cart.append(tmp1[0:3])
                species.append(atomtype.index(tmp[0]))
            #print(string.split())
            abene=float(string.split(':')[1])
    return period_table,cart,species,cell
    

def Energy_eann(period_table,cart,species,cell):
    # used for select a unoccupied GPU
    #gpu_sel()
    # gpu/cpu
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    # # same as the atomtype in the file input_density
    atomtype=['Cu','Ce','O','C']
    # #load the serilizable model
    pes=torch.jit.load("EANN_PES_DOUBLE.pt")
    # # FLOAT: torch.float32; DOUBLE:torch.double for using float/double in inference
    pes.to(device).to(torch.double)
    # # set the eval mode
    pes.eval()
    # # save the lattic parameters
    # cell=np.zeros((3,3),dtype=np.float)
    # period_table=torch.tensor([1,1,1],dtype=torch.float,device=device)   # same as the pbc in the periodic boundary condition
    # npoint=0
    # rmse=0
    # with open("./%s" %filename,'r') as f1:
    #     while True:
    #         string=f1.readline()
    #         if not string or npoint==1: break
    #         string=f1.readline()
    #         cell[0]=np.array(list(map(float,string.split())))
    #         string=f1.readline()
    #         cell[1]=np.array(list(map(float,string.split())))
    #         string=f1.readline()
    #         cell[2]=np.array(list(map(float,string.split())))
    #         string=f1.readline()
    #         species=[]
    #         cart=[]
    #         while True:
    #             string=f1.readline()
    #             if "abprop" in string: break
    #             tmp=string.split()
    #             tmp1=list(map(float,tmp[2:8]))
    #             cart.append(tmp1[0:3])
    #             species.append(atomtype.index(tmp[0]))
    #         #print(string.split())
    #         abene=float(string.split(':')[1])
    species=torch.from_numpy(np.array(species)).to(device)  # from numpy array to torch tensor
    cart=torch.from_numpy(np.array(cart)).to(device).to(torch.double)  # also float32/double
    tcell=torch.from_numpy(cell).to(device).to(torch.double)  # also float32/double
    energy,force=pes(period_table,cart,tcell,species)
    ene=energy.detach().cpu().numpy() # numpy() convert to the numpy array, cpu() from gpu to cpu
    #if abs(ene-abene) >=0.8:
    #    with open('err','a') as f:
    #        f.write('NN:%s   DFT: %s \n' %(ene,abene))
    force=force.detach().cpu().numpy() # numpy() convert to the numpy array, cpu() from gpu to cpu
    return ene, force

#print(np.sqrt(rmse/749))


