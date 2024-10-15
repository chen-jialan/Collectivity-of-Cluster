# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 16:43:42 2021

@author: DELL
"""

import os
import re
import copy
import math
import matplotlib.pyplot as plt
import numpy as np
from shutil import copyfile, rmtree
from tranfer import POSCAR_To_atom1
import pandas as pd
import random
# import random

def energy_all_str(path,filename):
    file_name = path+'/'+filename
    energy_all = []
    number_O = []
    number_CO = []
    gcmc_energy_all = []
    poscar = []
    O2 =  -10.71316/2
    CO =  -15.67154   
    with open(file_name,'r') as f:
        for line in f:
            number_O.append(float(line.strip().split()[0]))
            number_CO.append(float(line.strip().split()[1]))
            energy_all.append(float(line.strip().split()[2]))
            energy_gcmc = float(line.strip().split()[2]) + \
                                (20-float(line.strip().split()[0]))*O2+\
                                 (20-float(line.strip().split()[1]))*CO
            gcmc_energy_all.append(energy_gcmc)
            poscar.append(float(line.strip().split()[-1]))
    f.close()
    return energy_all, number_O, number_CO, gcmc_energy_all,poscar

def energy_POSCAR_all(path,filename):
    file_u = path + '/' + filename
    energy = []
    with open(file_u,'r') as f:
        for line in f:
            if 'real' in line.strip():
                energy.append(re.split(':', line.strip())[-1])
    f.close()
    return energy

def show(delta,P0):
    # x = [x0 for x0 in range(len(P0))]
    x = delta
    y = P0
    plt.rcParams['font.sans-serif'] = ['Times New Roman']
    plt.rcParams['axes.unicode_minus'] = False
    plt.rcParams['savefig.dpi'] = 1000
    plt.rcParams['figure.dpi'] = 1000
    params = {
    'figure.figsize': '6, 6'}
    plt.rcParams.update(params)
    
    font1 = {'family' : 'Times New Roman',
    'weight' : 'normal',
    'size'   : 14,
    }
    ax=plt.subplot(111)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    ax.set_xlabel(r'$\Delta$E/eV',fontsize=18)
    ax.set_ylabel('P',fontsize=18)
    plt.axhline(y=0.01,c="red",ls='--',lw=2)
    plt.xlim((0, 0.5))
    plt.xticks(size = 16)
    plt.yticks(size = 16)
    # model0 = make_interp_spline(x, y)
    # ys0 = model0(x)
    ax.scatter(x, y, color='b', marker='o', label="data")
    
##-------------read stuture on support----------------------
def support_atom(structure,Height):
    atom_dis = structure
    dictO = dict()
    dictCe = dict()
    H = []
    for pot in atom_dis.keys():
        if 'Cu-' in pot:
            b = atom_dis[pot][2]
            H.append(b)
    b = min(H)
        
    for pot in atom_dis.keys():
        if float(atom_dis[pot][2]) <= (Height):
            dictO[pot] = atom_dis[pot]    

    support = dict()
    support.update(dictCe)
    support.update(dictO)
    height = []
    for pot in dictO.keys():
        height.append(dictO[pot][2])
    max_H = max(height)
    return support, max_H

##-------------read cluster on stuture----------------------
def cluster_atom(structure, Height):
    atom_dis = structure
    dictCu = dict()
    dictO = dict()
    cluster_and_gas = dict()
    atom_dis_cluster = dict()
    cluster = dict()
#--------------confirm cluster&gas------------------
    for pot in atom_dis.keys():
        if 'Cu-' in pot:
            dictCu[pot] = atom_dis[pot]
            cluster[pot] = atom_dis[pot]
        else:
            if float(atom_dis[pot][2]) > Height:
                cluster[pot] = atom_dis[pot]
    return cluster
    
#--------------------Cu sur site ---------------------
#-------if Cu and other >= 2.5 ,then we define there are no site-------------
def site(support,cluster):
    support_dis = copy.deepcopy(support)
    cluster_dis = copy.deepcopy(cluster)
    site_all = dict()
    length_Cu_Cu = 2.34
    length_Cu_O = 1.98
    # print(length_Cu_Cu*1.2,length_Cu_O*1.2)
    for pot in cluster_dis.keys():
        if 'Cu' in pot:
            arr_pot = np.array(cluster_dis[pot])
            sur_pot = dict()
            sur_dis = []
            site_sur = []
            # site_all = []
            for pot1 in cluster_dis.keys():
                if not pot == pot1:
                    arr_pot1 = np.array(cluster_dis[pot1])
                    dis = float(np.linalg.norm(arr_pot1-arr_pot))
                    sur_pot[dis] = pot1
                    sur_dis.append(dis)
                    # site_all.append(pot1)
                    
            for pot2 in support_dis.keys():
                if not pot == pot2:
                    arr_pot2 = np.array(support_dis[pot2])
                    dis = float(np.linalg.norm(arr_pot2-arr_pot))
                    sur_pot[dis] = pot2
                    sur_dis.append(dis)
                    # site_all.append(pot1)
                 
            # sur_copy = copy.deepcopy(sur_dis)
            
            sur_dis.sort()
            
            if len(sur_dis) <=4:
                for sur in sur_dis:
                    if sur <= 2.5:
                        site_sur.append(sur_pot[sur])
                    
            else:
                for i in range(10):
                    sur = sur_dis[i]
                    if sur >  2.8:
                        continue
                    else:
                        if sur >  length_Cu_Cu*1.15 and 'Cu' in sur_pot[sur]:
                            continue
                        elif sur >  length_Cu_O*1.15 and 'O' in sur_pot[sur]:
                            continue
                        else:
                            site_sur.append(sur_pot[sur])
                    
            site_all[pot] = site_sur
    # print(site_all)
    
    return site_all

#-----------------determine the atom ------------------------
def sum_atom(site_all,atom_list,support,cluster):
    atom_type = dict()
    pot_cluster = [pot for pot in cluster.keys()]
    pot_support = [pot for pot in support.keys()]
    for pot in site_all.keys():
        atom_type_sur = dict()
        atom_type_new = dict()
        if len(site_all[pot]) != 0:
            for sur_pot in site_all[pot]:
                atom_name = re.findall(r'[A-Za-z]', sur_pot)
                if sur_pot in pot_cluster:
                    atom_name = str(''.join(atom_name))
                elif sur_pot in pot_support:
                    atom_name = str(''.join(atom_name)+'_s')
                if atom_name in atom_list:
                    if  atom_name in atom_type_sur.keys():
                        atom_type_sur[atom_name] += 1
                    else:
                        atom_type_sur[atom_name] = 1
        atom_type_new[pot] = atom_type_sur
        # if len(atom_type) == 0:
        atom_type[pot] = atom_type_sur
        # print(atom_type)
        # else:
        #     atom_type = new_site(atom_type_new,atom_type,0)
    return atom_type
            
#------------------new site------------------------------
def new_site(atom_type_new,atom_type_old,t):
    for pot_new in atom_type_new.keys():
        i = 0
        if atom_type_new[pot_new] != 0:
            for atom_new in atom_type_new[pot_new].keys():
                number_i_new = atom_type_new[pot_new][atom_new]
                for pot_old in atom_type_old.keys():
                    if  atom_new in atom_type_old[pot_old].keys():
                        number_i_old = atom_type_old[pot_old][atom_new]
                        if number_i_old == number_i_new:
                            i += 1
            if i != len(atom_type_new[pot_new]):
                atom_new_rename = '%s_%s' %(pot_new,t)
                atom_type_old[atom_new_rename] =  atom_type_new[pot_new]
        return atom_type_old

#--------------------get distance between site and next site----------- 
def distance_site_site(atom1,atom_cu,atomall):
    arr2 = np.array(atomall[atom1])
    dis = 1000
    atom_Cu = atom1
    for pot in atomall.keys():
        if pot in atom_cu and pot != atom1:
            arr1 = np.array(atomall[pot])
            dis_1 = float(np.linalg.norm(arr1-arr2))
            if dis_1 <= dis:
                dis = dis_1
                atom_Cu = pot
    return atom_Cu, dis


if __name__ == '__main__':
    path = os.getcwd()
    # filename = 'POSCAR-accept3'
    # fileall1 = ['energy_final1','energy_final2']
    # posfile = ['poscar1-1','poscar1-2']
    # fileall1 = ['energy_meta3','energy_meta_new','energy_meta_new2']
    fileall1 = ['energy_meta','energy_meta2'] #---------- energy-----------
    posfile = ['dataall1','dataall2'] #-----poscar in metastr---
    meta_str_file = 'file_metastr'
    if os.path.exists('./%s' %meta_str_file):
        rmtree('./%s' %meta_str_file)
    os.mkdir('./%s' %meta_str_file)
    atom_list = ['Cu','Ce(s)','O','C','O_s'] #----------sites---------
    H_fix = 1.6
    list_site_all = []
    list_site_num = []
    index_poscar = []
    # energy_ini = -806.91
    R = 8.618E-5
    temp = 400
    # energy_min = energy_ini
    x = 0
    # P = []
    energy_final = []
    file_final = []
    # finds = ['findstr3','findstr2','findstr1']
    finds = ['findstr1','findstr2'] #----------find structure---------
    for j in range(len(finds)):
        if os.path.exists('./%s' %finds[j]):
            rmtree('./%s' %finds[j])
        os.mkdir('./%s' %finds[j])
    # for i_all in range(3):
    #     P = []
    #     energy_all, number_O, number_CO, gcmc_energy_all,poscar = energy_all_str(path,fileall1[i_all])
    #     if min(energy_all) <= energy_min:
    #         energy_min = min(energy_all)
    for i_all in range(2):
        P = []
        energy_all, number_O, number_CO, gcmc_energy_all,poscar = energy_all_str(path,fileall1[i_all])
        energy_all_copy = copy.deepcopy(gcmc_energy_all)
        energy_all_copy.sort()
        # energy_min = float(min(energy_all_copy))
        # print(energy_all_copy)
        energy_min = min(energy_all_copy)
        for energy_e in energy_all_copy:
            if energy_e - energy_min <= 0.5 and energy_e >= energy_min:
                # print(i_all,energy_e-energy_min)
                P.append(math.exp(-(float(energy_e)-energy_min)/(temp*R)))
                x += 1
                # if len(energy_final) == 0:
                    # index_i = energy_all_copy.index(energy_e)
                    # energy_final.append(energy_e - energy_min)
                    # file_final.append(poscar[index_i])
                if energy_e - energy_min not in energy_final:
                    # abs((energy_e - energy_min)-energy_final[-1])>=0.01:
                    index_i = energy_all_copy.index(energy_e)
                    energy_final.append(energy_e - energy_min)
                    file_final.append(poscar[index_i])
            else:
                P.append(0)
        print(energy_final)
        # P = [x/np.sum(np.array(P)) for x in P]
        # P = [math.exp(-(float(energy_e)-energy_min)/(temp*R))  for energy_e in energy_all_copy]
        # detla = [-(float(energy_e)-energy_min) for energy_e in gcmc_energy_all]
        # show(detla,P)
        
        P_copy = copy.deepcopy(P)
        P_copy.sort(reverse=True)
        a = -1
        # print(P_copy)
        for i in range(len(P_copy)):
            if P_copy[i] != 0:
                a += 1
                index_i = P.index(P_copy[i])
                number_index = gcmc_energy_all.index(energy_all_copy[index_i])
                # print(abs(gcmc_energy_all[number_index]-energy_ini),a,P_copy[i]) 1
                # print(i_all,int(poscar[number_index]),energy_all_copy[i]-energy_min)
                if os.path.exists('./%s/POSCAR-%s' %(posfile[i_all],int(poscar[number_index]))):
                    copyfile('./%s/POSCAR-%s' %(posfile[i_all],int(poscar[number_index])), './%s/POSCAR-%s' %(finds[i_all],a))
                    copyfile('./%s/POSCAR-%s' %(posfile[i_all],int(poscar[number_index])), 'POSCAR')
                # os.system('cp pos/POSCAR-%s POSCAR' %number_index)
                # os.system('python3 dire2cart.py POSCAR;mv POSCAR_C POSCAR')
                allatom = POSCAR_To_atom1(path,'POSCAR')
                support,Height = support_atom(allatom,H_fix)
                cluster = cluster_atom(allatom,H_fix)
                site_all = site(support,cluster)
                atom_type = sum_atom(site_all,atom_list,support,cluster)
                
                atom_Cu_all = [pot for pot in allatom if 'Cu' in pot]
                
                # print(atom_type)
                for pot in atom_type.keys():
                    atom_Cu, dis = distance_site_site(pot,atom_Cu_all,allatom)
                    atom_type_new = [atom_type[pot],atom_type[atom_Cu]]
                    if not atom_type_new in list_site_all:
                        print(a,energy_all_copy[i],P_copy[i],i_all,int(poscar[number_index]))
                        if os.path.exists('./%s/POSCAR-%s' %(posfile[i_all],int(poscar[number_index]))):
                            copyfile('./%s/POSCAR-%s' %(posfile[i_all],int(poscar[number_index])), './%s/POSCAR-%s' %(meta_str_file,int(poscar[number_index])))
                        #     copyfile('%s/POSCAR-%s' %(posfile[i_all],int(poscar[number_index])), 'POSCAR')
                        # print(atom_type[pot])
                        list_site_all.append(atom_type_new)
                        list_site_num.append(P_copy[i])
                        index_poscar.append(int(poscar[number_index]))
                    else:
                        index_site = list_site_all.index(atom_type_new)
                        list_site_num[index_site] += P_copy[i]
                        
                        
        print(len(P),a)   
    

    # print(list_site_all,list_site_num)
    sum_of_site = []
    s = ['-','~']
    for i in range(len(list_site_all)):
        site_all = []
        # print(list_site_all[i],list_site_all[i][0])
        for j in range(len(list_site_all[i])):
            
            list_1 = []
            for pot in list_site_all[i][j].keys():
                name = '%s%s' %(list_site_all[i][j][pot],pot)
                list_1.append(name)
            site_e = s[0].join(list_1)
            site_all.append(site_e)
            
        # site_e_new = s[1].join(site_all)
            # print(list_site_all[i],site_e_new)
            # break
        # print(list_site_all[i],site_e_new)
        sum_of_site.append(site_all)
    # print(sum_of_site,list_site_num,len(sum_of_site),len(list_site_num))
    # T_rate = dict(zip(sum_of_site,list_site_num))
    # df = pd.DataFrame.from_dict(T_rate,orient='index')
    df = pd.DataFrame({'site1':[sum_of_site[i][0] for i in range(len(sum_of_site))],\
                        'site2':[sum_of_site[i][1] for i in range(len(sum_of_site))],\
                        'num':list_site_num,'poscar':index_poscar})
    # print(T_rate,df)
    print(df)
    if os.path.exists('./test_site6.xlsx'):
        os.remove('./test_site6.xlsx')
    df.to_excel('test_site6.xlsx')
    ten_la = [[] for i in range(3)]
    energy_ini = 0
    energy_copy = copy.deepcopy(energy_final)
    for i in range(len(energy_copy)):
        if i >=1:
            if energy_copy[i] - energy_ini <= 0.01:
                energy_final.remove(energy_copy[i])
            else:
                energy_ini = energy_copy[i]
        else:
            energy_ini = energy_copy[i]
    for i in range(len(energy_final)):
        ramdon_number = random.randint(0,2)
        ten_la[ramdon_number].append(energy_final[i])
    data_all = {'file':file_final}
    for i in range(3):
        data_all[i] = ten_la[i]
    df = pd.DataFrame.from_dict(data_all,orient='index')
    
    if os.path.exists('./test_energy6.xlsx'):
        os.remove('./test_energy6.xlsx')
    df.to_excel('test_energy6.xlsx')
    