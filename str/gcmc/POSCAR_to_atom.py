# -*- coding: utf-8 -*-
"""
Created on Sat May  8 19:07:28 2021

POSCAR to allatom and to POSCAR
@author: DELL
"""
import os
import re
import linecache


def POSCAR_To_atom1(path, filename):
    list1  =  list()
    with open('%s/%s' %(path,filename),'r') as fp:
        for  line in  fp:
            list1.append(line.strip())

    atom_all = re.split(r"[ ]+", list1[5])
    atom_num = re.split(r"[ ]+", list1[6])
    atom_and_num = dict(zip(atom_all,atom_num))
    s = 9
    at_dec1 = []
    # at_dec2 = []
    atom = []
    for kind_num in atom_and_num.keys():
        s_old = s
        s += int(atom_and_num[kind_num])
        for dis in range(s_old,s):
            distance = re.split(r"[ ]+", list1[dis])
            atom.append('%s-%s' %(kind_num, dis-6))
            # dis_all2 = [float(distance[0]),float(distance[1]),float(distance[2]),distance[3],distance[4],distance[5]]
            dis_all1 = [float(distance[0]),float(distance[1]),float(distance[2])]
            at_dec1.append(dis_all1)
            # at_dec2.append(dis_all2)
    structure1 = dict(zip(atom,at_dec1))
    # structure2 = dict(zip(atom,at_dec2))
    return structure1

def write_POSCAR(atom_all, filename, newfilename):
    for line in range(6):
        start = linecache.getline('./%s' %filename,line)
        with open('./%s' %newfilename, 'a') as fp1:
            fp1.write(start)
    
    atom_all = dict(sorted(atom_all.items(), key=lambda e: e[1]))
    
    atom = {list(atom_all.items())[0].keys():0}
    # print(atom_all)
    for pot in atom_all.keys():
        number = 1
        for pot1 in atom.keys():
            if pot1 in pot:
                atom[pot1] = 1 + float(atom[pot1])
                break
            else:
                number += 1
        if number > len(atom):
            atom_new = re.findall(r'[A-Za-z]', pot)
            atom_new = ''.join(atom_new)
            atom['%s' %atom_new] = 1
    
    with open('./%s' %newfilename, 'a') as fp1:
        for pot in atom.keys():
            fp1.write('%s \t' %pot)
        fp1.write('\n')
        for pot in atom.keys():
            fp1.write('%s \t' %int(atom[pot]))
        fp1.write('\n')
        
    for line in range(8,10):
        start = linecache.getline('./%s' %filename,line)
        with open('./%s' %newfilename, 'a') as fp1:
            fp1.write(start)
    a = 1
    
    for pot1 in atom_all.keys():
        dis = atom_all[pot1]
        with open('./%s' %newfilename, 'a') as fp1:
            fp1.write('%.5f  %.5f  %.5f    %s  %s  %s\n' %(dis[0],dis[1],dis[2],dis[3],dis[4],dis[5]))
        a +=1
    return atom

def Energy_VASP():
    exa=os.popen('grep "energy  without entropy" OUTCAR | tail -1 ','r')
    line=exa.readline()
    energy=float(line[64:])
    E=energy
    return E

def Energy_VASP_new():
    exa = os.popen("grep 'F=' OSZICAR | tail -1 | awk '{print $5}'").readline()
    energy=float(exa)
    E=energy
    return E

# def main():
#     path = os.getcwd()
#     filename = 'POSCAR'
#     all_atom = POSCAR_To_atom1(path,filename)
#     # os.system('python dire2cart.py POSCAR')
#     write_POSCAR(all_atom,'POSCAR','POSCAR_1')
    
# if __name__ == main():
#     main()