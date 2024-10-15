# -*- coding: utf-8 -*-
"""
Created on Sat May  8 19:07:28 2021

POSCAR to allatom and to POSCAR
@author: DELL
"""
import os
import re
import linecache
import copy

#--------------allstr.arc to atom(dict)---------------
def LASP_To_atom1(path, filename):
    list1  =  list()
    file = path + '/' + filename
    with open(file,'r') as fp:
        for  line in  fp:
            list1.append(line.strip())
    a = len(list1)
    at_dec = []
    atom = []
    for i in range(5,a-2):
        text = list1[i]
        line_all = re.findall(r"-?\d+\.?\d*e?-?\d*?",text)
        split1 =  '%s-%s' %(str(text.split(' ',1)[0]),i-4)
        atom.append(split1)

        x = float(line_all[0])
        y = float(line_all[1])
        z = float(line_all[2])
        
        at_dec.append([x,y,z])

    structure = dict(zip(atom,at_dec))
    return structure

#--------------POSCAR to atom(dict)---------------
def POSCAR_To_atom1(path, filename):
    list1  =  list()
    with open('%s/%s' %(path,filename),'r') as fp:
        for  line in  fp:
            list1.append(line.strip().replace('\n','').replace('\t',''))

    atom_all = re.split(r"[ ]+", list1[5])
    atom_num = re.split(r"[ ]+", list1[6])
    #print(atom_all,atom_num)
    atom_and_num = dict(zip(atom_all,atom_num))
    s = 9
    at_dec1 = []
    atom = []
    for kind_num in atom_and_num.keys():
        s_old = s
        s += int(atom_and_num[kind_num])
        for dis in range(s_old,s):
            distance = re.split(r"[ ]+", list1[dis])
            atom.append('%s-%s' %(kind_num, dis-6))
            dis_all1 = [float(distance[0]),float(distance[1]),float(distance[2])]
            at_dec1.append(dis_all1)
    structure1 = dict(zip(atom,at_dec1))
    return structure1

#---------------------lammps to atom------------
def lammps2atom(path, filename):
    list1  =  list()
    file = path + '/' + filename
    with open(file,'r') as fp:
        for  line in  fp:
            list1.append(line.strip())
    a = len(list1)
    at_dec = []
    atom = []
    atom_type = ['Cu','Ce','O','C']
    for i in range(9,a):
        text = list1[i]
        line_all = re.split(r"[ ]+", text)
        #print(line_all)
        atom_number = int(line_all[1])-1
        at = atom_type[atom_number]
        split1 =  '%s-%s' %(at,int(line_all[0]))
        atom.append(split1)
        x = float(line_all[2])
        y = float(line_all[3])
        z = float(line_all[4])
        at_dec.append([x,y,z])

    structure = dict(zip(atom,at_dec))
    return structure

#---------------------eann to atom------------
def eann2atom(path,filename):
    list1  =  list()
    a = 0
    line_all = []
    file1 = path + '/' + filename
    with open(file1,'r') as f:
        for  line in  f:
            a += 1
            if 'pbc' in line:
                line_all.append(a)
            elif 'abprop' in line:
                line_all.append(a)
                break
    atom_all = dict()
    x = 0
    #print(line_all)
    for line in range((int(line_all[0])+1),line_all[1]):
        x +=1
        txt = linecache.getline(file1,line)
        all_info = re.split(r"[ ]+", txt)
        atom = all_info[0]
        distance = [float(all_info[2]),float(all_info[3]),float(all_info[4])]
        atom_all['%s-%s'%(atom,x)] = distance
    linecache.clearcache()
    return atom_all


#---------------------eann to atom------------
def eann2force(path,filename):
    list1  =  list()
    a = 0
    line_all = []
    file1 = path + '/' + filename
    with open(file1,'r') as f:
        for  line in  f:
            a += 1
            if 'pbc' in line:
                line_all.append(a)
            elif 'abprop' in line:
                line_all.append(a)
                break
    atom_all = dict()
    x = 0
    #print(line_all)
    for line in range((int(line_all[0])+1),line_all[1]):
        x +=1
        txt = linecache.getline(file1,line)
        all_info = re.split(r"[ ]+", txt)
        atom = all_info[0]
        force = [float(all_info[-3]),float(all_info[-2]),float(all_info[-1])]
        atom_all['%s-%s'%(atom,x)] = force
    linecache.clearcache()
    return atom_all

        

#--------------atom to input.arc-------------------
def new_input(atom_all, filename, newfilename):
    for line in range(1,6):
        start = linecache.getline('./%s' %filename,line)
        with open('./%s' %newfilename, 'a') as fp1:
            fp1.write(start)
    a = 1
    for pot1 in atom_all.keys():
        atom = re.findall(r'[A-Za-z]', pot1)
        atom = ''.join(atom)
        with open('./%s' %newfilename, 'a') as fp1:
            fp1.write('%s  %15.9f%15.9f%15.9f CORE %s  %s %s 0.0000 %s\n' %(atom.ljust(3),float(atom_all[pot1][0]), float(atom_all[pot1][1]), float(atom_all[pot1][2]), str('%4.0f'%(a)), atom.ljust(3), atom.ljust(3), str('%4.0f'%(a))))
        a +=1
    with open('./%s' %newfilename, 'a') as f:
        f.write('end\n')
        f.write('end')
    linecache.clearcache()

#--------------atom to POSCAR---------------
def write_POSCAR(atom_all_copy, filename, newfilename,H):
    atom_all = copy.deepcopy(atom_all_copy)
    for line in range(6):
        start = linecache.getline('./%s' %filename,line)
        with open('./%s' %newfilename, 'a') as fp1:
            fp1.write(start)
    
    #atom_all = dict(sorted(atom_all.items(), key=lambda e: e[1]))
    dictCu = dict()
    dictCe = dict()
    dictC = dict()
    dictO = dict()
    atom_list = ['Cu','Ce','O','C-']
    atom_dis = [dictCu,dictCe,dictO,dictC]
    for i,atom in enumerate(atom_list):
        for pot in atom_all.keys():
            if atom in pot:
                atom_dis[i][pot] = atom_all[pot]
    atom_all_new = dict()
    for pot_atom in atom_dis:
        atom_all_new.update(pot_atom)
    atom_all = atom_all_new
    #print(atom_all)
    atom = {atom_list[0]:0}
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
    #print(atom) 
    for pot_fix in atom_all.keys():
        if atom_all[pot_fix][2] > H:
            for i in range(3):
                atom_all[pot_fix].append('T')
        else:
            for i in range(3):
                atom_all[pot_fix].append('F')
    with open('./%s' %newfilename, 'a') as fp1:
        for pot in atom.keys():
        #    print(pot)
            fp1.write('%s ' %pot)
        fp1.write('\n')
        for pot in atom.keys():
    #        print(pot,int(atom[pot]))
            fp1.write('%s ' %int(atom[pot]))
        fp1.write('\n')
    
    for line in range(8,9):
        start = linecache.getline('./%s' %filename,line)
        with open('./%s' %newfilename, 'a') as fp1:
            fp1.write(start)
    with open('./%s' %newfilename, 'a') as fp1:
        fp1.write('Cartesian\n')
    a = 1
    
    for pot1 in atom_all.keys():
        dis = atom_all[pot1]
        with open('./%s' %newfilename, 'a') as fp1:
            fp1.write('%.5f  %.5f  %.5f    %s  %s  %s\n' %(dis[0],dis[1],dis[2],dis[3],dis[4],dis[5]))
        a +=1
    linecache.clearcache()
    fp1.close()

#-------------------write For from allfor.arc--------------
def write_For(path,filename):
    list1  =  list()
    with open('%s/%s' %(path,filename),'r') as fp:
        for  line in  fp:
            list1.append(line.strip())
    at_F1 = []
    for lines in range(2,len(list1)):
        Force = re.split(r"[ ]+", list1[lines])
        if '' not in Force:
            at_F1.append([float(Force[0]),float(Force[1]),float(Force[2])])
    
    return at_F1

#-----------------write all to eann-train-------------------
def write_data_eann(path,atom_str,Force,E):
    file_train = path +'/train'
    with open(file_train,'a+') as fp:
        fp.write('point \n')
        for line in range(4):
            start = linecache.getline('./%s' %'lat',line)
            fp.write(start)
        fp.write('pbc   1    1    1\n')
        a = 0
        atom_all = dict()
        for pot1 in atom_str.keys():
            atom = ''.join(re.findall(r'[A-Za-z]', pot1))
            structure1 = atom_str[pot1]
            F = Force[a]
            if 'C' == atom:
                mass = 12.01
            elif 'Ce' == atom:
                mass = 140.1
            elif 'Cu' == atom:
                mass = 63.55
            elif 'O' == atom:
                mass = 16.00
            a +=1
            atom_num = '%s-%s' %(atom,a)
            atom_all[atom_num] = [float(mass), structure1[0],structure1[1],structure1[2],F[0],F[1],F[2]]
        # atom_all = dict(sorted(atom_all.items(), key=lambda e: e[1]))
        list_atom = ['Cu-','Ce-','C-','O-']
        for atom_e in list_atom:
            for pot1 in atom_all.keys():
                if atom_e in pot1:
                    atom = ''.join(re.findall(r'[A-Za-z]', pot1))
                    structure1 = atom_all[pot1]
                    fp.write('%s  %.4f  %.5f  %.5f  %.5f    %s  %s  %s\n' %(atom, structure1[0], structure1[1],structure1[2],structure1[3],structure1[4],structure1[5],structure1[6]))
        fp.write('abprop:%s\n' %E)
    linecache.clearcache()
    fp.close()

#------------------------------------write data without F------------------------
def write_withoutF_eann(path,atom_str,E):
    file_train = path +'/train'
    with open(file_train,'a+') as fp:
        fp.write('point \n')
        for line in range(4):
            start = linecache.getline('./%s' %'lat',line)
            fp.write(start)
            # print(start)
        fp.write('pbc   1    1    1\n')
        a = 0
        atom_all = dict()
        for pot1 in atom_str.keys():
            atom = ''.join(re.findall(r'[A-Za-z]', pot1))
            structure1 = atom_str[pot1]
            if 'C' == atom:
                mass = 12.01
            elif 'Ce' == atom:
                mass = 140.1
            elif 'Cu' == atom:
                mass = 63.55
            elif 'O' == atom:
                mass = 16.00
            a +=1
            atom_num = '%s-%s' %(atom,a)
            atom_all[atom_num] = [float(mass), structure1[0],structure1[1],structure1[2]]
        # atom_all = dict(sorted(atom_all.items(), key=lambda e: e[1]))
        list_atom = ['Cu-','Ce-','O-','C-']
        for atom_e in list_atom:
            for pot1 in atom_all.keys():
                if atom_e in pot1:
                    atom = ''.join(re.findall(r'[A-Za-z]', pot1))
                    structure1 = atom_all[pot1]
                    fp.write('%s  %.4f  %.5f  %.5f  %.5f\n' %(atom, structure1[0], structure1[1],structure1[2],structure1[3]))
        fp.write('abprop:%s\n' %E)
    fp.close()
    linecache.clearcache()

#------------------------get Energy from POSCAR--------------
def Energy_VASP():
    #exa=os.popen('grep "energy  without entropy" OUTCAR | tail -1 ','r')
    exa=os.popen('grep F= OSZICAR |tail -1','r')
    line=exa.readline()
    if len(line) != 0:
        energy=float(line.split()[4])
    else:
        energy=0
    return energy

#-----------------------get Energy from LASP-------------- 
def Energy_LASP(path,filename):
    file1 = path + '/%s' %filename
    list1 = linecache.getline('%s' %file1 , 3).strip()
    energy = float(re.split(r"[ ]+", list1)[-1])
    linecache.clearcache()
    return energy

#-----------------------get Energy from lammps (dict(number_E))--------------
def Energy_lammps(path,filename):
    file1 = path + '/%s' %filename
    a = 0
    num_E = dict()
    line_E = []
    #with open(file1,'r',encoding='gb18030', errors='ignore') as fp:
    with open(file1,'r') as fp:
        for  line in  fp:
            a += 1
            if 'Step TotEng' in line:
                line_E.append(a+1)
                continue
            elif 'Loop time of' in line:
                line_E.append(a-1)
                break
    #num_E = dict()
    #-----------update file --------------------
    linecache.checkcache(file1)
    for line in range((int(line_E[0])),int(line_E[1])):
        txt = linecache.getline(file1,line).strip() 
        number_e = re.split(r"[ ]+", txt)
        #print(number_e)
        num_E[number_e[0]] = number_e[1]
    fp.close() 
    #-------------------clear ----------------------
    linecache.clearcache() 
    return num_E
                
#-------------------abs(energy) value is too large, we estimate it is bad-str
def badstr_energy(energy):
    #energy = []
    energy_abs = []
    x = 0
    if len(energy) >= 2:
        #for number in Energy.keys():
        #    energy.append(float(Energy[number]))
        for i in range(len(energy)-1):
            b = abs(energy[i+1] - energy[i]) 	
            energy_abs.append(b)
        energy_max_min = abs(max(energy)-min(energy))
        if not energy_max_min>= 10:
            for i, energy_a in enumerate(energy_abs):
                if energy_a >= 2 and i >= 1:
                    if energy_abs[i-1] <= energy_a:
                        x += 1
                elif energy_a >= 5:
                    x += 1
        else:
            x += 1
    return x

def lat(path,filename):
    fileposcar = path + '/' + filename
    with open('./lat','a+',encoding='utf-8') as ff:
        for line in range(3,6):
            ff.write('%s' %(linecache.getline(fileposcar,line)))
    linecache.clearcache()

def main():
    path = os.getcwd()
    filename = 'out'
#    atom = lammps2atom(path, filename)
 #   print(atom)
#    H = 2
#    filename_lasp = 'allstr-1-2.arc'
#    filename_eann = 'train'
#    atom_all = LASP_To_atom1(path, filename_lasp) 
    #atom_all = eann2atom(path,filename_eann) 
    #E = Energy_LASP(path,filename_lasp)
#    filename = 'POSCAR_base_3_3'
#    newname = 'E-lammps'
#    newfilename = 'POSCAR-1'
    #print(atom_all)
#    write_POSCAR(atom_all, filename, newfilename,H)
    #write_data_eann(path,atom_all,Force,E)
    #write_withoutF_eann(path,atom_all,E)
    number_E = Energy_lammps(path,filename)
    print(number_E)
    #print(number_E)
    # filename = 'mini.lammpstrj'
    # atom = lammps2atom(path, filename)
    # print(atom)
    
#----------------------vasp2lammps--------------------
#    os.system('awk -f VASP-poscar2lammps.awk %s > datal' %newfilename)
#if __name__ == main():
#    main()
