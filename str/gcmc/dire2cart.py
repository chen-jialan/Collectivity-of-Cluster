#!/usr/bin/env python
# -*- coding: utf-8 -*-
#Convert direc coordiation to cartesian Writen By Qiang and Xu Nan
"""
The right input way to run this script is : 
      
python dire2cart.py CONTCAR 3 

CONTCAR will be converted from Direct to Cartesian and 3 layers will be fixed.  

If you do not want to fix layers and keep the same as before, run command:
    
python dire2cat.py CONTCAR    

"""

import sys
import os
import numpy as np 

#print ('Please read the head part of this script and get more information!')
#print ("""
###################################
#                                 #
#for VASP 5.2 or higher versions  #
#                                 #
###################################
#""")

if len(sys.argv) <= 1:
    print('\n' + ' Warning! ' * 3 + '\n')
    print('You did not select the inputfile to be converted. \n By defalut, we are going to convert your CONTCAR.\n')
    
    if not os.path.isfile("POSCAR") and not os.path.isfile("CONTCAR"):
        print("Error:" * 3 + "\n Can not find neither POSCAR nor CONTCAR!!\n")
        exit()
    else:
        if os.path.isfile("CONTCAR"):
            script = sys.argv[0]
            file_to_be_converted = "CONTCAR"
            print ("\n Convertion starts........\n")
    
        else:
            script = sys.argv[0]
            file_to_be_converted = "POSCAR"
            print ("\n There is no CONTCAR in your directory. \n \n POSCAR Convertion starts........\n")
            
if len(sys.argv) == 2:
    print("\n%s Convertion starts......" %sys.argv[1]) 
    script, file_to_be_converted = sys.argv[:2]
else:
    print("\n%s Convertion starts......\n" %sys.argv[1]) 
    script, file_to_be_converted, fixedlayer = sys.argv[:3]
    fixedlayer = int(fixedlayer)

def get_infor():
    f = open(file_to_be_converted, 'r')
    lines = f.readlines()
    f.close()
    num_atoms = sum([int(x) for x in lines[6].split()])
    if lines[7][0]  == 'S' or lines[7][0]  == 's':  # # With Selected T T T, coordination starts from line 9
        start_num = 9
        if lines[8][0]  == 'D' or lines[8][0]  == 'd':
            is_direct = True 
        else:
            is_direct = False
    else: 
        start_num = 8
        print ('----------------------------------------------------')
        print ('Pay Attetion! There is no TTT in  %s   ' %(file_to_be_converted))
        print ('----------------------------------------------------')
        if  lines[7][0]  == 'D' or lines[7][0]  == 'd':
            is_direct = True 
        else:
            is_direct = False
    a = []
    b = []
    c = []
    if is_direct:         
        for i in np.arange(2,5):
            line = [float(i) for i in lines[i].strip().split()]
            a.append(line[0])
            b.append(line[1])
            c.append(line[2])
        vector = np.array([a,b,c])
    else: 
        vector = np.array([[1, 0 , 0], [0, 1, 0], [0, 0, 1]])
    return vector, lines, start_num, num_atoms, is_direct

def determinelayers(z_cartesian):
    threshold = 0.5  
    seq = sorted(z_cartesian)
    min = seq[0]
    layerscount = {}
    sets = [min]
    for j in range(len(seq)):
        if abs(seq[j]-min) >= threshold:
            min = seq[j]
            sets.append(min)

    for i in range(1,len(sets)+1):
        layerscount[i] = []            
        for k in range(len(z_cartesian)):   
            if abs(z_cartesian[k]-sets[i-1]) <= threshold:
                layerscount[i].append(k)

    return layerscount

def convert():
    x_cartesian = []
    y_cartesian = []
    z_cartesian = []
    tf = []
    for i in range(start_num, num_atoms + start_num):
        line_data =  [float(ele) for ele in lines[i].split()[0:3]]
        line_data = np.array([line_data])
        x, y, z = [sum(k) for k in line_data * vector ]     
        x_cartesian.append(x)
        y_cartesian.append(y)
        z_cartesian.append(z)
         
        if start_num == 9 : # if  T T T exist, the start_num will be 9        
            tf.append((lines[i].split()[3:]))
        else:
            tf.append(' ')   # if there is no T T T, use space instead. 

    layerscount =determinelayers(z_cartesian)
    
    file_out = open(file_to_be_converted+'_C', 'w')
    for i in range(0,7):
        file_out.write(lines[i].rstrip() + '\n')  # first 7 lines are kept the same 
    
#    print ('\n Find %s layers!' * 3 %(len(layerscount), len(layerscount), len(layerscount)))
    
    if len(sys.argv) >= 3:  # This means that the nuber for fixing layers is given by the user.
        file_out.write('Selective\n')  
        file_out.write('Cartesian' + '\n') 
        for i in range(1,len(layerscount)+1):
            if i <= fixedlayer: 
                for j in layerscount[i]:
                    tf[j] = ['F','F','F']
            else:
                for k in layerscount[i]:
                    tf[k] = ['T','T','T']                    
    else:
        if start_num == 9:  # 9 means there are T T T or F F F in the file
            file_out.write('Selective\n')
            file_out.write('Cartesian' + '\n')        
        else:
            file_out.write('Cartesian' + '\n') 
    
    for i in range(0,len(x_cartesian)):
        file_out.write("\t%+-3.10f   %+-3.10f   %+-3.10f  %s\n" %(x_cartesian[i], y_cartesian[i], z_cartesian[i], ' '.join(tf[i])))
        
    file_out.close()
    
vector, lines, start_num, num_atoms, is_direct = get_infor()

if is_direct : 
#    print ("\n%s has Direct Coordinates, Contersion starts.... "  %(file_to_be_converted))
    convert()
else:
#    print ("\n%s has Cartesian Coordinates Already! We are going to fix layers only."  %(file_to_be_converted))
    convert()
    
#print ('-----------------------------------------------------\n')
#print ('\n %s with Cartesian Coordiates is named as  %s_C\n' %(file_to_be_converted, file_to_be_converted))
#print ('-----------------------------------------------------\n')

