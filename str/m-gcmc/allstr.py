# -*- coding: utf-8 -*-
"""
Created on Sat Mar  6 22:50:25 2021

@author: DELL
"""

import os
import shutil
import re
# import linecache

##allstr.acr 分解成单个的结构#####
def mkdir_file(path_name):
    path1 = path_name + '/' + 'pos'
    if not os.path.exists('%s' %path1):
        os.mkdir('%s' % path1)
        
    a = 0
    path11 = []
    files = path_name+'/str'
    for name in os.listdir(files):
        filename = '%s/%s' %(files,name)
        path11.append(filename)
        
    ##路径排序###    
    path11.sort()
    # print(path11)
    
    for i in range(len(path11)):
        file = path11[i]
        filename1 = os.path.basename(file)
        with open(file,'r') as fp:
            lines = fp.readlines()
    
        for line in lines:
            if 'all' in line:
                a += 1
            
            newfile = path1 + '/POSCAR-%s' %(a)
            with open(newfile,'a')as f:
                f.write(line)
                
        a = 0
            
    ##遍历文档####
    #for file2 in os.listdir(path1):
    #    file2 = path1 + '/%s' %file2
        # print(file2)
    #    with open(file2,'r+') as f:
    #        content = f.read()
    #        f.seek(0, 0)
    #        f.write('PBC=ON\n'+content)
    #    with open(file2,'r+') as f:
    #        content = f.read()
    #        f.seek(0, 0)
    #        f.write('!BIOSYM archive 2\n'+content)
            
    #####OSZICAR判断合理结构####
def mkdir_OSZ(path_name):
    path1 = path_name + '/OSZall'
    if not os.path.exists('%s' %path1):
        os.mkdir('%s' % path1)
        
    a = 0
    path11 = []
    files = path_name + '/OSZ'
    for name in os.listdir(files):
        filename = '%s/%s' %(files,name)
        path11.append(filename)
        
    ##路径排序###    
    path11.sort()
    # print(path11)
    
    for i in range(len(path11)):
        file1 = path11[i]
        filename1 = os.path.basename(file1)
        with open(file1,'r') as fp:
            lines = fp.readlines()
    
        for line in lines:
            if 'rms' in line:
                a += 1
            
            newfile = path1 + '/%s-%s' %(filename1,a)
            with open(newfile,'a')as f:
                f.write(line)
        a = 0
    

def True_OSZ(path_name):
    path1 = path_name + '/OSZall'
    path11 = []
    # files = path_name + '/OSZ'
    for name in os.listdir(path1):
        filename = '%s/%s' %(path1,name)
        path11.append(filename)
        
    ##路径排序###    
    path11.sort()
    # print(path11)
    C = []
    D0 = []
    D1 = []
    b = 0
    # c = 0
    for file in path11:
        # file = path11[i]
        # c +=1
        with open(file,'r') as fp:
            lines = fp.readlines()
    
        for line in lines:
            b += 1
        # print(b)
            
        if  b >= 62:
            b = 0
        else:
            # C.append(c)
            filename1 = os.path.basename(file)
            count1 = re.findall(r'([\d+\.]+)',filename1)
            D0.append(int(count1[0]))
            D1.append(int(count1[1]))
            C.append(count1)
            b = 0
            
    E0 = [D0[0]]
    E1 = [D1[0]]
    a = 0
    for i in range(1,len(D0)):
        for j in range(len(E0)):
            if D0[i] == E0[j]:
                if D1[i] > E1[j]:
                    a = 1
                    # break
                # else:
                #     E0.append(D0[i])
                #     E1.append(D1[i])
            else:
                continue
            
        if a == 0:
            E0.append(D0[i])
            E1.append(D1[i])
            a = 0
        else:
            a = 0

    
    G = (E0,E1)
    #print(G)
    return G


def get_goodstr(path_name,G):
    D0 = G[0]
    D1 = G[1]
    D2 = G[2]
    path1 = path_name + '/goodstr'
    path2 = path_name + '/allstr1'
    if not os.path.exists('%s' %path1):
        os.mkdir('%s' % path1)
    
    for i in range(len(D0)):
        #oldfile = path2 +'/allstr.arc-%s-%s.arc' %(D0[i],D1[i])
        oldfile = path2 +'/allstr.arc-%s.arc-%s.arc' %(D0[i],D2[i])
        newfile = path1 +'/allstr-%s-%s.arc' %(D1[i],D2[i])
        shutil.copyfile(oldfile, newfile)
    

    
def mkdir_allfor(path_name):
    path1 = path_name + '/allfor1'
    if not os.path.exists('%s' %path1):
        os.mkdir('%s' % path1)
        
    a = 0
    path11 = []
    files = path_name + '/For'
    for name in os.listdir(files):
        filename = '%s/%s' %(files,name)
        path11.append(filename)
        
    ##路径排序###    
    path11.sort()
    # print(path11)
    
    for i in range(len(path11)):
        file = path11[i]
        filename1 = os.path.basename(file)
        with open(file,'r') as fp:
            lines = fp.readlines()
    
        for line in lines:
            if 'For' in line:
                a += 1
            
            newfile = path1 + '/%s-%s' %(filename1,a)
            with open(newfile,'a')as f:
                f.write(line)
        a = 0
        
        
def get_goodsfor(path_name,G):
    path1 = path_name + '/goodfor'
    path2 = path_name + '/allfor1'
    if not os.path.exists('%s' %path1):
        os.mkdir('%s' % path1)
    
    D0 = G[0]
    D1 = G[1]
    for i in range(len(D0)):
        oldfile = path2 +'/allfor.arc-%s-%s' %(D0[i],D1[i])
        newfile = path1 +'/allfor-%s-%s.arc' %(D0[i],D1[i])
        shutil.copyfile(oldfile, newfile)
    
def get_tranfer(filename,C):
    list1  =  list()
    with open(filename,'r') as f:
        for  line in  f:
            list1.append(line.strip())
    number_i = []
    number_p = []
    for number_all in list1:
        num_all = re.split(r"[ ]+", number_all)
        number_i.append(int(num_all[0]))
        number_p.append(int(num_all[1])) 
     
    num_tran = dict(zip(number_i,number_p))           
    num_real = num_tran[C]
    return num_real
 
def main():
    Path = os.getcwd()
    mkdir_file(Path)
    #mkdir_allfor(Path)
    #mkdir_OSZ(Path)
    #C = True_OSZ(Path)
    #D0 = []
    #for i in range(len(C[0])):
    #     D0.append(get_tranfer('file_c',C[0][i]))
    #D = (C[0],D0,C[1])
    #D11 = (D0,C[1])
    #print(D)
    #get_goodstr(Path,D)
    #get_goodsfor(Path, D11)
    

if __name__ == main():
    main()
