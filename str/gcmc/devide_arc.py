# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 18:21:40 2021

@author: DELL
"""

import os
import shutil
import time
# import linecache

##allstr.acr 分解成单个的结构#####
def mkdir_file(path_name,filename,filename2):
    path1 = path_name + '/%s' %filename2
    if not os.path.exists('%s' %path1):
        os.mkdir('%s' % path1)
    a = 0
    # path = os.getcwd()
    file = path_name + '/%s' %filename
    
    with open(file,'r') as fp:
        lines = fp.readlines()
    
    for line in lines:
        if 'TIMESTEP' in line:
            a += 1
            
        newfile = path1 + '/%s' %a
        with open(newfile,'a')as f:
            f.write(line)
            
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
    #os.system("for i in %s/*.arc;do sed -i '5s/22/17/g' $i;done" %path1)

###放到文件夹中
# def put_in_file(path_name):
#     path1 = path_name + '\\str'
    
#     ####统计文件夹下多少个文件######
#     a = len([name for name in os.listdir(path1) if os.path.isfile(os.path.join(path1, name))])
#     # print(a)
#     for i in range(a):
#         if not os.path.exists('%s\%s' %(path1,i)):
#             os.mkdir('%s\%s' % (path1,i))
#         oldfile1 = '%s\%s.arc' %(path1,i)
#         # print(oldfile1)
#         newfile1 = '%s\%s\input.arc' %(path1,i)
#         shutil.move(oldfile1, newfile1)
        
###创建para放入KPOINTS、POTCAR、lasp.in、slurmlasp.sh、INCAR
def other_file(path_name, file_old):
    # path1 = file_old + '\str'
    path2 = file_old + '/para'
    path_new = path_name
    # a = len([name for name in os.listdir(path1) if os.path.isfile(os.path.join(path1, name))])
    for name in os.listdir(path2):
        # print(name)
        # for i in range(1,a):
        filename = '%s/%s' %(path2,name)
        filename_new = '%s/%s' %(path_new,name)
        shutil.copyfile(filename, filename_new)

def sbatch(path_name,path_old):
    Path = path_name
    father_path = os.getcwd()
    Path_arc = Path + '/str'
    Path_for = Path + '/For'
    Path_OSZ = Path + '/OSZ'
    # print(Path_arc)
    if not os.path.exists('%s' %Path_arc):
        os.mkdir('%s' %Path_arc)
    if not os.path.exists('%s' %Path_for):
        os.mkdir('%s' %Path_for)
    if not os.path.exists('%s' %Path_OSZ):
        os.mkdir('%s' %Path_OSZ)
    path1 = path_old + '/str'
    # print(path1)
    a = len([name for name in os.listdir(path1) if os.path.isfile(os.path.join(path1, name))])
    # print(a)
    for i in range(1,a):
        oldfile1 = '%s/%s.arc' %(path1,i)
        # print(oldfile1)
        newfile1 = '%s/input.arc' %Path
        # print(newfile1)
        shutil.copyfile(oldfile1, newfile1)
        os.chdir('%s' %Path)
        os.system('mpirun /public/home/zzucae09/lasp/INTER_pro2.3.0_intel12_2cw/Src/lasp > log')
	#os.system('mpirun /public/home/zzucae09/lasp/INTER_pro2.3.0_intel12_2cw/Src/lasp')
        os.chdir('%s' %father_path)
        oldfile2 = '%s/allstr.arc' %Path
        newfile2 = '%s/allstr.arc-%s' %(Path_arc,i)
        shutil.move(oldfile2, newfile2)
        oldfile3 = '%s/allfor.arc' %Path
        newfile3 = '%s/allfor.arc-%s' %(Path_for,i)
        shutil.move(oldfile3, newfile3)
        oldfile4 = '%s/OSZICAR' %Path
        newfile4 = '%s/OSZICAR-%s' %(Path_OSZ,i)
        shutil.move(oldfile4, newfile4)
        time.sleep(2)
    
 
def main():
    ##单个目录
    #file_name = str(input('allstr绝对路径:'))
    filepathZ = os.getcwd()
    filec = 9
    filepath_1 = '%s/%s' %(filepathZ,filec)
    for root,dirs,files in os.walk(filepath_1):
        if 'allstr.arc' in files:
            filepath = root
            mkdir_file(filepath)
            # put_in_file(filepath)
            file = os.getcwd()
            if not os.path.exists('%s-comp_1' %(filepath)):
                os.mkdir('%s-comp_1' %(filepath))
            file_comp = '%s-comp_1' %(filepath)
            other_file(file_comp,file)
            sbatch(file_comp,filepath)
        else:
            continue 
    ###所有目录
    # file_name = str(input('allstr绝对路径:'))
    # # file_name2 = str(input('参数的绝对路径:'))
    # for file in os.listdir(file_name):
    #     filepath = '%s\%s' %(file_name, file)
    #     for root,dirs,files in os.walk(filepath):
    #         if 'allstr.arc' in files:
    #             filepath = root
    #             mkdir_file(filepath)
    #             # put_in_file(filepath)
    #             if not os.path.exists('%s-1' %file):
    #                 os.mkdir('%s-1' %file)
    #             file_comp = '%s-1' %file
    #             other_file(file_comp,file_name)
    #             sbatch(file_comp,filepath)
    #             print('------%s is end------' %(filepath))
    #         else:
    #             continue
        
        

    
    
if __name__ == main():
    main()
            
