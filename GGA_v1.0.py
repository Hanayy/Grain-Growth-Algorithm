#########################################################################
# Grain-grow algorithm v1.0
# code by Cheng-han Zhang
# E-mail: zhangchsr@163.com
# Jan,10,2024
#########################################################################


#########################################################################
# NOTE: After the code execution is complete, commands that 
#       can be directly used in the PFC2D 5.0 will be outputted
#       for assigning groups and extra information to each particle. 

# NOTE: Variable with the suffix 'vol' should be considered as area in 2D.

# NOTE: GGA need the particle information output by PFC.
## This code requires PFC to output specific content in a fixed format.
## In PFC2D v5.0, a .txt file can be output using the following code:
## ↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
## [textname = "filename.txt"]
## def savefile
##     ball_num = ball.num
##     file_infp = array.create(ball_num+2)
##     counting = 1
##     loop foreach bp ball.list
##         file_infp(counting) = string.build("%1 %2 %3 %4",ball.id(bp),ball.pos.x(bp),ball.pos.y(bp),ball.radius(bp))
##         counting += 1
##     endloop
##     outflag = file.open(textname,1,1)
##     file.write(file_infp,ball_num)
##     file.close()
## end
## @savefile
## ↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑
#########################################################################

#########################################################################
# 0. Import some packages
#########################################################################
import math
import random
import win32ui
import numpy as np
from shapely import geometry

#########################################################################
# 1. Set the parameters
#########################################################################
# For proportional change the number of grains.
multiplier = 1.0

# File that store the commands for PFC2D by GGA
commandFileName = 'Result-1'

# Random seed
random.seed(2)
np.random.seed(2)

# Number of grains
pl_num  = 15 
bt_num  = 22
kfs_num = 30
qtz_num  = 60

# The growth threshold, range [0,100].
pl_gt = 30
bt_gt = 30
kfs_gt = 30

# The angle range of growth direction.
# For random growth directions, this range is [0,180]. 
# Narrow this range to get a layer texture, for example, [44.9,45.1].
angle_range = [0,180]

# The aspect ratio of the assumed rectangle corresponding to the euhedral morphology.
pl_ar = 1.5
bt_ar = 4.5              
kfs_ar = 1.5          
qtz_ar = 1.0         

# The volume fraction of each type of mineral.
pl_vf = 0.3
bt_vf = 0.05
kfs_vf = 0.3
qtz_vf = 1-bt_vf-pl_vf-kfs_vf

# Initial length and growth rate. 
# These two values are related to the unit system of the basic model.
# For this code version, the values need to be relatively small to ensure a stable growth process. 
# In this case, the basic unit of length is meters.
ini_length = 0.0003
pl_grow_rate = 0.00008 

# As this parameter increases, the dispersion of particle sizes of similar mineral particles becomes larger. 
# The corresponding upper and lower quartiles are as follows:
# 1.1→[0.36,1.64] 1.3→[0.46,1.54] 1.5→[0.53,1.47] 2.0→[0.65,1.35] 4.0→[0.82,1.18]
# 10.0→[0.93,1.07] 20.0→[0.96,1.04] 50.0→[0.986,1.014] 100.0→[0.993,1.007] 
pl_dis = 2.0
bt_dis = 4.0 
kfs_dis = 2.0
qtz_dis = 2.0

# Width of the model
modelWidth = 0.05

# ↓↓↓Some variables depend on the above parameters.↓↓↓
number_of_pl  = int(pl_num*multiplier) 
number_of_bt  = int(bt_num*multiplier)
number_of_kfs = int(kfs_num*multiplier)
number_of_qtz  = int(qtz_num*multiplier)

bt_grow_rate = pl_grow_rate*number_of_pl/number_of_bt
kfs_grow_rate = bt_grow_rate*number_of_bt/number_of_kfs
qtz_grow_rate = kfs_grow_rate*number_of_kfs/number_of_qtz

lower_limit = angle_range[0]*np.pi/180       
upper_limit = angle_range[1]*np.pi/180    

randomNumber_pl = list(np.random.randn(1,number_of_pl)[0])
randomNumber_bt = list(np.random.randn(1,number_of_bt)[0])
randomNumber_kfs = list(np.random.randn(1,number_of_kfs)[0])
randomNumber_qtz = list(np.random.randn(1,number_of_qtz)[0])
randomNumber_list = [randomNumber_pl,randomNumber_bt,randomNumber_kfs,randomNumber_qtz]
randomMulti_random = [[],[],[],[]]

for i in range(len(randomNumber_list)):
    if i == 0:
        cr = pl_dis         
    elif i == 1:
        cr = bt_dis       
    elif i == 2:
        cr = kfs_dis
    elif i == 3:
        cr = qtz_dis
    min_num = cr*min(randomNumber_list[i])
    for rN in randomNumber_list[i]:
        new_num = (rN-min_num)/abs(min_num)
        randomMulti_random[i].append(new_num)

np.set_printoptions(suppress=True)   

#########################################################################
# 2. Select and read the PFC data file.
#########################################################################
dlg = win32ui.CreateFileDialog(1) 
dlg.SetOFNInitialDir('D:/') 
dlg.DoModal()
filename = dlg.GetPathName() 
ParticleData = open(filename,'r')
ParticleData_Clean = []
for line in ParticleData:
    data0 = line.replace('\n'," ")
    data1 = data0.split(' ')
    data2 = [int(data1[0]),float(data1[1]),float(data1[2]),float(data1[3])]
    ParticleData_Clean.append(data2)
ParticleData.close()
num_of_particles = len(ParticleData_Clean)

# This list is used to store all the balls that have been used.
usedList = []  
# This list is used to store all the balls that have not been used.
unusedList = []

# Total volume (area) and porosity of the basic model.
total_vol = 0.0
for i in ParticleData_Clean:
    vol = (i[-1])**2*np.pi
    total_vol += vol

ParticleInfo = []
for p in ParticleData_Clean:
    particleID = p[0]
    particleVol = (p[-1])**2*np.pi
    coor_of_ball = [p[1],p[2]]
    radius = p[-1]
    ParticleInfo.append([particleID,particleVol,coor_of_ball,radius])

#########################################################################
# 3. Build functions
#########################################################################
# Find the difference between the two lists
def diff_of_two(info,todelete):
    diff_list = []
    for item in info:
        if item not in todelete:
            diff_list.append(item)
    return diff_list

# Calculate the four vertices of the search region
def rectangle_vertex(ini_length,grow_rate,AR,ori_coor,direction,cycle):
    length = ini_length + grow_rate * cycle
    ori_x,ori_y = ori_coor[0],ori_coor[1]
    arc1,arc2 = direction[0],direction[1]
    ori_pointA = [ori_x+length*np.cos(arc1),ori_y+length*np.sin(arc1)]
    ori_pointB = [ori_x+length/AR*np.cos(arc2),ori_y+length/AR*np.sin(arc2)]
    point_a = [length*np.cos(arc1),length*np.sin(arc1)]
    point_b = [length/AR*np.cos(arc2),length/AR*np.sin(arc2)]
    ori_pointI = [ori_pointA[0]+point_b[0],ori_pointA[1]+point_b[1]]
    ori_pointJ = [-point_a[0]+ori_pointB[0],-point_a[1]+ori_pointB[1]]
    ori_pointK = [ori_x-(point_a[0]+point_b[0]),ori_y-(point_a[1]+point_b[1])]
    ori_pointL = [ori_pointA[0]-point_b[0],ori_pointA[1]-point_b[1]]
    return [ori_pointI,ori_pointJ,ori_pointK,ori_pointL]  

# Determine which particle are included in the growing grain
def graingrow(mineral_list,ggMulti_list,infor_list,ini_length, \
              the_grow_rate,AR,this_cycle,unusedList, \
              save_list,usedList,acc_vol):
    for i in range(len(mineral_list)):
        ggMulti = ggMulti_list[i]
        four_vertex = rectangle_vertex(ini_length,ggMulti*the_grow_rate,AR,
            [mineral_list[i][0],mineral_list[i][1]],[mineral_list[i][2],mineral_list[i][3]],this_cycle)
        tri1 = [tuple(four_vertex[0]),tuple(four_vertex[1]),
                tuple(four_vertex[2]),tuple(four_vertex[0])]    
        tri2 = [tuple(four_vertex[2]),tuple(four_vertex[3]),
                tuple(four_vertex[0]),tuple(four_vertex[2])]                
        check_tri1 = geometry.Polygon(tri1)                          
        check_tri2 = geometry.Polygon(tri2)                          
        ori_co = [mineral_list[i][0],mineral_list[i][1]]
        for pp in unusedList:
            permit = False                                       
            p_ID = infor_list[pp-1][0]
            p_vol = infor_list[pp-1][1]
            p_x = infor_list[pp-1][2][0]
            p_y = infor_list[pp-1][2][1]
            p_cen = geometry.Point(tuple(infor_list[pp-1][2]))
            p_radius = infor_list[pp-1][3]
            in_tri1 = check_tri1.covers(p_cen)
            in_tri2 = check_tri2.covers(p_cen)
            if in_tri1 == True or in_tri2 == True:
                if p_ID not in usedList:   
                    if len(save_list[i]) == 0:
                        permit = True                      
                    else:    
                        for ob in save_list[i]:
                            sg_ParticleRadius = infor_list[ob-1][3]
                            sg_Particle_x = infor_list[ob-1][2][0]
                            sg_Particle_y = infor_list[ob-1][2][1]
                            distance = ((sg_Particle_x-p_x)**2+(sg_Particle_y-p_y)**2)**0.5
                            if distance < sg_ParticleRadius+p_radius:
                                permit = True 
                    if permit == True:     
                        save_list[i].append(p_ID)
                        usedList.append(p_ID) 
                        unusedList.remove(p_ID)
                        acc_vol += p_vol
    return acc_vol

# Coordinates and growth direction of grains
def seedAndDirection(unusedList,number_of_mineral,ParticleInfo,\
                     lower_limit,upper_limit,mineral_list):
    #print('sam',unusedList,number_of_mineral)
    the_random = random.sample(unusedList, number_of_mineral)  
    the_seed = []
    for pp in the_random:
        coor = ParticleInfo[pp-1][2]
        arc1 = random.uniform(lower_limit,upper_limit)
        if arc1 >= 0.5*np.pi:
            arc2 = arc1 - 0.5*np.pi
        else:
            arc2 = arc1 + 0.5*np.pi
        the_seed.append([coor[0],coor[1],arc1,arc2])  
    for i in range(len(the_seed)):
        mineral_list.append([])
    return the_random,the_seed

#########################################################################
# 4. Main process of grain growth
#########################################################################
all_particle_ID = list(np.linspace(1,num_of_particles,num_of_particles,dtype=int))
unusedList = diff_of_two(all_particle_ID,usedList)

# Used to store information about contained particles.
pl_list = []              
bt_list = []               
kfs_list = []              
qtz_list = []               

# Actual growth rate.
pl_ggMulti = randomMulti_random[0]
bt_ggMulti = randomMulti_random[1]
kfs_ggMulti = randomMulti_random[2]
qtz_ggMulti = randomMulti_random[3]

# Initialize variables.
pl_vol = 0.
bt_vol = 0.
kfs_vol = 0.
qtz_vol = 0.

# Assign the value to the mineral first to grow.
pl_random,pl_seed = seedAndDirection(unusedList,number_of_pl,ParticleInfo,\
                                     lower_limit,upper_limit,pl_list)

# Initialize growth state.
pl_can_grow = True
bt_can_grow = False
kfs_can_grow = False
qtz_can_grow = False
bt_random = []
kfs_random = []
qtz_random = []

# Store the number of unused particles after each loop to determine whether to break the loop.
judgeIfStop = []

# Initialize cycle number.
bt_start_cycle = 0
kfs_start_cycle = 0
qtz_start_cycle = 0
total_cycle = 0

# Grain growth progress
while len(unusedList) > 0:
    judgeIfStop.append(len(unusedList))
    total_cycle += 1
    pl_acc = pl_vol/total_vol*100
    bt_acc = bt_vol/total_vol*100
    kfs_acc = kfs_vol/total_vol*100
    qtz_acc = qtz_vol/total_vol*100
    total_acc = (pl_vol+bt_vol+kfs_vol+qtz_vol)/total_vol*100
    print('Pl: '+str("%6.2f" %pl_acc)+'% 丨 '+'Bt: '+str("%6.2f" %bt_acc)+'% 丨 ' \
          'Kfs: '+str("%6.2f" %kfs_acc)+'% 丨 '+'Qtz: '+str("%6.2f" %qtz_acc)+'% 丨 ' \
          'total: '+str("%6.2f" %total_acc)+'%  ')
    if pl_can_grow == True:
        if pl_vol < pl_vf*total_vol:
            pl_vol = graingrow(pl_seed,pl_ggMulti,ParticleInfo,ini_length, \
                                pl_grow_rate,pl_ar,total_cycle,unusedList, \
                                pl_list,usedList,pl_vol)
            if pl_vol < pl_gt/100*pl_vf*total_vol:
                bt_start_cycle = total_cycle
            else:
                if len(bt_random) == 0:
                    bt_random,bt_seed = seedAndDirection(unusedList,number_of_bt,ParticleInfo,\
                                         lower_limit,upper_limit,bt_list)
                bt_can_grow = True
        else:
            pl_can_grow = False
    if bt_can_grow == True:
        bt_cycle = total_cycle - bt_start_cycle
        if bt_vol < bt_vf*total_vol:
            bt_vol = graingrow(bt_seed,bt_ggMulti,ParticleInfo,ini_length, \
                                bt_grow_rate,bt_ar,bt_cycle,unusedList, \
                                bt_list,usedList,bt_vol)
            if bt_vol < bt_gt/100*bt_vf*total_vol:
                kfs_start_cycle = total_cycle
            else:
                if len(kfs_random) == 0:
                    kfs_random,kfs_seed = seedAndDirection(unusedList,number_of_kfs,ParticleInfo,\
                                         lower_limit,upper_limit,kfs_list)
                kfs_can_grow = True
        else:
            bt_can_grow = False
    if kfs_can_grow == True:
        kfs_cycle = total_cycle - kfs_start_cycle
        if kfs_vol < kfs_vf*total_vol:
            kfs_vol = graingrow(kfs_seed,kfs_ggMulti,ParticleInfo,ini_length, \
                                kfs_grow_rate,kfs_ar,kfs_cycle,unusedList, \
                                kfs_list,usedList,kfs_vol)
            if kfs_vol < kfs_gt/100*kfs_vf*total_vol:
                qtz_start_cycle = total_cycle
            else:
                if len(qtz_random) == 0:
                    qtz_random,qtz_seed = seedAndDirection(unusedList,number_of_qtz,ParticleInfo,\
                                         lower_limit,upper_limit,qtz_list)
                qtz_can_grow = True
        else:
            kfs_can_grow = False
    if qtz_can_grow == True:
        qtz_cycle = total_cycle - qtz_start_cycle
        if qtz_vol < qtz_vf*total_vol:
            qtz_vol = graingrow(qtz_seed,qtz_ggMulti,ParticleInfo,ini_length, \
                                qtz_grow_rate,qtz_ar,qtz_cycle,unusedList, \
                                qtz_list,usedList,qtz_vol)

# Break the cycle if the number of unused particles does not decrease for 20 consecutive cycles
    if len(judgeIfStop) >= 20:        
        if judgeIfStop[-1] == judgeIfStop[-20]:
            break   

#########################################################################
# 5. Grouping of the unused particles
#########################################################################
while len(unusedList) != 0:
    print('There are '+str(len(unusedList))+' particles still unuse.')
    if len(unusedList) == 0:
        break
    if pl_vol < pl_vf*total_vol:
        if len(unusedList) == 0:
            break
        while pl_vol < pl_vf*total_vol:
            pl_re_cycle = 0
            if len(unusedList) == 0:
                break
            pl_list_re = []
            pl_random_re,pl_seed_re = seedAndDirection(unusedList,1,ParticleInfo,\
                                                       lower_limit,upper_limit,pl_list_re)            
            num_in_pl = []
            for i in pl_list:
                num_in_pl.append(len(i))
            ave_pl_grain = int(np.mean(num_in_pl))
            the_length = 0.
            while the_length < ave_pl_grain:
                pl_re_cycle += 1
                pl_vol = graingrow(pl_seed_re,[1],ParticleInfo,ini_length, \
                                    pl_grow_rate,pl_ar,pl_re_cycle,unusedList, \
                                    pl_list_re,usedList,pl_vol)
                if pl_re_cycle >= modelWidth*0.15/pl_grow_rate:  
                    pl_list.append(pl_list_re[0])
                    break
    if bt_vol < bt_vf*total_vol:
        if len(unusedList) == 0:
            break
        while bt_vol < bt_vf*total_vol:
            bt_re_cycle = 0
            if len(unusedList) == 0:
                break
            bt_list_re = []
            bt_random_re,bt_seed_re = seedAndDirection(unusedList,1,ParticleInfo,\
                                                       lower_limit,upper_limit,bt_list_re)            
            num_in_bt = []
            for i in bt_list:
                num_in_bt.append(len(i))
            ave_bt_grain = int(np.mean(num_in_bt))
            the_length = 0.
            while the_length < ave_bt_grain:
                bt_re_cycle += 1
                bt_vol = graingrow(bt_seed_re,[1],ParticleInfo,ini_length, \
                                    bt_grow_rate,bt_ar,bt_re_cycle,unusedList, \
                                    bt_list_re,usedList,bt_vol)
                if bt_re_cycle >= modelWidth*0.15/bt_grow_rate:   
                    bt_list.append(bt_list_re[0])
                    break                    
    if kfs_vol < kfs_vf*total_vol:
        if len(unusedList) == 0:
            break
        while kfs_vol < kfs_vf*total_vol:
            kfs_re_cycle = 0
            if len(unusedList) == 0:
                break
            kfs_list_re = []
            kfs_random_re,kfs_seed_re = seedAndDirection(unusedList,1,ParticleInfo,\
                                                       lower_limit,upper_limit,kfs_list_re)            
            num_in_kfs = []
            for i in kfs_list:
                num_in_kfs.append(len(i))
            ave_kfs_grain = int(np.mean(num_in_kfs))
            the_length = 0.
            while the_length < ave_kfs_grain:
                kfs_re_cycle += 1
                kfs_vol = graingrow(kfs_seed_re,[1],ParticleInfo,ini_length, \
                                    kfs_grow_rate,kfs_ar,kfs_re_cycle,unusedList, \
                                    kfs_list_re,usedList,kfs_vol)
                if kfs_re_cycle >= modelWidth*0.15/kfs_grow_rate:
                    kfs_list.append(kfs_list_re[0])
                    break
    if qtz_vol < qtz_vf*total_vol:
        if len(unusedList) == 0:
            break
        while qtz_vol < qtz_vf*total_vol:
            qtz_re_cycle = 0
            if len(unusedList) == 0:
                break
            qtz_list_re = []
            qtz_random_re,qtz_seed_re = seedAndDirection(unusedList,1,ParticleInfo,\
                                                       lower_limit,upper_limit,qtz_list_re)            
            num_in_qtz = []
            for i in qtz_list:
                num_in_qtz.append(len(i))
            ave_qtz_grain = int(np.mean(num_in_qtz))
            the_length = 0.
            while the_length < ave_qtz_grain:
                qtz_re_cycle += 1
                qtz_vol = graingrow(qtz_seed_re,[1],ParticleInfo,ini_length, \
                                    qtz_grow_rate,qtz_ar,qtz_re_cycle,unusedList, \
                                    qtz_list_re,usedList,qtz_vol)
                if qtz_re_cycle >= modelWidth*0.15/qtz_grow_rate:  
                    qtz_list.append(qtz_list_re[0])
                    break

#########################################################################
# 6. Output the result of the particle grouping.
#########################################################################

with open((str(commandFileName)+'.txt'), 'w') as file:
    for g in range(len(pl_list)):
        for balls in pl_list[g]:
            file.write("ball group pl_"+str(g+1)+" range id "+str(balls)+ '\n')
    for g in range(len(bt_list)):
        for balls in bt_list[g]:
            file.write("ball group bt_"+str(g+1)+" range id "+str(balls)+ '\n')
    for g in range(len(kfs_list)):
        for balls in kfs_list[g]:
            file.write("ball group kfs_"+str(g+1)+" range id "+str(balls)+ '\n')
    for g in range(len(qtz_list)):
        for balls in qtz_list[g]:
            file.write("ball group qtz_"+str(g+1)+" range id "+str(balls)+ '\n')
    for l in range(len(pl_list)):
        file.write("ball extra 1 pl range group 'pl_"+str(l+1)+"'"+ '\n')
        file.write("ball extra 2 '"+str(l+1)+"' range group 'pl_"+str(l+1)+"'"+ '\n')    
    for l in range(len(bt_list)):
        file.write("ball extra 1 bt range group 'bt_"+str(l+1)+"'"+ '\n')
        file.write("ball extra 2 '"+str(l+1)+"' range group 'bt_"+str(l+1)+"'"+ '\n')    
    for l in range(len(kfs_list)):
        file.write("ball extra 1 kfs range group 'kfs_"+str(l+1)+"'"+ '\n')
        file.write("ball extra 2 '"+str(l+1)+"' range group 'kfs_"+str(l+1)+"'"+ '\n')    
    for l in range(len(qtz_list)):
        file.write("ball extra 1 qtz range group 'qtz_"+str(l+1)+"'"+ '\n')
        file.write("ball extra 2 '"+str(l+1)+"' range group 'qtz_"+str(l+1)+"'"+ '\n')    


print('All commands generated by GGA for PFC2D have been stored in '+str(commandFileName)+'.txt')