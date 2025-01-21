#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 20 18:29:31 2024

@author: Matthew
"""


"""
Created on Fri May  3 21:23:10 2024

@author: Matthew
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  3 19:39:26 2024

@author: Matthew
"""
'''imports necessary functions'''
import canalizing_function_toolbox_v2_0 as can2
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import networkx as nx
import canalizing_function_toolbox_v13 as can13
import random
import csv
import pandas as pd



'''Modify Networks'''
def SortList(List): #Sort a list from lowest to highest number

    
    #Intial Variables
    j = 0
    k = len(List)
    SortedList = list()
    
    
    while j < k: #Do an iteration for every value in List
    
        min_value = 1000 #Set min_value to some large number
        i = 0 #Set index to 0
        
        
        while i < len(List): #For every value in List
            if List[i] < min_value: #Find minimum value
                min_value = List[i]
            i += 1
            
            
        List.remove(min_value) #Remove that value from the List
        SortedList.append(min_value) #Add that value to the sorted list
        
        j += 1 #Iterate
        
    return SortedList

def PercentOfOnes(F, I):
    k = can2.num_of_attractors(F, I, len(I), EXACT = True)
    i = 0
    List_of_attractors = k[0]
    TotalDigits = 0
    TotalOnes = 0
    Temp = 0
    
    while i < len(k[0]):
        Attractors = List_of_attractors[i]
        j = 0
        while j < len(List_of_attractors[i]):
            binList = can13.dec2bin(Attractors[j])
            l = 0
            while l < len(binList):
                if binList[l] == 0:
                    Temp = 0
                else:
                    TotalOnes += 1
                l += 1
            TotalDigits += len(I)
            j += 1
        i += 1
    return TotalOnes/TotalDigits

def DrawWiringDiagram(NetworkI): #Draw a wiring diagram
    matrix = can2.adjacency_matrix(NetworkI)
    G = nx.from_numpy_array(matrix, parallel_edges = False, create_using = nx.MultiDiGraph())
    layout = nx.drawing.layout.fruchterman_reingold_layout
    pos_G=layout(G)
    node_size=100
    nodes= np.array([pos_G[node] for node in list(G.nodes)]).T
    title = 'wiring diagram'
    f,ax = plt.subplots()
    nx.draw_networkx(G,pos=pos_G,with_labels=True,ax=ax,node_size=node_size,edge_color='#C1C1C1')
    [x1,x2] = ax.get_xlim()
    [y1,y2] = ax.get_ylim()
    ax.set_title(title)
    ax.plot(nodes[0],nodes[1],'o',color=cm.tab10(0),markersize=5)
    
def ILIST(I):
    #Creates list with all values of I
    i = 0
    I_TempList = list()
    
    while i < len(I):
        I_TempList = np.concatenate((I[i], I_TempList), axis = None)
        i += 1
    I_TempList = list(I_TempList)
    
    
    #Removes duplicates from I_TempList
    I_list = list()
    i = 0
    
    while i < len(I_TempList):
        if I_TempList[i] in I_list:
            i += 1
        else:
            I_list.append(int(I_TempList[i]))
            i += 1
    
    #Outputs I_list
    return I_list
    
def ShiftNetworkI(N, Network1, startIndex, endIndex):
    index = 0
    network1 = Network1.copy() #Ensure Network1 isn't changed
    while index < endIndex - startIndex + 1: #Shift NetworkI between Start and End Indices
        network1[startIndex + index] = network1[startIndex + index] + N*(np.ones(len(network1[index]), dtype = int)) #Shift each value in network[index] by N
        index += 1
    return network1

def JoinTwoNetworks(Network1, Network2, Type): #Type 0 is F, Type 1 is I, Type 0 is if it has already been used

    network1 = Network1.copy()
    network2 = Network2.copy()
    
    if Type == 0:
        return network1 + network2 #Just combine two functions
    
    
    elif Type == 1:
        return network1 + ShiftNetworkI(len(network1), network2, 0, len(network2) - 1) #Combines two I's but shifts the second one automatically to ensure no overlap
    
def RandomlyConnectNetworks(I_1, I_2, num_of_edges): #Connections flow from I_1 to I_2
    
    NewI_1 = I_1.copy()
    NewI_2 = I_2.copy() #Makes sure original I_1 and I_2 are not changed
    
    
    I_1list = ILIST(NewI_1) #Get a list of all the nodes in I
    i = 0
    
    
    if num_of_edges == -1:
        while i < len(NewI_2):
            k = random.randint(1, len(I_1list)) #Take a random sample of nodes from I_1
            j = SortList(random.sample(I_1list, k))
            NewI_2[i] = np.concatenate((NewI_2[i], j), axis = None) #Connect those nodes to a specific node in I_2
            i += 1
    else:
        
        
        
        while i < num_of_edges:
            
            
            k = random.randint(0, len(NewI_2)-1)
            j = random.sample(I_1list, 1) #Pick a random node
            l = list(NewI_2[k].copy())
            
            if l.count(j[0]) < 1: #check to make sure the node isn't already on NewI_2[k]
            
            
                NewI_2[k] = np.concatenate((NewI_2[k], j), axis = None) #Connect those nodes to a specific node in I_2
                i += 1
                
                
            else:
                i = i
        
    return JoinTwoNetworks(NewI_2, NewI_1, Type = 0) #Join I_1 and I_2

def AddOnesToF(F, nodeInF, numOfOnes): 
    i = 0
    if numOfOnes > 0: #make sure numOfOnes is not 0 to ensure no errors
        while i < numOfOnes:
            F[nodeInF] = np.concatenate((F[nodeInF], [1]), axis = None) #add a bunch of ones, specifically numOfOnes to the end of F[index]
            i += 1
    else:
        F[nodeInF] = F[nodeInF] #Just don't do anything if numOfOnes is 0
    return F

def UpdateF(F, I):  #Designed for canalyzing functions, F_combined in the later part of the code
    i = 0
    if len(F) == len(I): #make sure that the length of F and I are the same
        while i < len(I): #while loop for every value in F
                if 2**len(I[i]) - len(F[i]) > 0: #Checks to see if len(F[index]) is less than 2^len(I[index]) 
                    F = AddOnesToF(F, i, 2**len(I[i]) - len(F[i])) #Add enough ones to make len(F[index]) equal to 2^(len(I[index]))
                    i += 1
                    
                else:
                    i +=1 #If len(F[index]) = 2^(lenI[index]) then just continue
    else:
        i += 1 #Don't do anything if len(F) is not len(I)

def SortCanalyzingFunctions(F_1_Attractors, F_2_Attractors, AttractorList, List_of_F_1, List_of_F_2, List_of_I_1, List_of_I_2):
    #index 
    i = 0
    
    #Lists to return
    F_1 = list()
    F_2 = list()
    
    #Attractor List should be structured [Attractors of F_1, Attractors of F_2]
    
    #Iterate through Attractor List
    while i < 1:
        TempList = AttractorList[i]
        #index
        i = 0 
        
        #Lists to return
        F_1 = list()
        F_2 = list()
        I_1 = list()
        I_2 = list()

        #Attractor List should be structured with each entry being [Attractors of F_1, Attractors of F_2]
        #To allow F_1_Attractors or F_2_Attractors to be anything, let them be -1

        
        #Iterate through all values of AttractorList
        while i < len(AttractorList) - 1:
            
            #Create temp list to access entry in Attractor List
            TempList = AttractorList[i]
            
            #Check to see if either F_1_Attractors or F_2_Attractors is equal to -1
            if F_1_Attractors > 0 and F_2_Attractors > 0:
                if TempList[0] == F_1_Attractors and TempList[1] == F_2_Attractors:
                    F_1.append(List_of_F_1[i])
                    F_2.append(List_of_F_2[i])
                    I_1.append(List_of_I_1[i])
                    I_2.append(List_of_I_2[i])
                    
            elif F_1_Attractors == -1 and F_2_Attractors > 0:
                if TempList[1] == F_2_Attractors:
                    F_1.append(List_of_F_1[i])
                    F_2.append(List_of_F_2[i])
                    I_1.append(List_of_I_1[i])
                    I_2.append(List_of_I_2[i])
                
            if F_1_Attractors > 0 and F_2_Attractors == -1:
                if TempList[0] == F_1_Attractors:
                    F_1.append(List_of_F_1[i])
                    F_2.append(List_of_F_2[i])
                    I_1.append(List_of_I_1[i])
                    I_2.append(List_of_I_2[i])
            i += 1
        
    return F_1, I_1, F_2, I_2

    



'''Generate Networks'''

def GenerateCanalyzingNetworks(N_1, n_1, N_2, n_2):
    #Intital variables
    i = 0 
    l_Num_of_Attractors = list()
    List_of_Num_of_Attractors = list()
    List_of_F_1 = list()
    List_of_I_1 = list()
    List_of_F_2 = list()
    List_of_I_2 = list()
    List_of_F_combined = list()
    List_of_I_combined = list()
        
    #Loop to find number of attractors 100 times
    while i < 100:

        #Create Two Random Networks aand Find Their Attractors

        l_Num_of_Attractors = list()
        [F_1, I_1, degree_1] = can13.random_BN(N = N_1, n = n_1)
        k = 0
        k = can2.num_of_attractors(F_1, I_1, len(F_1), EXACT = True)
        l_Num_of_Attractors.append(k[1])
        List_of_F_1.append(F_1)
        List_of_I_1.append(I_1)

        
        
        [F_2, I_2, degree_2] = can13.random_BN(N = N_2 , n = n_2)
        k = 0
        k = can2.num_of_attractors(F_2, I_2, len(F_2), EXACT = True)
        l_Num_of_Attractors.append(k[1])
        List_of_F_2.append(F_2)
        List_of_I_2.append(I_2)
        

        
        #Join the Networks Randomly
        
        i_1 = I_1.copy()
        i_1 = ShiftNetworkI(len(I_2), i_1, 0, len(i_1) - 1)
        F_combined = JoinTwoNetworks(F_2, F_1, 0)
        I_combined = RandomlyConnectNetworks(i_1, I_2, 1)  #Connections flow from I_1 to I_2
        UpdateF(F_combined, I_combined)
        List_of_F_combined.append(F_combined)
        List_of_I_combined.append(I_combined)
        
        
        #Find the Attractors
        
        k = 0
        k = can2.num_of_attractors(F_combined, I_combined, len(F_combined), EXACT = True)
        l_Num_of_Attractors.append(k[1])
        
        
        #Repeat
        
        List_of_Num_of_Attractors.append(l_Num_of_Attractors)
        i += 1 
        
    
   # Returns List of Attractors
    
    
    return List_of_Num_of_Attractors, List_of_F_1, List_of_I_1, List_of_F_2, List_of_I_2, List_of_F_combined, List_of_I_combined
        
def DefinedCanFunction(Num_of_Connections): #Uses Preset F_1 and F_2, connections flow from F_1 to F_2, 
    F_1 = [[0, 0, 1, 1, 0, 1, 1, 0],
           [0, 1, 1, 1, 1, 0, 1, 0],
           [1, 1, 1, 0, 0, 0, 0, 1],
           [1, 0, 0, 1, 1, 0, 0, 0],
           [1, 1, 1, 0, 1, 1, 1, 1],
           [1, 1, 0, 1, 1, 0, 1, 1],
           [1, 1, 1, 0, 0, 0, 0, 1],
           [0, 1, 1, 1, 1, 0, 1, 0],
           [1, 1, 0, 1, 1, 1, 1, 1],
           [0, 0, 1, 0, 1, 1, 0, 0]]
    I_1 = [[2, 4, 8],
           [0, 7, 9], 
           [6, 7, 9], 
           [0, 5, 8], 
           [0, 5, 7], 
           [0, 2, 8], 
           [2, 3, 7],
           [2, 3, 8], 
           [0, 1, 7],
           [3, 4, 5]]
    F_2 = [[0, 0, 0, 1], 
           [0, 1, 0, 0], 
           [1, 1, 0, 1], 
           [1, 0, 1, 1], 
           [0, 1, 1, 0]]
    
    I_2 = [[2, 3], 
           [0, 3], 
           [1, 4], 
           [1, 2], 
           [2, 3]]
    
    I_1 = ShiftNetworkI(len(I_2), I_1, 0, len(I_1) - 1) #Shift I_1
    F_combined = JoinTwoNetworks(F_2, F_1, 0) #Join the F's
    I_combined = RandomlyConnectNetworks(I_1, I_2, Num_of_Connections)  #Connections flow from I_1 to I_2
    UpdateF(F_combined, I_combined) #Update F to be canalyzing
   
    
    return can2.num_of_attractors(F_combined, I_combined, 15, EXACT = True), F_combined, I_combined

def SetCanFunction(f_1, f_2, i_1, i_2, Num_of_Connections): #Uses Preset F_1 and F_2, connections flow from F_1 to F_2, 
    F_1 = f_1.copy()
    I_1 = i_1.copy()
    F_2 = f_2.copy()
    I_2 = i_1.copy()
    
    I_1 = ShiftNetworkI(len(I_2), I_1, 0, len(I_1) - 1) #Shift I_1
    F_combined = JoinTwoNetworks(F_2, F_1, 0) #Join the F's
    I_combined = RandomlyConnectNetworks(I_1, I_2, Num_of_Connections)  #Connections flow from I_1 to I_2
    UpdateF(F_combined, I_combined) #Update F to be canalyzing
   
    
    return can2.num_of_attractors(F_combined, I_combined, 15, EXACT = True), F_combined, I_combined

def RandomCanFunction(N_1, n_1, N_2, n_2, Num_of_Connections, specific):
    #Intital variables
    
    #Temp Variable, used to add stuff to lists
    k = 0
    Num_of_Attractors = list()
    
    #Percennt of Ones List
    BN_1_Percent_Of_Ones = 0

        

    #Create Random BN-1 And Find Its Attractors

    
    #Define BN-1
    [F_1, I_1, degree_1] = can13.random_BN(N = N_1, n = n_1) 
    
    #Calculate attractors of BN-1
    k = can2.num_of_attractors(F_1, I_1, len(F_1), EXACT = True)
    
    #Add attractors to temp list
    Num_of_Attractors.append(k[1])
    Specific_Attractors_F1 = k[0]
    
    
    #Add Percent of Ones to List
    BN_1_Percent_Of_Ones = PercentOfOnes(F_1, I_1)
    
    
    
    


    #Create Random BN-2 And Find Its Attractors

    
    #Define BN-2
    [F_2, I_2, degree_2] = can13.random_BN(N = N_2 , n = n_2)
    
    #Calculate attractors for BN-2
    k = can2.num_of_attractors(F_2, I_2, len(F_2), EXACT = True)
    
    #Add attractors to temp list
    Num_of_Attractors.append(k[1])
    Specific_Attractors_F2 = k[0]
        
       
        




    #Join the Networks Randomly

    
    #Ensure I_1 is not altered
    i_1 = I_1.copy()
    
    #Shift I_1 so it can be properly joined at the end of I_2
    i_1 = ShiftNetworkI(len(I_2), i_1, 0, len(i_1) - 1)
    
    #Join F-1 and F-2
    F_combined = JoinTwoNetworks(F_2, F_1, 0)
    
    #Join I_1 and I_2
    I_combined = RandomlyConnectNetworks(i_1, I_2, Num_of_Connections)  #Connections flow from I_1 to I_2
    
    #Correct the number of 1's to make F_combined canalyzing
    UpdateF(F_combined, I_combined)
    
    #Calulate the number of attractors of BN-combined
    k = can2.num_of_attractors(F_combined, I_combined, len(F_combined), EXACT = True)
    
    #Add attractors to list
    Num_of_Attractors.append(k[1])
    Specific_Attractors_F_combined = k[0]
        
        
        
        

    #Return

    if specific == True:
        return F_1, F_2, I_1, I_2, F_combined, I_combined, Num_of_Attractors, BN_1_Percent_Of_Ones, Specific_Attractors_F1, Specific_Attractors_F2, Specific_Attractors_F_combined
    else:
        return Num_of_Attractors, BN_1_Percent_Of_Ones




'''General SImulations'''
def SetCanSimulation(f_1, f_2, i_1, i_2):
    j = 1
    Average_List = []
    
    while j < 21:
        Temp_List = []
        i = 0 
        l = []
        while i < 1000:
            [n, F, I] = SetCanFunction(f_1, i_1, f_2, i_2, j)
            l.append(n[1])
            i += 1
            
        sum_values = sum(l)
        average = sum_values /1000
        Temp_List.append(j)
        Temp_List.append(average)
        Average_List.append(Temp_List)
        j += 1

    # Write Average_List to CSV
    with open('simulation_results.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Index', 'Average'])  # Header row
        writer.writerows(Average_List)  # Data rows

    return Average_List


def DefinedCanSimulation():
    j = 1
    Average_List = []
    
    while j < 21:
        Temp_List = []
        i = 0 
        l = []
        while i < 1000:
            [n, F, I] = DefinedCanFunction(j)
            l.append(n[1])
            i += 1
            
        sum_values = sum(l)
        average = sum_values /1000
        Temp_List.append(j)
        Temp_List.append(average)
        Average_List.append(Temp_List)
        j += 1

    # Write Average_List to CSV
    with open('simulation_results.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Index', 'Average'])  # Header row
        writer.writerows(Average_List)  # Data rows

    return Average_List
 
    
 
'''  #Stuff

def RandomCanSimulation(N_1, n_1, N_2, n_2):
    #Initial Variables
    
    #index for the number of trials. 1000 trials per # of connections
    j = 1
    
    #Average List to return at the end of the program
    Average_List = list()

    #Run 20 trials, 1 trial per # of connections
    while j < 21:
        
        
        #Initial Variables for subtrials
        
        #index
        i = 0 
        
        #Used for adding multiple values to Average_Lists
        Temp_List = []
        
        #DIfferent Lists for Percent of Ones in upstream
        list_0_to_10 = []
        list_10_to_20 = []
        list_20_to_30 = []
        list_30_to_40 = []
        list_40_to_50 = []
        list_50_to_60 = []
        list_60_to_70 = []
        list_70_to_80 = []
        list_80_to_90 = []
        list_90_to_100 = []
        
        #Start of 1000 subtrials
        while i < 1000:
            
            #Run RandomCanFunction with two networks of size 6
            [F_1, F_2, I_1, I_2, F_combined, I_combined, Attractors, Percent] = RandomCanFunction(N_1, n_1, N_2, n_2, j, True)
            Percent *= 100
            
            #Sort Attractors into proper list depending on the percent of ones in the upstream
            
            if Percent >= 0 and Percent < 10:
                 list_0_to_10.append(Attractors[2]) 
                 
            elif Percent >= 10 and Percent < 20:
                      list_10_to_20.append(Attractors[2]) 
                      
            elif Percent >= 20 and Percent < 30:
                      list_20_to_30.append(Attractors[2]) 
                      
            elif Percent >= 30 and Percent < 40:
                      list_30_to_40.append(Attractors[2]) 
                      
            elif Percent >= 40 and Percent < 50:
                      list_40_to_50.append(Attractors[2])   
                    
            elif Percent >= 50 and Percent < 60:
                      list_50_to_60.append(Attractors[2])   
                     
            elif Percent >= 60 and Percent < 70:
                      list_60_to_70.append(Attractors[2])    
                      
            elif Percent >= 70 and Percent < 80:
                      list_70_to_80.append(Attractors[2])  
                      
            elif Percent >= 80 and Percent < 90:
                      list_80_to_90.append(Attractors[2])   
                      
            elif Percent >= 90 and Percent < 100:
                      list_90_to_100.append(Attractors[2])    
                      
                      
                      
                      
            #Repeat  
                      
            i += 1
            
        #Take average of each of the percent-lists
        Temp_List.append(j)



        #Average 0-10
        #Find average
        if len(list_0_to_10) > 0:
            sum_values = sum(list_0_to_10)
            average = sum_values / len(list_0_to_10)
            
            #Add average to temp list
            Temp_List.append(average)
        else:
            Temp_List.append(-1)
        
        
        
        #Average 10-20
        if len(list_10_to_20) > 0:
            #Find average
            sum_values = sum(list_10_to_20)
            average = sum_values / len(list_10_to_20)
            
            #Add average to temp list
            Temp_List.append(average)
        else:
            Temp_List.append(-1)
        
        
        
        #Average 20-30
        if len(list_20_to_30) > 0:
            #Find average
            sum_values = sum(list_20_to_30)
            average = sum_values / len(list_20_to_30)
            
            #Add average to temp list
            Temp_List.append(average)
        else:
            Temp_List.append(-1)



        #Average 30-40
        if len(list_30_to_40) > 0:
            #Find average
            sum_values = sum(list_30_to_40)
            average = sum_values / len(list_30_to_40)
            
            #Add average to temp list
            Temp_List.append(average)
        else:
            Temp_List.append(-1)
        
        
        #Average 40-50
        if len(list_40_to_50) > 0:
            #Find average
            sum_values = sum(list_40_to_50)
            average = sum_values / len(list_40_to_50)
            
            #Add average to temp list
            Temp_List.append(average)
        else:
            Temp_List.append(-1)
            
        
        
        #Average 50-60
        if len(list_50_to_60) > 0:
            #Find average
            sum_values = sum(list_50_to_60)
            average = sum_values / len(list_50_to_60)
            
            #Add average to temp list
            Temp_List.append(average)
        else:
            Temp_List.append(-1)
        
        
        
        #Average 60-70
        if len(list_60_to_70) > 0:
            #Find average
            sum_values = sum(list_60_to_70)
            average = sum_values / len(list_60_to_70)
            
            #Add average to temp list
            Temp_List.append(average)
        else:
            Temp_List.append(-1)
        
        
        
        #Average 70-80
        if len(list_70_to_80) > 0:
            #Find average
            sum_values = sum(list_70_to_80)
            average = sum_values / len(list_70_to_80)
            
            #Add average to temp list
            Temp_List.append(average)
        else:
            Temp_List.append(-1)
        
        
        #Average 80-90
        if len(list_80_to_90) > 0:
            #Find average
            sum_values = sum(list_80_to_90)
            average = sum_values / len(list_80_to_90)
            
            #Add average to temp list
            Temp_List.append(average)
        else:
            Temp_List.append(-1)
            
        
        
        #Average 90-100
        if len(list_90_to_100) > 0:
            #Find average
            sum_values = sum(list_90_to_100)
            average = sum_values / len(list_90_to_100)
            
            #Add average to temp list
            Temp_List.append(average)
        else:
            Temp_List.append(-1)
            

        
        Average_List.append(Temp_List)
        j += 1

    # Write Average_List to CSV
    with open('simulation_results.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Connections', '0-10', '10 - 20', '20 - 30', '30 - 40', '40 - 50', '50 - 60', '60 - 70', '70 - 80', '80 - 90', '90 - 100'])  # Header row
        writer.writerows(Average_List)  # Data rows

    return Average_List








def RandomCanSimulation2(N_1, n_1, N_2, n_2):
    # Initial Variables
    max_connections = 2  # Number of trials for connections
    num_trials = 10      # Number of subtrials per connection

    # Dictionary to hold attractor lists per connection
    connections_dict = {}
    
    # Percent range labels
    percent_ranges = {
        '0_to_10': [],
        '10_to_20': [],
        '20_to_30': [],
        '30_to_40': [],
        '40_to_50': [],
        '50_to_60': [],
        '60_to_70': [],
        '70_to_80': [],
        '80_to_90': [],
        '90_to_100': []
    }

    # Run trials for each connection count
    for j in range(1, max_connections + 1):
        # Initialize a new list for each connection count in the dictionary
        connections_dict[j] = {
            '0_to_10': [],
            '10_to_20': [],
            '20_to_30': [],
            '30_to_40': [],
            '40_to_50': [],
            '50_to_60': [],
            '60_to_70': [],
            '70_to_80': [],
            '80_to_90': [],
            '90_to_100': []
        }

        # Start of subtrials
        for i in range(num_trials):
            # Run RandomCanFunction with two networks of size N_1 and N_2
            [F_1, I_1, F_2, I_2, F_combined, I_combined, Attractors, Percent] = RandomCanFunction(N_1, n_1, N_2, n_2, j, True)

            # Convert Percent to a percentage value
            Percent *= 100

            # Append attractors along with F's and I's to the appropriate list based on the percent value
            data_row = {
                'Attractors': Attractors,
                'F_1': F_1, 'I_1': I_1,
                'F_2': F_2, 'I_2': I_2,
                'F_combined': F_combined, 'I_combined': I_combined
            }

            if 0 <= Percent < 10:
                connections_dict[j]['0_to_10'].append(data_row)
            elif 10 <= Percent < 20:
                connections_dict[j]['10_to_20'].append(data_row)
            elif 20 <= Percent < 30:
                connections_dict[j]['20_to_30'].append(data_row)
            elif 30 <= Percent < 40:
                connections_dict[j]['30_to_40'].append(data_row)
            elif 40 <= Percent < 50:
                connections_dict[j]['40_to_50'].append(data_row)
            elif 50 <= Percent < 60:
                connections_dict[j]['50_to_60'].append(data_row)
            elif 60 <= Percent < 70:
                connections_dict[j]['60_to_70'].append(data_row)
            elif 70 <= Percent < 80:
                connections_dict[j]['70_to_80'].append(data_row)
            elif 80 <= Percent < 90:
                connections_dict[j]['80_to_90'].append(data_row)
            elif 90 <= Percent < 100:
                connections_dict[j]['90_to_100'].append(data_row)

    # Write each percent list per connection count to separate CSV files
    for num_connections, percent_lists in connections_dict.items():
        for key, data_list in percent_lists.items():
            if data_list:  # Only save if the list contains data
                # Create a dataframe from the data
                df = pd.DataFrame([{
                    'Attractor 1': data['Attractors'][0],
                    'Attractor 2': data['Attractors'][1],
                    'Attractor 3': data['Attractors'][2],
                    'F_1': data['F_1'], 'I_1': data['I_1'],
                    'F_2': data['F_2'], 'I_2': data['I_2'],
                    'F_combined': data['F_combined'], 'I_combined': data['I_combined']
                } for data in data_list])

                # Save to CSV
                df.to_csv(f'connections_{num_connections}_{key}.csv', index=False)

    return connections_dict


def RandomCanSimulation3(N_1, n_1, N_2, n_2):
    # Initial Variables
    max_connections = 2  # Number of trials for connections
    num_trials = 10    # Number of subtrials per connection

    # Dictionary to hold attractor lists per connection
    connections_dict = {}
    
    # Percent range labels
    percent_ranges = {
        '0_to_10': [],
        '10_to_20': [],
        '20_to_30': [],
        '30_to_40': [],
        '40_to_50': [],
        '50_to_60': [],
        '60_to_70': [],
        '70_to_80': [],
        '80_to_90': [],
        '90_to_100': []
    }

    # Run trials for each connection count
    for j in range(1, max_connections + 1):
        # Initialize a new list for each connection count in the dictionary
        connections_dict[j] = {
            '0_to_10': [],
            '10_to_20': [],
            '20_to_30': [],
            '30_to_40': [],
            '40_to_50': [],
            '50_to_60': [],
            '60_to_70': [],
            '70_to_80': [],
            '80_to_90': [],
            '90_to_100': []
        }

        # Start of subtrials
        for i in range(num_trials):
            # Run RandomCanFunction with two networks of size N_1 and N_2
            [Attractors, Percent] = RandomCanFunction(N_1, n_1, N_2, n_2, j, False)

            # Convert Percent to a percentage value
            Percent *= 100

            # Append attractors to the appropriate list based on the percent value
            if 0 <= Percent < 10:
                connections_dict[j]['0_to_10'].append(Attractors)
            elif 10 <= Percent < 20:
                connections_dict[j]['10_to_20'].append(Attractors)
            elif 20 <= Percent < 30:
                connections_dict[j]['20_to_30'].append(Attractors)
            elif 30 <= Percent < 40:
                connections_dict[j]['30_to_40'].append(Attractors)
            elif 40 <= Percent < 50:
                connections_dict[j]['40_to_50'].append(Attractors)
            elif 50 <= Percent < 60:
                connections_dict[j]['50_to_60'].append(Attractors)
            elif 60 <= Percent < 70:
                connections_dict[j]['60_to_70'].append(Attractors)
            elif 70 <= Percent < 80:
                connections_dict[j]['70_to_80'].append(Attractors)
            elif 80 <= Percent < 90:
                connections_dict[j]['80_to_90'].append(Attractors)
            elif 90 <= Percent < 100:
                connections_dict[j]['90_to_100'].append(Attractors)

    # Write each percent list per connection count to separate CSV files
    for num_connections, percent_lists in connections_dict.items():
        for key, attractor_list in percent_lists.items():
            if attractor_list:  # Only save if the list contains data
                df = pd.DataFrame(attractor_list, columns=['Attractor 1', 'Attractor 2', 'Attractor 3'])
                df.to_csv(f'connections_{num_connections}_{key}.csv', index=False)

    return connections_dict
'''

'''

def RandomCanSimulation4(N_1, n_1, N_2, n_2):
    # Initial Variables
    max_connections = 2  # Number of trials for connections
    num_trials = 10      # Number of subtrials per connection

    # List to store all data for a single CSV file
    all_data = []

    # Percent range labels
    percent_ranges = {
        '0_to_10': [],
        '10_to_20': [],
        '20_to_30': [],
        '30_to_40': [],
        '40_to_50': [],
        '50_to_60': [],
        '60_to_70': [],
        '70_to_80': [],
        '80_to_90': [],
        '90_to_100': []
    }

    # Run trials for each connection count
    for j in range(1, max_connections + 1):
        # Start of subtrials
        for i in range(num_trials):
            # Run RandomCanFunction with two networks of size N_1 and N_2
            [F_1, I_1, F_2, I_2, F_combined, I_combined, Attractors, Percent] = RandomCanFunction(N_1, n_1, N_2, n_2, j, True)

            # Convert Percent to a percentage value
            Percent *= 100

            # Determine the percent range
            if 0 <= Percent < 10:
                percent_range = '0_to_10'
            elif 10 <= Percent < 20:
                percent_range = '10_to_20'
            elif 20 <= Percent < 30:
                percent_range = '20_to_30'
            elif 30 <= Percent < 40:
                percent_range = '30_to_40'
            elif 40 <= Percent < 50:
                percent_range = '40_to_50'
            elif 50 <= Percent < 60:
                percent_range = '50_to_60'
            elif 60 <= Percent < 70:
                percent_range = '60_to_70'
            elif 70 <= Percent < 80:
                percent_range = '70_to_80'
            elif 80 <= Percent < 90:
                percent_range = '80_to_90'
            elif 90 <= Percent < 100:
                percent_range = '90_to_100'
            else:
                continue  # Skip any out-of-range percentages

            # Append the data row to the all_data list
            all_data.append({
                'Percent Range': percent_range,
                'Attractor 1': Attractors[0],
                'Attractor 2': Attractors[1],
                'Attractor 3': Attractors[2],
                'F_1': F_1, 'I_1': I_1,
                'F_2': F_2, 'I_2': I_2,
                'F_combined': F_combined, 'I_combined': I_combined
            })

    # Create a dataframe from all_data and save to a single CSV file
    df = pd.DataFrame(all_data)
    df.to_csv('combined_connections.csv', index=False)

    return df
'''





def RandomCanSimulation(N_1, n_1, N_2, n_2):
    # Initial Variables
    max_connections = 2  # Number of trials for connections
    num_trials = 10      # Number of subtrials per connection

    # List to store all data for a single CSV file
    all_data = []

    # Run trials for each connection count
    for j in range(1, max_connections + 1):
        # Start of subtrials
        for i in range(num_trials):
            # Run RandomCanFunction with two networks of size N_1 and N_2
            [F_1, I_1, F_2, I_2, F_combined, I_combined, Attractors, Percent, Specific_Attractors_F1, Specific_Attractors_F2, Specific_Attractors_F_combined] = RandomCanFunction(N_1, n_1, N_2, n_2, j, True)

            # Convert Percent to a percentage value
            Percent *= 100

            # Append the data row to the all_data list, including the number of connections
            all_data.append({  
                'Connections': j,  # Include the number of connections
                'Percent of 1\'s': Percent,
                'Attractor 1': Attractors[0],
                'Attractor 2': Attractors[1],
                'Attractor 3': Attractors[2],
                'F_1 Attractors': Specific_Attractors_F1,
                'F_2 Attractors': Specific_Attractors_F2,
                'F_combined Attractors': Specific_Attractors_F_combined,
                'F_1': F_1, 'I_1': I_1,
                'F_2': F_2, 'I_2': I_2,
                'F_combined': F_combined, 'I_combined': I_combined
            })

    # Create a dataframe from all_data
    df = pd.DataFrame(all_data)

    # Save to CSV
    csv_filename = 'combined_connections.csv'
    df.to_csv(csv_filename, index=False)

    # Return the DataFrame for further use in Python
    return df

# Example usage:
# df = RandomCanSimulation4(N_1, n_1, N_2, n_2)




'''Special Simulations'''

def AttractorsSimulations(BN_1_Attractors, BN_2_Attractors, Num_of_Connections): #Put -1 if you don't want restrictions on the attractors
    #index
    i = 0 
    k = 0
    TotalAttractors = list()
    
    [List_of_Num_of_Attractors, List_of_F_1, List_of_I_1, List_of_F_2, List_of_I_2, List_of_F_combined, List_of_I_combined] = GenerateCanalyzingNetworks(6, 3, 6, 3)
    [F_1, I_1, F_2, I_2] = SortCanalyzingFunctions(BN_1_Attractors, BN_2_Attractors, List_of_Num_of_Attractors, List_of_F_1, List_of_F_2, List_of_I_1, List_of_I_2)
    
    while i < len(F_1):
        I_1[i] = ShiftNetworkI(len(I_2[i]), I_1[i], 0, len(I_1[i]) - 1) #Shift I_1
        F_combined = JoinTwoNetworks(F_2[i], F_1[i], 0) #Join the F's
        I_combined = RandomlyConnectNetworks(I_1[i], I_2[i], Num_of_Connections)  #Connections flow from I_1 to I_2
        UpdateF(F_combined, I_combined) #Update F to be canalyzing
        k = can2.num_of_attractors(F_combined, I_combined, len(I_combined))
        TotalAttractors.append(k[1])
        i += 1 
    
    return sum(TotalAttractors)/len(F_1)
    
    
'''Visualizing Data'''

def CreateLineGraph(df, x, y):

    # Assuming df is your DataFrame and you want to plot 'Percent of 1's' against 'Connections'
    plt.figure(figsize=(8, 6))
    
    # Plot the data
    plt.plot(df[x], df[y], marker='o', linestyle='-', color='b')
    
    # Add labels and title
    plt.xlabel(x)
    plt.ylabel(y)
    plt.title(y + ' vs ' +  x)
    
    # Show the plot
    plt.show()
    
def CreateBarGraph(df):
    # Group by number of connections and calculate the mean percent of 1's
    mean_percent = df.groupby('Connections')['Percent of 1\'s'].mean()

    # Create the bar graph
    plt.bar(mean_percent.index, mean_percent.values, color='skyblue')

    # Add title and labels
    plt.title('Mean Percent of 1\'s by Number of Connections')
    plt.xlabel('Number of Connections')
    plt.ylabel('Mean Percent of 1\'s')

    # Show the plot
    plt.xticks(mean_percent.index)  # Show ticks for each connection
    plt.grid(axis='y', linestyle='--', alpha=0.7)  # Add grid lines for better visibility
    plt.show()

    

#x =  'Number of Connections'
#y = 'Percent of 1\'s'
