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
import json
import ast



'''Modify Networks'''

def PercentOfOnes(F, I):
    #determine specific attractors
    Attractors = can2.num_of_attractors(F, I, len(I), EXACT = True)[0]
    
    #Put all attractors into 1 list
    Temp_Variable = list()
    for index in range(len(Attractors)):
        Temp_Variable += Attractors[index]
    
    Attractors = Temp_Variable
    
    #Define total binary digits and 1's
    TotalDigits = 0
    TotalOnes = 0
    
    #For each attractor
    for index in range(len(Attractors)): 

           
        #Convert attractor into binary number
        binList = can13.dec2bin(Attractors[index]) 
        
        
        #Iterate throughh binList (binary number)
        for index_1 in range(len(binList)): 
            
             #determine number of ones and 0's
             if binList[index_1] != 0: 
                 TotalOnes += 1

        #add length of I (max length of binary num, ensures extra 0s are counted)
        TotalDigits += len(I) 
           
    return TotalOnes/TotalDigits #returns decimal


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

    
    
  
def RandomlyConnectNetworks(i_1, i_2, connections): #Connections flow from I_1 to I_2
    #Make sure original I_1 and I_2 are not changed
    I_1 = i_1.copy()
    I_2 = i_2.copy() 
    
    for index in range(len(I_1)): 
        I_1[index] += (len(I_2))*(np.ones(len(I_1[index]), dtype = int)) #Shift each value in network[index] by N

    




    #Define list to return
    I_1List = list()
    
    for index in range(len(I_1)):
        I_1List += list(I_1[index]) #Create a list with ever value in I (possible duplicates)
        
    #Remove duplicates from list and put values in order
    #Set removes duplicates, list converts to list, sorted orders list
    I_1List = sorted(list(set(I_1List)))
    
    
    
    #Define list to return
    I_2List = list()
    
    for index in range(len(I_2)):
        I_2List += list(I_2[index]) #Create a list with ever value in I (possible duplicates)
        
    #Remove duplicates from list and put values in order
    #Set removes duplicates, list converts to list, sorted orders list
    I_2List = sorted(list(set(I_2List)))


        
    #Iterate 
    index = 0
    while index < connections:
        
        I_1_Random_Node = random.sample(I_1List, 1) #Pick a node from I_1 to add to I_2
        I_2_Random_Node = random.sample(I_2List, 1)[0] #Pick a variable from I_2

        
        if (list(I_2[I_2_Random_Node].copy())).count(I_1_Random_Node[0]) < 1: #Checks to see if there is already a connection 
        
        
            I_2[I_2_Random_Node] = np.concatenate((I_2[I_2_Random_Node], I_1_Random_Node), axis = None) #Connect those nodes to a specific node in I_2

            index += 1
            
        else:
            index = index #Don't progress, repeat
        
    return I_2 + I_1

def UpdateF(F, I):  #Designed for canalyzing functions, F_combined in the later part of the code
    for nodeInF in range(len(I)): #Iterate through every value in F
        if len(F[nodeInF]) < 2**len(I[nodeInF]): #Checks to see if F has enough 1's

            #Iterate enough times to add enough 1's to F
            for index_1 in range(2**len(I[nodeInF]) - len(F[nodeInF])):
                     F[nodeInF] = np.concatenate((F[nodeInF], [1]), axis = None) #add a bunch of ones



def CombineTwoNetworks(f_1, i_1, f_2, i_2, connections): #Uses Preset F_1 and F_2, connections flow from F_1 to F_2, 

    #Make copies of networks so that originals are not altered
    [F_1, I_1, F_2, I_2] = [f_1.copy(), i_1.copy(), f_2.copy(), i_2.copy()]
    
    #Join F-1 and F-2 and join I_1 and I_2
    [F_combined, I_combined] = [F_2 + F_1, RandomlyConnectNetworks(I_1, I_2, connections)]
    
    #Correct the number of 1's to make F_combined canalyzing
    UpdateF(F_combined, I_combined)
    

    return F_combined, I_combined, can2.num_of_attractors(F_combined, I_combined, len(F_combined), EXACT = True) #Return network and attractors


def RandomCanFunction(N_1, n_1, N_2, n_2, Num_of_Connections, specific):
    #Intital Lists
    Num_of_Attractors = list()
    Specific_Attractors = list()
    
    #Percennt of Ones List
    BN_1_Percent_Of_Ones = 0

    
    #Define BN-1 and BN-2
    [F_1, I_1, degree_1] = can13.random_BN(N = N_1, n = n_1) 
    [F_2, I_2, degree_2] = can13.random_BN(N = N_2 , n = n_2)
    
    #Calculate attractors of BN-1 and BN-2
    [BN_1_Attractors, BN_2_Attractors] = [can2.num_of_attractors(F_1, I_1, len(F_1), EXACT = True), can2.num_of_attractors(F_2, I_2, len(F_2), EXACT = True)]
    
    #Determine Percent of Ones
    BN_1_Percent_Of_Ones = PercentOfOnes(F_1, I_1)
    
        
    #Make combined network aand determine attractors
    [F_combined, I_combined, BN_Combined_Attractors] = CombineTwoNetworks(F_1, I_1, F_2, I_2, Num_of_Connections)
    
    
    #Add attractors to lists
    Specific_Attractors.append(BN_1_Attractors[0])
    Specific_Attractors.append(BN_2_Attractors[0])   
    Specific_Attractors.append(BN_Combined_Attractors[0])
    
    Num_of_Attractors.append(BN_1_Attractors[1])
    Num_of_Attractors.append(BN_2_Attractors[1])
    Num_of_Attractors.append(BN_Combined_Attractors[1])    


        

    #Return

    if specific == True:
        return F_1, F_2, I_1, I_2, F_combined, I_combined, Num_of_Attractors, Specific_Attractors, BN_1_Percent_Of_Ones
    else:
        return Num_of_Attractors, BN_1_Percent_Of_Ones





#Rewrite
def Experiment_1_Generate_Networks(Num_Of_Networks): 
    all_data = []
    for num in range(Num_Of_Networks):
        
        #Random N (Nodes) and n (outgoing connections from each node) for networks
        N = random.randint(3, 5)
        n = random.randint(1, N-1)
        
        #Generate Networks
        [F, I, degree] = can13.random_BN(N, n)
        
        #Make copies, ensures orignials are not changed
        [f, i] = [F.copy(), I.copy()]
        
        #Determine Percent of Ones of Upstream Network
        Percent = PercentOfOnes(f, i)
        Percent *= 100
        
        #Determine Attractors of Each Network
        Attractors = can2.num_of_attractors(f, i, len(f), EXACT = True)[1]

        
        #Append Data
        all_data.append({  
            'Percent of 1\'s': Percent,
            'Attractors': Attractors,
            'F': list(F), 'I': list(I), 

        })
        
    # Create a dataframe from all_data
    df = pd.DataFrame(all_data)

    # Save to CSV
    csv_filename = 'Experiment 1 Networks.csv'
    df.to_csv(csv_filename, index=False)

    # Return the DataFrame for further use in Python
    return df


def Connect_Networks(df, Num_of_Networks):
    #Initialize Dataframe
    all_data = []

    
    #repeat Num_of_Networks times
    for index in range(Num_of_Networks):
        Network1 = random.randint(0, 10)
        Network2 = random.randint(0, 10)
        
        
        # Generate the networks
        [F_1, I_1, F_2, I_2] = [
            df.loc[Network1, 'F'],
            df.loc[Network1, 'I'],
            df.loc[Network2, 'F'],
            df.loc[Network2, 'I']
        ]
        
        # Convert from string to list
        [F_1_processed, I_1_processed, F_2_processed, I_2_processed] = [
            F_1.replace("array(", "").replace(")", ""),
            I_1.replace("array(", "").replace(")", ""),
            F_2.replace("array(", "").replace(")", ""),
            I_2.replace("array(", "").replace(")", "")
        ]
        
        [F_1, I_1, F_2, I_2] = [
            ast.literal_eval(F_1_processed),
            ast.literal_eval(I_1_processed),
            ast.literal_eval(F_2_processed),
            ast.literal_eval(I_2_processed)
        ]


        
        for connections in range(1, len(F_1) * len(F_2)): #Repeat for connections between min and max connections (1 and N^2. len(F_1) = N)
        
                #Make copies, ensures orignials are not changed
                [f_1, f_2, i_1, i_2] = [F_1.copy(), F_2.copy(), I_1.copy(), I_2.copy()]
                
                #Connect Networks
                [F_combined, I_combined, BN_Combined_Attractors] = CombineTwoNetworks(f_1, i_1, f_2, i_2, connections)
               
                all_data.append({  
                    'Connections': connections,  # Include the number of connections
                    'BN Combined_Attractors': BN_Combined_Attractors[1], #Attractors
                    'F_combined': F_combined, 'I_combined': I_combined #Specific F and I combined
                })
                

    # Create a dataframe from all_data
    df1 = pd.DataFrame(all_data)

    # Save to CSV
    csv_filename = 'Experiment 1 Combined Networks.csv'
    df1.to_csv(csv_filename, index=False)

    # Return the dataFrame for further use in Python
    return df1
            

def Experiment2():
    all_data = []
    for i in range(1000):
        N1 = random.randint(3, 6)
        n1 = random.randint(1, N1-1)
        N2 = random.randint(3, 6)
        n2 = random.randint(1, N2-1)
        [F_1, I_1, degree] = can13.random_BN(N1, n1)
        [F_2, I_2, degree] = can13.random_BN(N2, n2)
        print(i)
        for connections in range(1, (N1 * N2)):
            f_1 = F_1.copy()
            f_2 = F_2.copy()
            i_1 = I_1.copy()
            i_2 = I_2.copy()
            Percent = PercentOfOnes(F_1, I_1)
            Percent *= 100
            
            BN1_Attractors = can2.num_of_attractors(F_1, I_1, len(F_1), EXACT = True)
            BN2_Attractors = can2.num_of_attractors(F_1, I_1, len(F_1), EXACT = True)
            
            BN1_Attractors = BN1_Attractors[1]
            BN2_Attractors = BN2_Attractors[1]
            
            for j in range(1000):
                [F_combined, I_combined, BN_Combined_Attractors] = CombineTwoNetworks(f_1, i_1, f_2, i_2, connections)
                
                


    
                # Append the data row to the all_data list, including the number of connections
                all_data.append({  
                    'Connections': connections,  # Include the number of connections
                    'Percent of 1\'s': Percent,
                    'BN1 Attractors': BN1_Attractors,
                    'BN2 Attractors': BN2_Attractors,
                    'BN Combined_Attractors': BN_Combined_Attractors[1],
                    'F_1': F_1, 'I_1': I_1,
                    'F_2': F_2, 'I_2': I_2,
                    'F_combined': F_combined, 'I_combined': I_combined
                })
    
    # Create a dataframe from all_data
    df = pd.DataFrame(all_data)

    # Save to CSV
    csv_filename = 'Experiment 2.csv'
    df.to_csv(csv_filename, index=False)

    # Return the DataFrame for further use in Python
    return df


def RandomCanSimulation(N_1, n_1, N_2, n_2):
    # Initial Variables
    max_connections = 20  # Number of trials for connections
    num_trials = 1000      # Number of subtrials per connection

    # List to store all data for a single CSV file
    all_data = []

    # Run trials for each connection count
    for j in range(1, max_connections + 1):
        # Start of subtrials
        for i in range(num_trials):
            # Run RandomCanFunction with two networks of size N_1 and N_2
            [F_1, F_2, I_1, I_2, F_combined, I_combined, Num_of_Attractors, Specific_Attractors, BN_1_Percent_Of_Ones] = RandomCanFunction(N_1, n_1, N_2, n_2, j, True)

            # Convert Percent to a percentage value
            BN_1_Percent_Of_Ones *= 100

            # Append the data row to the all_data list, including the number of connections
            all_data.append({  
                'Connections': j,  # Include the number of connections
                'Percent of 1\'s': BN_1_Percent_Of_Ones,
                'Num of BN_1 Attractors': Num_of_Attractors[0],
                'Num of BN_2 Attractors': Num_of_Attractors[1],
                'Num of BN_Combined Attractors': Num_of_Attractors[2],
                'BN_1 Attractors': Specific_Attractors[0],
                'BN_2 Attractors': Specific_Attractors[1],
                'BN_Combined Attractors': Specific_Attractors[2],
                'F_1': F_1, 'I_1': I_1,
                'F_2': F_2, 'I_2': I_2,
                'F_combined': F_combined, 'I_combined': I_combined
            })

    # Create a dataframe from all_data
    df = pd.DataFrame(all_data)

    # Save to CSV
    csv_filename = 'New_combined_connections_17.csv'
    df.to_csv(csv_filename, index=False)

    # Return the DataFrame for further use in Python
    return df







    
'''Visualizing Data'''











def CreateLineGraph(df, df_x, df_y, x_values, y_values, x_label, y_label, plot_title):
    plt.figure(figsize=(8, 6))
    
    # Plot the data
    plt.plot(x_values, y_values, marker='o', linestyle='-', color='b', label='Data')
    
    # Fit a line to the data
    coefficients = np.polyfit(x_values, y_values, 1)
    fitted_line = np.polyval(coefficients, x_values)
    plt.plot(x_values, fitted_line, 'r--', label='Best fit line')

    # Calculate R² value manually
    y_mean = np.mean(y_values)
    ss_total = np.sum((y_values - y_mean) ** 2)  # Total sum of squares
    ss_residual = np.sum((y_values - fitted_line) ** 2)  # Residual sum of squares
    r_squared = 1 - (ss_residual / ss_total)
    
    
    #Calculate the Correlation Coefficient
    x = df[df_x]
    y = df[df_y]
    
    correlation_matrix = np.corrcoef(x, y)
    correlation_coefficient = correlation_matrix[0, 1]

    # Display R² value on the plot
    plt.annotate(f'R² = {r_squared:.4f}', 
                 xy=(0.70, 0.10),  # Adjusted position to bottom-right
                 xycoords='axes fraction', 
                 fontsize=14, 
                 color='black', 
                 bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))
    
    #Display the Corrrelation Coeefficient on the plot
    plt.annotate(f'r = {correlation_coefficient:.4f}', 
                 xy=(0.70, 0.18),  # Position it above R²
                 xycoords='axes fraction', 
                 fontsize=14, 
                 color='black', 
                 bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))

    # Add labels and title
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(plot_title)
    plt.legend(loc='upper left')  # Move legend to upper left to avoid overlap
    
    # Save and show the plot
    plt.savefig(plot_title, dpi=300, bbox_inches='tight')
    plt.show()

# Example usage
# CreateLineGraph([1, 2, 3, 4], [10, 20, 15, 25], 'Connections', 'Percent of 1\'s', restrict_x_int=True)


    

def CreateBarGraph(df, x, y, plot_title, xlabel, ylabel, connections=None, operation=None, min_percent=None, max_percent=None, dummy_x = None):
    """
    Creates a bar graph from the given DataFrame with optional filtering based on connections
    and a 'Percent_of_1\'s' range.
    
    Parameters:
        df (pd.DataFrame): The input DataFrame.
        x (str): The column for the x-axis.
        y (str): The column for the y-axis.
        connections (int): Filter for 'Connections' column (0 for no filtering).
        operation (str): Operation to perform ('subtract' or 'divide').
        min_percent (float): Minimum percent to filter 'Percent of 1\'s'.
        max_percent (float): Maximum percent to filter 'Percent of 1\'s'.
        
    Returns:
        list: The computed mean values as a list.
    """
    
    #Sorts Data based on initial conditions
    # Filter by 'Connections' if needed
    if connections is not None:
        df_filtered = df.copy()
        df_filtered = df_filtered[df_filtered['Connections'] == connections]
    else:
        df_filtered = df.copy()


    # Filter by 'Percent_of_1\'s' if min_percent and max_percent are specified
    if min_percent is not None and max_percent is not None:
        df_filtered = df_filtered[
            (df_filtered['Percent of 1\'s'] >= min_percent) &
            (df_filtered['Percent of 1\'s'] <= max_percent)
        ]

    # If no data remains after filtering, return early
    if df_filtered.empty:
        print("No data points fall within the specified filters.")
        return []



    #Change operation between subtract, divide, and noormal
    if operation == "subtract":
        # Subtract x from y, take the average
        df_filtered["Difference"] = df_filtered[y] - df_filtered[dummy_x]
        mean_result = df_filtered.groupby("Connections")["Difference"].mean()
        
    elif operation == "divide":
        # Divide y by x, take the average
        df_filtered["Ratio"] = df_filtered[y]/df_filtered[dummy_x]
        mean_result = df_filtered.groupby("Connections")["Ratio"].mean()

    else:
        # Default behavior for mean of y grouped by x
        mean_result = df_filtered.groupby(x)[y].mean()        


    # Create the bar graph
    plt.bar(mean_result.index, mean_result.values, color='skyblue')
    plt.title(plot_title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(plot_title, dpi=300, bbox_inches='tight')



    # Show the plot
    plt.xticks(mean_result.index)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.show()

    return mean_result.tolist()

               



def CreateViolinPlot(df, x, y, connections=None, operation=None, dummy_x=None, plot_title=None):
  """
  Creates a violin plot from the given DataFrame with optional filtering and operations.
  
  Parameters:
      df (pd.DataFrame): The input DataFrame.
      x (str): The column for the x-axis.
      y (str): The column for the y-axis.
      connections (int): Filter for 'Connections' column (0 for no filtering).
      operation (str): Operation to perform ('subtract' or 'divide').
      dummy_x (str): The column to use in the operation (required for 'subtract' or 'divide').
      plot_title (str): Title for the plot.
      
  Returns:
      None
  """
  df_filtered = df.copy()

  # Filter by 'Connections' if specified
  if connections is not None:
      df_filtered = df_filtered[df_filtered['Connections'] == connections]

  # If no data remains after filtering, return early
  if df_filtered.empty:
      print("No data points fall within the specified filters.")
      return

  # Perform the specified operation
  if operation == "subtract":
      if dummy_x is None:
          print("dummy_x must be specified for 'subtract' operation.")
          return
      df_filtered["Difference"] = df_filtered[y] - df_filtered[dummy_x]
      plot_data = df_filtered.groupby(x)["Difference"].apply(list)

  elif operation == "divide":
      if dummy_x is None:
          print("dummy_x must be specified for 'divide' operation.")
          return
      df_filtered["Ratio"] = df_filtered[y] / df_filtered[dummy_x]
      plot_data = df_filtered.groupby(x)["Ratio"].apply(list)

  else:
      # Default behavior: use y values grouped by x
      plot_data = df_filtered.groupby(x)[y].apply(list)

  # Create the violin plot
  plt.figure(figsize=(10, 6))
  plt.violinplot(plot_data, showmeans=True, showmedians=True)

  # Set x-ticks and labels
  x_categories = plot_data.index
  plt.xticks(ticks=np.arange(1, len(x_categories) + 1), labels=x_categories)
  plt.title(plot_title if plot_title else f'Violin Plot of {y} vs {x}')
  plt.xlabel(x)
  plt.ylabel(operation.capitalize() if operation else y)
  plt.grid(axis='y', linestyle='--', alpha=0.7)
  plt.savefig(plot_title if plot_title else f'Violin Plot of {y} vs {x}', dpi=300, bbox_inches='tight')
  plt.show()
  
    

def CountPercentIntervalsList(df, percent_column='Percent of 1\'s', bin_size=10):
    """
    Counts the number of items within percent intervals and returns a list.

    Parameters:
        df (pd.DataFrame): The DataFrame containing the data.
        percent_column (str): The column name for percent values.
        bin_size (int): The size of each interval (e.g., 10 for 0-10%, 10-20%).

    Returns:
        list: A list of counts for each interval.
    """
    # Define the bin edges from 0% to 100% with the specified bin size
    bin_edges = list(range(0, 101, bin_size))
    
    # Use pd.cut to create bins and include lowest boundary in the first interval
    counts = pd.cut(df[percent_column], bins=bin_edges, right=True, include_lowest=True).value_counts().sort_index()
    
    # Return the counts as a list
    return counts.tolist()

def count_networks_by_attractors(df, attractor):
    """
    Counts the number of networks with 1 attractor, 2 attractors, etc., for F1.

    Parameters:
        df (pd.DataFrame): The DataFrame containing the data. It should have a column named 'F_1 Attractors'.

    Returns:
        list: A list where the i-th element represents the count of networks with (i+1) attractors.
    """
    # Ensure the column is numeric
    df[attractor] = pd.to_numeric(df[attractor], errors='coerce')
    df = df.dropna(subset=[attractor])
    df[attractor] = df[attractor].astype(int)

    # Get the maximum number of attractors in the data
    max_attractors = df[attractor].max()
    
    # Initialize a list to store counts for each number of attractors
    counts = [0] * max_attractors

    # Iterate through each row and count occurrences
    for num in df[attractor]:
        counts[num - 1] += 1

    return counts

df = pd.read_csv('AllData.csv')

#Generates the average y values on the Line graph
y_values = CreateBarGraph(df, 'Attractor 2', 'Attractor 3', 'Mean Total Attractors vs Downstream Attractors', 'Downstream Attractors', 'Mean Total Attractors')

#Remove values after 6th value of y_valuess
i = 0
Temp = list()
while i < 6:
    Temp.append(y_values[i])
    i += 1

#Make x_values

x_values = [1, 2, 3, 4, 5, 6]

#Make line graph

CreateLineGraph(df, 'Attractor 2', 'Attractor 3', x_values, y_values, 'Downstream Attractors', 'Mean Total Attractors',  'Mean Total Attractors vs Downstream Attractors')

#pd.read_csv()





