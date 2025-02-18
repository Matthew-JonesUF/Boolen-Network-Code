#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 12 15:53:40 2019

@author: ckadelka

#Version 1.9

This toolbox is for the investigation of Boolean functions and Boolean networks and
particular targeted at canalization. Throughout, a Boolean function in n variables is
stored as an array of length 2^n (the condensed truth table, i.e. the right side).

0) some basic useful tricks
a) convert a decimal number to a binary number in array form
b) convert a binary number in array form into a decimal number
c) find the indices of all matches in an arrray (list or np.array) 
d) convert an edgelist [(a,b),(c,d),...] of a graph into a list of inputs per node
e) generate the condensed truth table of a Boolean function specified as text in publications, etc.

1) determine if a Boolean function is 
a) degenerated
b) constant
c) canalizing
d) k-canalizing
e) k-set canalizing (collectively canalizing)

determine for a Boolean function,
f) the number of essential variables
g) the symmetry groups (x and y are symmetric if f(x,y) = f(y,x), irrespective of other inputs)
h) the k-set canalizing proportion (the probability that a random k-set canalizes the function)
i) all canalizing variables
j) absolute bias of a Boolean function
k) average sensitivity of a Boolean function

2) determine the 
a) canalizing depth of a Boolean function 
b) layer structure of a canaliizing Boolean function

3) randomly generate (uniform distribution)
a) non-degenerated Boolean functions
a) non-canalizing Boolean functions
c) non-canalizing non-degenerated Boolean functions
d) k-canalizing Boolean functions 
e) k-canalizing Boolean functions with a given layer structure
f) Boolean functions with exact canalizing depth k
g) Boolean functions with exact canalizing depth k with a given layer structure
h) a Boolean network with specified in-degree distribution, canalizing depth, etc.

4) obtain some basic estimates of the magnitude of various subclasses of Boolean functions

5) determine for a Boolean network, 
a) the constant nodes (external parameters)
b) the next state
c) the attractors (exact or approximative)
d) the Derrida value
e) all feed-forward loops and their type
f) all feed-back loops up to size m and their type

6) some plots
a) function that generates networkx.Graph from an edgelist (a list of inputs per node), plot this using nx.draw(G)
"""

#1.9: new functionality added: calculate feed forward loops and feedback loops
#1.5: added is_collectively_canalizing
#1.4: fixed issue with k==0 and EXACT_DEPTH==True in random_BN
#1.3: Python3.7 compatible, kis passed to random_BN can also be [0] for random networks
#1.2: added functionality to randomly create and analyze Boolean networks based on random or canalizing functions
#1.1: fixed a couple issues in is_k_canalizing, is_k_canalizing_return_inputs_outputs_corefunction and get_layer_structure_given_outputs_corefunction

##Imports

import numpy as np
import matplotlib.pyplot as plt
import itertools
import networkx as nx

import sympy
import pandas as pd

from collections import defaultdict
from matplotlib import colors as mcolors

## 0) Basics

def tobin(x):
    '''returns the binary representation (in array form) of a decimal number'''
    return tobin(x//2) + [x%2] if x > 1 else [x]

def dec2bin(x,n=[]):
    '''returns the binary representation (in array form) of a decimal number.
    Input can be an array itself, in which case each decimal number in the 
    array is separately converted into a binary array.
    The second input n, if specified, describes the number of positions in the
    binary representation array for each number.
    Example: dec2bin(10)=dec2bin(10,4)=[1,0,1,0]
             dec2bin(10,6)=[0,0,1,0,1,0]'''
    if type(x) in [list,np.ndarray]:
        return [dec2bin(el,n) for el in x]
    if n==[]:
        return tobin(x)
    else:
        help=tobin(x)
        res=[0]*(n-len(help))
        res.extend(help)
        return res

def bin2dec(state):
    n = len(state)
    b = [2**i for i in range(n)]
    return sum([state[n-i-1]*b[i] for i in range(n)])

def dim(A):
    try:
        return (np.size(A,0),np.size(A,1),np.size(A,2),np.size(A,3),np.size(A,4))
    except IndexError:
        try:
            return (np.size(A,0),np.size(A,1),np.size(A,2),np.size(A,3))
        except IndexError:
            try:
                return (np.size(A,0),np.size(A,1),np.size(A,2))
            except IndexError:
                try:
                    return (np.size(A,0),np.size(A,1))
                except IndexError:
                    try:
                        return (np.size(A,0))
                    except IndexError: 
                        return (np.size(A))

def find_all_indices(array,el):
    res=[]
    for i,a in enumerate(array):
        if a==el:
            res.append(i)
    if res==[]:
        raise ValueError('The element is not in the array at all')
    return res

def edgelist_to_I(edgelist):
    regulators = np.array(edgelist)[:,0]
    targets = np.array(edgelist)[:,1]
    var = list(set(regulators)|set(targets))
    n_var = len(var)
    dict_var = dict(zip(var,range(n_var)))
    I = [[] for _ in range(n_var)]
    for i in range(len(regulators)):
        I[dict_var[targets[i]]].append(dict_var[regulators[i]])
    return I,var

def f_from_expression(expr):
    expr = expr.replace('(',' ( ').replace(')',' ) ')
    expr_split = expr.split(' ')
    var = []
    dict_var = dict()
    n_var = 0
    for i,el in enumerate(expr_split):
        if el not in ['',' ','(',')','and','or','not','AND','OR','NOT','&','|','~','+','-','*','%','>','>=','==','<=','<'] and not el.isdigit():
            try:
                new_var = dict_var[el]
            except KeyError:
                new_var = 'x[%i]' % n_var
                dict_var.update({el:new_var})
                var.append(el)
                n_var += 1
            expr_split[i] = new_var
    expr = ' '.join(expr_split)
    F = []
    for x in itertools.product([0, 1], repeat = n_var):
        F.append(int(eval(expr))) #x is used here "implicitly"
    return F,var

## 1) Methods to analyze Boolean functions

def is_degenerated(F):
    len_F = len(F)
    n=int(np.log2(len_F))
    for i in range(n):
        dummy_add=(2**(n-1-i))
        dummy=np.arange(2**n)%(2**(n-i))//dummy_add
        depends_on_i=False
        for j in range(2**n):
            if dummy[j]==1:
                continue
            else:
                if F[j]!=F[j+dummy_add]:
                    depends_on_i=True
                    break
        if depends_on_i==False:
            return True
    return False

def nr_essential_variables(F):
    len_F = len(F)
    n=int(np.log2(len_F))
    nr_essential_variables = n
    for i in range(n):
        dummy_add=(2**(n-1-i))
        dummy=np.arange(2**n)%(2**(n-i))//dummy_add
        depends_on_i=False
        for j in range(2**n):
            if dummy[j]==1:
                continue
            else:
                if F[j]!=F[j+dummy_add]:
                    depends_on_i=True
                    break
        if depends_on_i==False:
            nr_essential_variables -= 1
    return nr_essential_variables

def is_constant(F):
    return sum(F) in [0,len(F)]

def is_monotonic(F,GET_DETAILS=False):
    n=int(np.log2(len(F)))
    F = np.array(F)
    monotonic = []
    for i in range(n):
        dummy_add=(2**(n-1-i))
        dummy=np.arange(2**n)%(2**(n-i))//dummy_add
        diff = F[dummy==1]-F[dummy==0]
        min_diff = min(diff)
        max_diff = max(diff)
        if min_diff==0 and max_diff==0:
            monotonic.append('not essential')
        elif min_diff==-1 and max_diff==1:
            monotonic.append('not monotonic')
        elif min_diff>=0 and max_diff==1:
            monotonic.append('increasing')            
        elif min_diff==-1 and max_diff<=0:
            monotonic.append('decreasing')   
    if GET_DETAILS:
        return ('not essential' not in monotonic,monotonic)
    else:
        return 'not essential' not in monotonic

def get_symmetry_groups(F,bool_list=[]):
    len_F = len(F)
    n=int(np.log2(len_F))
    if bool_list==[] or bool_list.shape[0]!=len_F:
        bool_list = np.array(list(itertools.product([0, 1], repeat=n)))
    symmetry_groups = []
    left_to_check = np.ones(n)
    for i in range(n):
        if left_to_check[i]==0:
            continue
        else:
            symmetry_groups.append([i])
            left_to_check[i]=0
        for j in range(i+1,n):
            diff = sum(2**np.arange(n-i-2,n-j-2,-1))
            for ii,x in enumerate(bool_list):
                if x[i]!=x[j] and x[i]==0 and F[ii]!=F[ii+diff]:
                    break
            else:
                left_to_check[j] = 0
                symmetry_groups[-1].append(j)
    return symmetry_groups

def is_canalizing(F,n):
    if type(F) == list:
        F = np.array(F)
    desired_value = 2**(n-1)
    T = np.array(list(itertools.product([0, 1], repeat=n))).T
    A = np.r_[T,1-T]
    AtimesF = np.dot(A,F)
    if np.any(AtimesF==desired_value):
        return True
    elif np.any(AtimesF==0):
        return True
    else:
        return False

def is_collectively_canalizing(F,k,n):
    '''computationally poor check'''
    if type(F) == list:
        F = np.array(F)
    if k==0:
        return is_constant(F)
    desired_value = 2**(n-k)
    T = np.array(list(itertools.product([0, 1], repeat=n))).T
    A = np.r_[T,1-T]
    Ak = []
    for i in range(2*n):
        for j in range(i+1,2*n):
            if j-i == n:
                continue
            else:
                Ak.append( np.bitwise_and(A[i,:],A[j,:]) )
                
    Ak = []                
    for indices in itertools.combinations(range(2*n),k):
        dummy = np.sum(A[np.array(indices),:],0)==k
        if sum(dummy)==desired_value:
            Ak.append(dummy)
    Ak = np.array(Ak)
    
    AktimesF = np.dot(Ak,F)
    is_there_canalization = 0 in AktimesF or desired_value in AktimesF
    return is_there_canalization
    
def get_proportion_of_collectively_canalizing_input_sets(F,k,n,bool_list=[]):
    if type(F) == list:
        F = np.array(F)
    if k==0:
        return float(is_constant(F))
    desired_value = 2**(n-k)
    if bool_list == []:
        T = np.array(list(itertools.product([0, 1], repeat=n))).T
    else:
        T = np.array(bool_list).T        
    A = np.r_[T,1-T]
    Ak = []                
    for indices in itertools.combinations(range(2*n),k):
        dummy = np.sum(A[np.array(indices),:],0)==k
        if sum(dummy)==desired_value:
            Ak.append(dummy)
    Ak = np.array(Ak)
    is_there_canalization = np.in1d(np.dot(Ak,F),[0,desired_value])    
    return sum(is_there_canalization)/len(is_there_canalization)

def get_canalizing_strength(F,bool_list=[]):
    nfloat = np.log2(len(F))
    n = int(nfloat)
    assert abs(n-nfloat)<1e-10, "F needs to be of length 2^n for some n>1"
    assert n>1, "Canalizing strength is only defined for Boolean functions with n>1 inputs"
    res = []
    for k in range(1,n):
        res.append(get_proportion_of_collectively_canalizing_input_sets(F,k,n,bool_list=bool_list))
    return np.mean(np.multiply(res,2**np.arange(1,n)/(2**np.arange(1,n)-1))),res

def is_k_canalizing(F,k,n):
    if k>n:
        return False
    if k==0:
        return True
    w = sum(F)
    if w == 0 or w == 2**n: #constant F
        return False
    if type(F) == list:
        F = np.array(F)
    desired_value = 2**(n-1)
    T = np.array(list(itertools.product([0, 1], repeat=n))).T
    A = np.r_[T,1-T]
    try: #is 1 canalized output for one of the variables
        index = list(np.dot(A,F)).index(desired_value)
        newF = F[np.where(A[index]==0)[0]]
        return is_k_canalizing(newF,k-1,n-1)
    except ValueError:
        try: #is 0 canalized output for one of the variables
            index = list(np.dot(A,1-F)).index(desired_value)
            newF = F[np.where(A[index]==0)[0]]
            return is_k_canalizing(newF,k-1,n-1)
        except ValueError:
            return False

def is_k_canalizing_return_inputs_outputs_corefunction(F,k,n,can_inputs=np.array([],dtype=int),can_outputs=np.array([],dtype=int)):
    if k==0:
        return (True,can_inputs,can_outputs,F)
    w = sum(F)
    if w == 0 or w == 2**n: #constant F
        return (False,can_inputs,can_outputs,F)
    if type(F) == list:
        F = np.array(F)
    desired_value = 2**(n-1)
    T = np.array(list(itertools.product([0, 1], repeat=n))).T
    A = np.r_[T,1-T]
    try: #is 1 canalized output for one of the variables
        index = list(np.dot(A,F)).index(desired_value)
        newF = F[np.where(A[index]==0)[0]]
        return is_k_canalizing_return_inputs_outputs_corefunction(newF,k-1,n-1,np.append(can_inputs,int(index<n)),np.append(can_outputs,1))
    except ValueError:
        try: #is 0 canalized output for one of the variables
            index = list(np.dot(A,1-F)).index(desired_value)
            newF = F[np.where(A[index]==0)[0]]
            return is_k_canalizing_return_inputs_outputs_corefunction(newF,k-1,n-1,np.append(can_inputs,int(index<n)),np.append(can_outputs,0))
        except ValueError:
            return (False,can_inputs,can_outputs,F)     

def is_k_canalizing_return_inputs_outputs_corefunction_order(F,k,n,can_inputs=np.array([],dtype=int),can_outputs=np.array([],dtype=int),can_order=np.array([],dtype=int),variables=[]):
    if k==0:
        return (True,can_inputs,can_outputs,F,can_order)
    w = sum(F)
    if w == 0 or w == 2**n: #constant F
        return (False,can_inputs,can_outputs,F,can_order)
    if type(variables)==np.ndarray:
        variables = list(variables)
    if variables == []:
        variables = list(range(n))
    if type(F) == list:
        F = np.array(F)
    desired_value = 2**(n-1)
    T = np.array(list(itertools.product([0, 1], repeat=n))).T
    A = np.r_[T,1-T]
    try: #is 1 canalized output for one of the variables
        index = list(np.dot(A,F)).index(desired_value)
        newF = F[np.where(A[index]==0)[0]]
        variable = variables.pop(index%n)
        return is_k_canalizing_return_inputs_outputs_corefunction_order(newF,k-1,n-1,np.append(can_inputs,int(index<n)),np.append(can_outputs,1),np.append(can_order,variable),variables)
    except ValueError:
        try: #is 0 canalized output for one of the variables
            index = list(np.dot(A,1-F)).index(desired_value)
            newF = F[np.where(A[index]==0)[0]]
            variable = variables.pop(index%n)
            return is_k_canalizing_return_inputs_outputs_corefunction_order(newF,k-1,n-1,np.append(can_inputs,int(index<n)),np.append(can_outputs,0),np.append(can_order,variable),variables)
        except ValueError:
            return (False,can_inputs,can_outputs,F,can_order)    

def get_all_canalizing_variables_of_boolean_function(rule,or_symbol = 'or', and_symbol = 'and', not_symbol='not',can_inputs=np.array([],dtype=int),can_outputs=np.array([],dtype=int),can_variables=np.array([],dtype=int)):
    #requires the rule as a generalized DNF or CNF, several smaller parts separated by or (DNF) or and (CNF)
    #rule='(   (  x24  &   (   (   (  x4   |  x13   |  x10   |  x14   |  x9   |  x22   |  x25   |  x7   |  x3   |  x21   |  x17   |  x6   |  x18   |  x15   |  x26  )   )   )       )   &   ~   (  x19   )   )    |   (   (  x8   )   &   ~   (  x19   )   )    |   (   (  x5  &   (   (   (  x4   |  x13   |  x10   |  x14   |  x9   |  x22   |  x25   |  x7   |  x3   |  x21   |  x17   |  x6   |  x18   |  x15   |  x26  )   )   )       )   &   ~   (  x19   )   )    |   (   (  x12  &   (   (   (  x4   |  x13   |  x10   |  x14   |  x9   |  x22   |  x25   |  x7   |  x3   |  x21   |  x17   |  x6   |  x18   |  x15   |  x26  )   )   )       )   &   ~   (  x19   )   )    |   (   (  x23   )   &   ~   (  x19   )   )    |   (   (  x1  &   (   (   (  x4   |  x13   |  x10   |  x14   |  x9   |  x22   |  x25   |  x7   |  x3   |  x21   |  x17   |  x6   |  x18   |  x15   |  x26  )   )   )       )   &   ~   (  x19   )   )    |   (   (  x16  &   (   (   (  x4   |  x13   |  x10   |  x14   |  x9   |  x22   |  x25   |  x7   |  x3   |  x21   |  x17   |  x6   |  x18   |  x15   |  x26  )   )   )       )   &   ~   (  x19   )   )    |   (   (  x0  &   (   (   (  x4   |  x13   |  x10   |  x14   |  x9   |  x22   |  x25   |  x7   |  x3   |  x21   |  x17   |  x6   |  x18   |  x15   |  x26  )   )   )       )   &   ~   (  x19   )   )    |   (   (  x20  &   (   (   (  x4   |  x13   |  x10   |  x14   |  x9   |  x22   |  x25   |  x7   |  x3   |  x21   |  x17   |  x6   |  x18   |  x15   |  x26  )   )   )       )   &   ~   (  x19   )   )    |   (   (  x11  &   (   (   (  x4   |  x13   |  x10   |  x14   |  x9   |  x22   |  x25   |  x7   |  x3   |  x21   |  x17   |  x6   |  x18   |  x15   |  x26  )   )   )       )   &   ~   (  x19   )   )    |   (   (  x2  &   (   (   (  x4   |  x13   |  x10   |  x14   |  x9   |  x22   |  x25   |  x7   |  x3   |  x21   |  x17   |  x6   |  x18   |  x15   |  x26  )   )   )       )   &   ~   (  x19   )   )'
    python_or = '|'
    python_and = '&'
    python_not = '~'
    mod_rule = rule.replace(or_symbol,' %s ' % or_symbol).replace(and_symbol,' %s ' % and_symbol).replace(not_symbol,' %s ' % not_symbol).replace(or_symbol,python_or).replace(and_symbol,python_and).replace(not_symbol,python_not).replace('[','').replace(']','')
    INRULE = False
    parts = []
    count_open_parantheses=0
    symbols_between_rules = []
    for char in mod_rule:
        if char=='(':
            count_open_parantheses+=1
            if count_open_parantheses==1 and INRULE==False:
                parts.append('')
                INRULE = True
            else:
                parts[-1]+='('
        elif char==')':
            count_open_parantheses-=1
            if count_open_parantheses==0:
                if parts[-1][0]==python_not:
                    parts[-1]+=')'
                INRULE=False
            else:
                parts[-1]+=')'
        elif INRULE:
            parts[-1] += char
        elif char!=' ':
            if char==python_not:
                parts.append(python_not+' ')
                INRULE = True
            else:
                symbols_between_rules.append(char)
    IS_DNF = python_or in set(symbols_between_rules) and python_and not in set(symbols_between_rules)
    IS_CNF = python_and in set(symbols_between_rules) and python_or not in set(symbols_between_rules)
    assert IS_DNF==True or IS_CNF==True
    var = list(set(mod_rule.split(' '))-set([python_or,python_and,python_not,')','(','']))
    dict_var = dict(zip(var,range(len(var))))
    
    #check for each part of the rule if a variable canalizes it
    potential_canalizing_inputs = []
    potential_canalizing_vars = []
    potential_canalizing_outputs = []
    count_vars = [0]*len(var)
    keep_part=np.ones(len(parts),dtype=bool)
    for i,part in enumerate(parts):
        var_in_part = list(set(part.split(' '))-set([python_or,python_and,python_not,')','(','','false','true']))
        if len(var_in_part)==0:
            keep_part[i] = False
            
        for v in var_in_part:
            count_vars[dict_var[v]] += 1
            for el in [' false ',' true ']:
                simplified_expression = sympy.simplify(part.replace(' '+v+' ',el))
                if simplified_expression in [False,True]:
                    print(v,el,simplified_expression)
                    potential_canalizing_inputs.append(1 if 'true' else 0)
                    potential_canalizing_vars.append(v)
                    potential_canalizing_outputs.append(1 if simplified_expression==True else 0)
    parts = np.array(parts)[keep_part]
    
    count_potential_canalizing_inputs = pd.value_counts(potential_canalizing_vars)
    canalizing_variables = []
    canalizing_inputs = []
    canalizing_outputs = []
    subrule = mod_rule
    if IS_DNF:
        #if one variable makes at least one part True, this variable is canalizing
        for i in range(len(potential_canalizing_outputs)):
            if potential_canalizing_outputs[i]==True:
                if potential_canalizing_vars[i] not in canalizing_variables:
                    canalizing_variables.append(potential_canalizing_vars[i])
                    canalizing_inputs.append(potential_canalizing_inputs[i])
                    canalizing_outputs.append(potential_canalizing_outputs[i])
                    subrule = subrule.replace(' '+canalizing_variables[-1]+' ',' false ' if canalizing_inputs[-1]==1 else ' true ')
        #otherwise, if one variable occurs in every part and makes it False using the same input every time, this variable is canalizing
        for i in range(len(count_potential_canalizing_inputs)):
            if count_potential_canalizing_inputs[i]<len(parts):
                break
            elif count_potential_canalizing_inputs.index[i] in canalizing_variables: #already found this variable to be canalizing
                continue
            inputs = set(np.array(potential_canalizing_inputs)[np.array(potential_canalizing_vars)==count_potential_canalizing_inputs.index[i]])
            if len(inputs)==1:
                canalizing_variables.append(count_potential_canalizing_inputs.index[i])
                canalizing_inputs.append(int(list(inputs)[0]))
                canalizing_outputs.append(0)
                subrule = subrule.replace(' '+canalizing_variables[-1]+' ',' false ' if canalizing_inputs[-1]==1 else ' true ')
    #print(canalizing_variables)
    rule=subrule
    if canalizing_variables == []:
        return rule, can_variables, can_inputs, can_outputs
    else:
        return get_all_canalizing_variables_of_boolean_function(rule,can_inputs=np.append(can_inputs,canalizing_inputs),can_outputs=np.append(can_outputs,canalizing_outputs),can_variables=np.append(can_variables,canalizing_variables))

def average_sensitivity_old_wrong(F,nsim=10000): 
    #equals Derrida value D(F,1) if all n update rules in F are chosen from the same sampling space
    if type(F)==list:
        F = np.array(F)
    n = int(np.log2(len(F)))
    num_values = 2**n
    X = np.random.randint(num_values,size=nsim)
    add = 2**np.random.randint(n, size=nsim)
    Y = (X + add )%num_values
    return sum(np.bitwise_xor(F[X],F[Y]))*1./nsim

def average_sensitivity(F,nsim=10000,EXACT=False,NORMALIZED=True): 
    #equals Derrida value D(F,1) if all n update rules in F are chosen from the same sampling space
    if type(F)==list:
        F = np.array(F)
    n = int(np.log2(len(F)))
    num_values = 2**n
    s = 0
    if EXACT:
        bool_list = list(map(np.array,list(itertools.product([0, 1], repeat = n))))
        for ii,X in enumerate(bool_list):
            for i in range(n):
                Y=X.copy()
                Y[i] = 1-X[i]
                Ydec = bin2dec(Y)
                s += int(F[ii]!=F[Ydec])
        if NORMALIZED:
            return s/num_values/n
        else:
            return s/num_values
    else:
        for i in range(nsim):
            Xbin = np.random.randint(num_values)
            Y = dec2bin(Xbin,n)
            index = np.random.randint(n)
            Y[index] = 1-Y[index]
            Ybin = bin2dec(Y)
            s += int(F[Xbin]!=F[Ybin])
        if NORMALIZED:
            return s/nsim
        else:
            return n*s/nsim


def absolute_bias(F,n=None):
    if n==None:
        n = int(np.log2(len(F)))
    return abs(sum(F)*1./2**(n-1)-1)
    
## 2) Put everything together to obtain canalizing depth, layer structure, canalized outputs, canalizing inputs as well as core function (could also calculate canalizing variables in future versions but I don't see a need)

def get_canalizing_depth_inputs_outputs_corefunction(F):
    n = int(np.log2(len(F)))
    (NESTED_CANALIZING,can_inputs,can_outputs,corefunction) = is_k_canalizing_return_inputs_outputs_corefunction(F,n,n)
    return (n,len(can_inputs),can_inputs,can_outputs,corefunction)
   
def get_canalizing_depth_inputs_outputs_corefunction_order(F,variables = []):
    n = int(np.log2(len(F)))
    (NESTED_CANALIZING,can_inputs,can_outputs,corefunction,can_order) = is_k_canalizing_return_inputs_outputs_corefunction_order(F,n,n,variables=variables)
    return (n,len(can_inputs),can_inputs,can_outputs,corefunction,can_order)    

def get_layer_structure_given_outputs_corefunction(can_outputs_F,corefunction_F,n_F):
    k = len(can_outputs_F)
    if k == 0:
        return []
    if k == n_F and n_F>1: #The last layer of Boolean NCFs has size >=2
        can_outputs_F[-1] = can_outputs_F[-2]
    elif is_constant(corefunction_F) and k>1: #Exceptional case, again last layer here needs to be size >=2
        can_outputs_F[-1] = can_outputs_F[-2]
    kis = []
    ki = 1
    for i in range(1,k):
        if can_outputs_F[i]==can_outputs_F[i-1]:
            ki+=1
        else:
            kis.append(ki)
            ki=1
    kis.append(ki)
    return kis
    
def kindoflayer(k,w):
    '''For NCFs (n-canalizing functions) there is a bijection between the Hamming weight (assuming w is equivalent to 2^n-w) and the layer structure 
    
    Input: k - number of nodes that influence the output,
       w - odd Hamming weight of specific NCF
    Output: r - number of layers of this rule
       ki - vector of kis, [k1, ..., kr]'''
         
    if w==1:
        r=1
        ki=[k]
    else:
        assert type(w) == int or type(w) == np.int64
        assert 1<=w<=2**k-1
        assert w%2==1
        w_bin=dec2bin(w,k)
        current_el=w_bin[0]
        ki=[1]
        for el in w_bin[1:-1]:
            if el==current_el:
                ki[-1]+=1
            else:
                ki.append(1)
                current_el=el
        ki[-1]+=1
        r=len(ki)
    return (r,ki)


## 3) Methods to randomly generate Boolean functions (uniform distribution) and Boolean networks
def random_function(n):
    return np.random.randint(2, size = 2**n)    

def random_non_degenerated_function(n):
    while True: #works because most functions are non-degenerated
        F = np.random.randint(2, size = 2**n) 
        if not is_degenerated(F):
            return F

def random_non_canalizing_function(n):
    assert n>1
    while True: #works because most functions are non-canalizing
        F = np.random.randint(2, size = 2**n) 
        if not is_canalizing(F,n):
            return F

def random_non_canalizing_non_degenerated_function(n):
    assert n>1
    while True: #works because most functions are non-canalizing and non-degenerated
        F = np.random.randint(2, size = 2**n) 
        if not is_canalizing(F,n) and not is_degenerated(F):
            return F

def random_k_canalizing(n, k, EXACT_DEPTH_K=False, x=[]):
    try:
        assert (n-k!=1 or EXACT_DEPTH_K==False)
    except AssertionError:
        print('There are no functions of exact canalizing depth n-1.\nEither set EXACT_DEPTH_K=False or ensure k!=n-1')
        return
    try:
        assert 0<=k and k<=n
    except AssertionError:
        print('Error:\nEnsure 0 <= k <= n.')
        return
    if x==[]: #to reduce run time, this should be calculated once and then passed as argument
        x = list(itertools.product([0, 1], repeat = n))
    num_values = 2**n
    aas = np.random.randint(2, size = k)  # inputs
    bbs = np.random.randint(2, size = k)  # outputs
    can_vars = np.random.choice(n, k, replace = False)
    F = np.zeros(num_values, dtype = int)
    
    if k<n:
        if EXACT_DEPTH_K:
            core_polynomial = random_non_canalizing_non_degenerated_function(n-k)
        else:
            core_polynomial = random_non_degenerated_function(n-k)    
    else:
        core_polynomial = [1-bbs[-1]]
        
    counter_non_canalized_positions = 0
    
    for i in range(num_values):
        for j in range(k):
            if x[i][can_vars[j]] == aas[j]:
                F[i] = bbs[j]
                break
        else:
            F[i] = core_polynomial[counter_non_canalized_positions]
            counter_non_canalized_positions += 1
    return F

def random_k_canalizing_with_specific_weight(n, kis, EXACT_DEPTH_K=False, x=[]):
    k=sum(kis)
    if k==0:
        kis = [0]
    try:
        assert (n-k!=1 or EXACT_DEPTH_K==False)
    except AssertionError:
        print('Error:\nThere are no functions of exact canalizing depth n-1.\nEither set EXACT_DEPTH_K=False or ensure k=sum(kis)!=n.')
        return
    try:
        assert 0<=k and k<=n
    except AssertionError:
        print('Error:\nEnsure 0 <= k=sum(kis) <= n.')
        return
    try:
        assert k<n or kis[-1]>1 or n==1
    except AssertionError:
        print('Error:\nThe last layer of an n-canalizing function (NCF) has to have size >=2 for n>1.\nIf k=sum(kis)=n, ensure that kis[-1]>=2.')
        return
    try:
        assert min(kis)>=1
    except AssertionError:
        print('Error:\nThere needs to be at least one variable in each layer, i.e., each element of kis must be >=1.')
        return
    if x==[]: #to increase run time, this should be calculated once and then passed as argument
        x = list(itertools.product([0, 1], repeat = n))
    num_values = 2**n
    aas = np.random.randint(2, size = k)  # inputs
    b0 = np.random.randint(2)
    bbs = [b0]*kis[0]  # outputs
    for i in range(1,len(kis)):
        if i%2==0:
            bbs.extend([b0]*kis[i])
        else:
            bbs.extend([1-b0]*kis[i])
    can_vars = np.random.choice(n, k, replace = False)
    F = np.zeros(num_values, dtype = int)
    
    if k<n:
        if EXACT_DEPTH_K:
            core_polynomial = random_non_canalizing_non_degenerated_function(n-k)
        else:
            core_polynomial = random_non_degenerated_function(n-k)    
    else:
        core_polynomial = [1-bbs[-1]]
        
    counter_non_canalized_positions = 0
    
    for i in range(num_values):
        for j in range(k):
            if x[i][can_vars[j]] == aas[j]:
                F[i] = bbs[j]
                break
        else:
            F[i] = core_polynomial[counter_non_canalized_positions]
            counter_non_canalized_positions += 1
    return F

def random_adj_matrix(N,ns,NO_SELF_REGULATION=True,STRONGLY_CONNECTED=False): #recursive function definition
    matrix = np.zeros((N, N), dtype = int)
    indices = []
    for i in range(N):
        if NO_SELF_REGULATION:
            indexes = np.random.choice(np.append(np.arange(i),np.arange(i+1,N)), ns[i], replace=False)
        else:
            indexes = np.random.choice(np.arange(N), ns[i],replace=False)         
        indexes = np.sort(indexes)
        indices.append(indexes)
        for index in indexes:
            matrix[i][index] = 1
    if STRONGLY_CONNECTED:
        G = nx.from_numpy_array(matrix, parallel_edges = False, create_using = nx.MultiDiGraph())
        if not nx.is_strongly_connected(G):
            return random_adj_matrix(N,ns,NO_SELF_REGULATION,STRONGLY_CONNECTED)
    return (matrix, indices)

def random_edge_list(N,ns,NO_SELF_REGULATION,AT_LEAST_ONE_REGULATOR_PER_GENE=False):
    if AT_LEAST_ONE_REGULATOR_PER_GENE==False:
        edge_list = []
        for i in range(N):
            if NO_SELF_REGULATION:
                indices = np.random.choice(np.append(np.arange(i),np.arange(i+1,N)), ns[i], replace=False)
            else:
                indices = np.random.choice(np.arange(N), ns[i],replace=False)
            edge_list.extend( list(zip(indices,i*np.ones(ns[i],dtype=int))) )
    else:
        edge_list = []
        outdegree = np.zeros(N,dtype=int)
        sum_ns = sum(ns) #total number of regulations
        for i in range(N):
            if NO_SELF_REGULATION:
                indices = np.random.choice( np.append(np.arange(i),np.arange(i+1,N)), ns[i], replace=False)
            else:
                indices = np.random.choice(np.arange(N), ns[i],replace=False)
            outdegree[indices] += 1
            edge_list.extend( list(zip(indices,i*np.ones(ns[i],dtype=int))) )
        while min(outdegree)==0:
            index_sink = np.where(outdegree==0)[0][0]
            index_edge = int(random.random()*sum_ns)
            if NO_SELF_REGULATION:
                while edge_list[index_edge][1] == index_sink:
                    index_edge = int(random.random()*sum_ns)
            outdegree[index_sink] += 1
            outdegree[edge_list[index_edge][0]] -= 1
            edge_list[index_edge] = (index_sink,edge_list[index_edge][1])
    return edge_list

def random_BN(N, n = 2, k = 0, STRONGLY_CONNECTED = True, indegree_distribution = 'constant', list_x=[], kis = None, EXACT_DEPTH=False,NO_SELF_REGULATION=True):    
    #need to also accept vectors for k and kis
    if indegree_distribution in ['constant','dirac','delta']:
        if type(n) in [list,np.array]:
            try:
                assert type(n) in [list,np.array]
                assert np.all([type(el) in [int,np.int64] for el in n])
                assert len(n)==N
                assert min(n)>=1
                assert max(n)<=N
                ns = np.array(n[:])
            except AssertionError:
                print('Error: A vector n was submitted.\nTo use a user-defined in-degree vector, ensure that n is an N-dimensional vector where each element of n is an integer between 1 and N.')
                return
        else:
            try:
                assert type(n) in [int,np.int64]
                assert n>=1
                assert n<=N
                ns = np.ones(N,dtype=int)*n
            except AssertionError:
                print('Error: n must be a single integer (or N-dimensional vector of integers) between 1 and N when using a constant degree distribution.')
                return            
    elif indegree_distribution == 'uniform':
        if type(n) in [list,np.array]:
            try:
                assert type(n) in [list,np.array]
                assert np.all([type(el) in [int,np.int64] for el in n])
                assert len(n)==N
                assert min(n)>=1
                assert max(n)<=N
                ns = np.array(n[:])
            except AssertionError:
                print('Error: A vector n was submitted.\nEnsure that n is an N-dimensional vector where each element of n is an integer between 1 and N representing the upper bound of a uniform degree distribution (lower bound==1).')
                return
        else:
            try:
                assert type(n) in [int,np.int64]
                assert n>=1
                assert n<=N
                ns = np.ones(N,dtype=int)*n
            except AssertionError:
                print('Error: n must be a single integer (or N-dimensional vector of integers) between 1 and N representing the upper bound of a uniform degree distribution (lower bound==1).')
                return
    elif indegree_distribution == 'poisson':
        if type(n) in [list,np.array]:
            try:
                assert type(n) in [list,np.array]
                assert np.all([type(el) in [int,np.int64,float,np.float64] for el in n])
                assert len(n)==N
                assert min(n)>=1
                assert max(n)<=N
                ns = np.array(n[:])
            except AssertionError:
                print('Error: A vector n was submitted.\nEnsure that n is an N-dimensional vector where each element of n is > 0 and represents the Poisson parameter.')
                return
        else:
            try:
                assert type(n) in [int,np.int64,float,np.float64]
                assert n>=1
                assert n<=N
                ns = np.ones(N,dtype=int)*n
            except AssertionError:
                print('Error: n must be a single number (or N-dimensional vector) > 0 representing the Poisson parameter.')
                return            
    else:
        print('None of the predefined indegree distributions were chosen.\nTo use a user-defined in-degree vector, use the input n to submit an N-dimensional vector where each element of n must be between 1 and N')
        return

    if kis==None:
        if type(k) in [int,np.int64]:
            try:
                assert k>=0
                assert k<=N
                max_k = k
            except AssertionError:
                print('Error: k must be an integer between 0 and N.')
                return
        elif type(k) in [list,np.array]:
            try:
                assert len(k)==N
                assert np.all([type(el) in [int,np.int64] for el in k])
                max_k = max(k)
                assert min(k)>=0
                assert max_k<=N
            except AssertionError:
                print('Error: A vector k was submitted.\nTo use a user-defined vector k, ensure that k is an N-dimensional vector where each element of k is an integer between 0 and N.')
                return
        else:
            print('Error: Wrong input format for k.\nk must be a single integer (or N-dimensional vector of integers) between 0 and N')
            return
    else: #kis provided
        if np.all([type(el) in [int,np.int64] for el in kis]):
            try:
                assert np.all([type(el) in [int,np.int64] for el in kis])
                assert np.all([el>=1 for el in kis])
                max_k = sum(kis)
                assert max_k<=N #could be more restrictive but if n is also a vector this is tedious to program
            except AssertionError:
                print('Error: the layer structure kis must be a vector of positive integers with 0<= k=sum(kis) <= n.')
                return
        elif np.all([type(el) in [list,np.array] for el in kis]):
            try:
                assert len(kis)==N
                assert type(kis[0][0]) in [int,np.int64]
                max_k = max([sum(el) for el in kis])
                assert min([min(el) for el in kis])>=0
                assert max_k<=N
            except AssertionError:
                print('Error: A vector of kis was submitted.\nTo use a user-defined vector of kis, ensure that kis is an N-dimensional vector where each element represents a layer structure and is a vector of positive integers with 1<= k=sum(kis[i]) <= n.')
                return
        else:
            print('Error: Wrong input format for kis.\nkis must be a single vector (or N-dimensional vector of layer structures) where each the sum of each element must be between 0 and N.')
            return

    while True: # Keep generating until we have a strongly connected graph
        # create a new degree vector
        # Can be constant, random or a statistical distribution
        if indegree_distribution == 'uniform':            
            ns = 1 + np.random.randint(n-1, size = N)
        elif indegree_distribution == 'poisson':
            ns = np.random.poisson(lam = n , size = N)
            ns[ns==0] = 1
            ns[ns>N-int(NO_SELF_REGULATION)] = N-int(NO_SELF_REGULATION)
            
        #A, I = random_adj_matrix(N, ns, NO_SELF_REGULATION)
        E = random_edge_list(N,ns,NO_SELF_REGULATION)
        
        # If we care whether the graph is strongly connected (STRONGLY_CONNECTED),
        # we will check the connectedness of the graph using nx.is_strongly_connected(G)
        if STRONGLY_CONNECTED:
            #G = nx.from_numpy_array(A, parallel_edges = False, create_using = nx.MultiDiGraph())
            G = nx.from_edgelist(E, create_using = nx.MultiDiGraph())
            strongly_connected = nx.is_strongly_connected(G)
            if not strongly_connected:
                continue
        break

    max_n = max(ns)
    if max_k>0 and list_x==[] or len(list_x)<max_n: #list_x probably not correct, or fully generated
        #create list_x
        list_x = [[[0],[1]]]
        list_x.extend([list(itertools.product([0, 1], repeat = nn)) for nn in range(2, max_n+1)])

        
    F = []
    for i in range(N):
        if k>0 and kis==None:
            if type(k) in [int,np.int64]:
                F.append(random_k_canalizing(ns[i], k, EXACT_DEPTH_K = EXACT_DEPTH, x=list_x[ns[i]-1]))
            else:
                F.append(random_k_canalizing(ns[i], k[i], EXACT_DEPTH_K = EXACT_DEPTH, x=list_x[ns[i]-1]))
        elif kis!=None: #value of k is ignored if a layer structure is provided
            if np.all([type(el) in [int,np.int64] for el in kis]):
                F.append(random_k_canalizing_with_specific_weight(ns[i], kis, EXACT_DEPTH_K = EXACT_DEPTH, x=list_x[ns[i]-1]))                
            else:
                F.append(random_k_canalizing_with_specific_weight(ns[i], kis[i], EXACT_DEPTH_K = EXACT_DEPTH, x=list_x[ns[i]-1]))                                
        else:
            if EXACT_DEPTH==True: #i.e. if k==0
                F.append(random_non_canalizing_non_degenerated_function(ns[i]))                   
            else:
                F.append(random_non_degenerated_function(ns[i]))   
    
    I = [[] for _ in range(N)]
    for edge in E:
        I[edge[1]].append(edge[0])
    for i in range(N):
        I[i] = np.sort(I[i])
    
    return F, I, ns



## 4) Enumeration methods
def nr_non_canalizing_by_weight_exact(n):
    assert n<=4
    
    nr_values = 2**n
    nr_boolean_fcts = 2**(nr_values)

    F = np.array(list(itertools.product([0, 1], repeat=nr_values))).T
    
    ws_can,ws_non_can = [],[]
    for i in range(nr_boolean_fcts):
        if is_canalizing(F[:,i],n):
            ws_can.append(sum(F[:,i]))
        else:
            ws_non_can.append(sum(F[:,i]))
    return [ws_non_can,ws_can]

def nr_non_canalizing_by_weight_simulation(n,nsim=10000):
    nr_values = 2**n
    
    ws_non_can,ws_can = [],[]
    
    for i in range(nsim):
        F = np.random.randint(2, size = nr_values)
        if is_canalizing(F,n):
            ws_can.append(sum(F))
        else:
            ws_non_can.append(sum(F))
    return [ws_non_can,ws_can]

def stratify_Boolean_fcts_by_canalization_ns(n,nsim=10000):
    nr_values = 2**n  
    depths,kiss = [],[]
    for i in range(nsim):
        F = np.random.random(nr_values)>0.5
        (n_F,depth_F,can_inputs_F,can_outputs_F,corefunction_F) = get_canalizing_depth_inputs_outputs_corefunction(F)
        depths.append(depth_F)
        kis = get_layer_structure_given_outputs_corefunction(can_outputs_F,corefunction_F,n)
        kiss.append(kis)
    return (depths,kiss)


## 5) Analysis methods
def get_constant_nodes(I,degree,N):
    return np.array([i for i in range(N) if degree[i]==1 and I[i][0]==i])

def rnd_edge_deletion(F,I,N,degree,nsim=100,bool_lists=[]):
    res = []
    for i in range(nsim):
        rnd_var = np.random.randint(N)
        rnd_input = np.random.randint(degree[rnd_var])
        res.append( sum(F[rnd_var][np.array(bool_lists[degree[rnd_var]-1])[:,rnd_input]==0]!=F[rnd_var][np.array(bool_lists[degree[rnd_var]-1])[:,rnd_input]==1]) )
    return res

def update(F, I, N, X):
    Fx = np.zeros(N, dtype = int)
    for i in range(N):
        Fx[i] = F[i][bin2dec(X[I[i]])]
    return Fx    

def num_of_attractors(F, I, N, nsim = 500, EXACT = False, bool_list = []):
    dictF = dict()
    attractors = []
    basin_sizes = []
    attr_dict = dict()
    
    if EXACT and bool_list == []:
        bool_list = list(map(np.array,list(itertools.product([0, 1], repeat = N))))
    
    for i in range(nsim if not EXACT else 2**N):
        x = np.random.randint(2, size = N) if not EXACT else bool_list[i]
        xbin = bin2dec(x)
        queue = [xbin]        
        while True: #as long as we haven't reached an attractor state, keep updating
            try:
                fxbin = dictF[xbin]
            except KeyError:
                fx = update(F, I, N, x)
                fxbin = bin2dec(fx)
                dictF.update({xbin:fxbin})
                x = fx #works only because if we don't know fx now we also won't know F[fx] 
            try: # check if this state is a known attractor
                index_attr = attr_dict[fxbin] #returns a KeyError if we haven't identified fxbin as an attractor before
                basin_sizes[index_attr] += 1
                break
            except KeyError:
                try: #check if fxbin is part of a new attractor
                    index = queue.index(fxbin) #returns a ValueError if fxbin is not part of queue
                    #new attractor
                    attr_dict.update( list(zip( queue[index:] , [len(attractors)]*(len(queue)-index) )) )
                    attractors.append(queue[index:])
                    basin_sizes.append(1)
                    break
                except ValueError:
                    pass
            queue.append(fxbin)
            xbin = fxbin
    return (attractors, len(attractors), basin_sizes)

def num_of_attractors_v2(F, I, N, nb = 500): #should be faster if we have longer average path lengths, not significantly though
    dictF = dict()
    attractors = []
    basin_sizes = []
    attr_dict = dict()
    
    for i in range(nb):
        x = np.random.randint(2, size = N)
        xbin = bin2dec(x)
        queue = [xbin]        
        while True: #as long as we haven't reached an attractor state, keep updating
            try:
                fxbin = dictF[xbin]
            except KeyError:
                fx = update(F, I, N, x)
                fxbin = bin2dec(fx)
                dictF.update({xbin:fxbin})
                x = fx #works only because if we don't know fx now we also won't know F[fx] 
            try: # check if this state has a known attractor
                index_attr = attr_dict[fxbin] #returns a KeyError if we haven't identified fxbin as an attractor before
                basin_sizes[index_attr] += 1
                attr_dict.update( list(zip( queue , [index_attr]*len(queue) )) )
                break
            except KeyError:
                try: #check if fxbin is part of a new attractor
                    index = queue.index(fxbin) #returns a ValueError if fxbin is not part of queue
                    #new attractor
                    attr_dict.update( list(zip( queue[index:] , [len(attractors)]*(len(queue)-index) )) )
                    attractors.append(queue[index:])
                    basin_sizes.append(1)
                    break
                except ValueError:
                    pass
            queue.append(fxbin)
            xbin = fxbin
    return (attractors, len(attractors), basin_sizes)
        
def basin_size_largest(basin_sizes):
    return max(basin_sizes)*1./sum(basin_sizes)

def entropy(basin_sizes):
    total = sum(basin_sizes)
    return sum([-np.log(el)*el for el in [size*1./total for size in basin_sizes]])

def d(x, y):            
    return sum(x!=y) #works as long as x,y are numpy arrays

def derrida_value(F, I, N, m, nsim = 500):
    total = 0
    for i in range(nsim):
        X = np.random.randint(2, size = N)
        ones = np.zeros(N, dtype = int)
        ones[np.random.choice(N, m, replace = False)] = 1
        Y = np.bitwise_xor(X, ones)
        total += d(update(F, I, N, X), update(F, I, N, Y))
    return total*1./nsim

def adjacency_matrix(I,constants=[],IGNORE_SELFLOOPS=False,IGNORE_CONSTANTS=True):
    n = len(I)
    n_constants = len(constants)
    if IGNORE_CONSTANTS:
        m = np.zeros((n-n_constants,n-n_constants),dtype=int)
        for i,regulators in enumerate(I):
            for j in regulators:
                if j<n-n_constants and (IGNORE_SELFLOOPS==False or i!=j):
                    m[j,i] = 1
        return m
    else:
        return adjacency_matrix(I,[],IGNORE_CONSTANTS=True)

def get_ffls(adjacency_matrix,F=None,I=None):
    n = len(adjacency_matrix)
    ffls = []
    for i in range(n):
        for j in range(n):
            if i==j or adjacency_matrix[i,j]==0:
                continue
            #else: find other nodes regulated by i 
            for k in range(n):
                if j==k or i==k or adjacency_matrix[i,k]==0:
                    continue
                #else: check if that other node k also regulates i
                if adjacency_matrix[k,j]!=0: #found FFL
                    ffls.append([i,k,j])
    if F==None or I==None:
        return ffls
    dict_monotonic = dict()
    types = []
    for [j,k,i] in ffls:
        try:
            monotonic = dict_monotonic[k]
        except KeyError:
            monotonic = is_monotonic(F[k],True)[1]
            dict_monotonic.update({k : monotonic})
        indirect1 = monotonic[list(I[k]).index(j)]
        try:
            monotonic = dict_monotonic[i]
        except KeyError:
            monotonic = is_monotonic(F[i],True)[1]
            dict_monotonic.update({i : monotonic})
        indirect2 = monotonic[list(I[i]).index(k)]
        direct = monotonic[list(I[i]).index(j)]
        types.append([direct,indirect1,indirect2])
    return (ffls,types)

def get_ffls_from_I(I,types_I=None):
    all_tfs = list(range(len(I)))#list(set.union(*[] + [set(el) for el in I]))
    n_tfs = len(all_tfs)
    all_tfs_dict = dict(zip(all_tfs,list(range(n_tfs))))
    I_inv = [[] for _ in all_tfs]
    for target,el in enumerate(I):
        for regulator in el:
            I_inv[all_tfs_dict[regulator]].append(target)
    ffls = []
    types = []
    for i in range(n_tfs): #master regulators
        for j in range(n_tfs):
            if i==j or all_tfs[j] not in I_inv[i]:
                continue
            #else: find common nodes regulated by both
            common_targets = list(set(I_inv[i]) & set(I_inv[j]))
            for k in common_targets:
                if all_tfs[j]==k or all_tfs[i]==k:
                    continue
                ffls.append([i,j,k])
                if types_I != None:
                    direct = types_I[k][I[k].index(all_tfs[i])]
                    indirect1 = types_I[all_tfs[j]][I[all_tfs[j]].index(all_tfs[i])]
                    indirect2 = types_I[k][I[k].index(all_tfs[j])]
                    types.append([direct,indirect1,indirect2])
    if types_I != None:
        return (ffls,types)
    else:
        return ffls
    
def get_ffl_type_number(types_vector):
    if not set(types_vector).issubset(set(['decreasing', 'increasing'])):
        return -1 if 'not essential' not in types_vector else -2
    else:
        dummy = np.array([1 if el=='increasing' else 0 for el in types_vector])
        nr_type = np.dot(dummy,2**np.arange(len(types_vector)))
    return nr_type

def is_ffl_coherent(types_vector):
    if not set(types_vector).issubset(set(['decreasing', 'increasing'])):
        return np.nan
    else:
        dummy = np.array([1 if el=='increasing' else 0 for el in types_vector])
        COHERENT = sum(dummy)%2==1
    return COHERENT

def simple_cycles(G,max_len=4):
    def _unblock(thisnode, blocked, B):
        stack = set([thisnode])
        while stack:
            node = stack.pop()
            if node in blocked:
                blocked.remove(node)
                stack.update(B[node])
                B[node].clear()
    
    # Johnson's algorithm requires some ordering of the nodes.
    # We assign the arbitrary ordering given by the strongly connected comps
    # There is no need to track the ordering as each node removed as processed.
    # Also we save the actual graph so we can mutate it. We only take the
    # edges because we do not want to copy edge and node attributes here.
    subG = type(G)(G.edges())
    sccs = [scc for scc in nx.strongly_connected_components(subG)
            if len(scc) > 1]
    
    # Johnson's algorithm exclude self cycle edges like (v, v)
    # To be backward compatible, we record those cycles in advance
    # and then remove from subG
    for v in subG:
        if subG.has_edge(v, v):
            yield [v]
            subG.remove_edge(v, v)
    
    while sccs:
        scc = sccs.pop()
        sccG = subG.subgraph(scc)
        # order of scc determines ordering of nodes
        startnode = scc.pop()
        # Processing node runs "circuit" routine from recursive version
        path = [startnode]
        len_path = 1
        blocked = set()  # vertex: blocked from search?
        closed = set()   # nodes involved in a cycle
        blocked.add(startnode)
        B = defaultdict(set)  # graph portions that yield no elementary circuit
        stack = [(startnode, list(sccG[startnode]))]  # sccG gives comp nbrs
        while stack:
            thisnode, nbrs = stack[-1]
            if nbrs and len_path<=max_len:
                nextnode = nbrs.pop()
                if nextnode == startnode:
                    yield path[:]
                    closed.update(path)
                elif nextnode not in blocked:
                    path.append(nextnode)
                    len_path+=1
                    stack.append((nextnode, list(sccG[nextnode])))
                    closed.discard(nextnode)
                    blocked.add(nextnode)
                    continue
            # done with nextnode... look for more neighbors
            if not nbrs or len_path>max_len:  # no more nbrs
                if thisnode in closed:
                    _unblock(thisnode, blocked, B)
                else:
                    for nbr in sccG[thisnode]:
                        if thisnode not in B[nbr]:
                            B[nbr].add(thisnode)
                stack.pop()
                path.pop()
                len_path-=1
        # done processing this node
        H = subG.subgraph(scc)  # make smaller to avoid work in SCC routine
        sccs.extend(scc for scc in nx.strongly_connected_components(H)
                    if len(scc) > 1)

def get_type_of_loop(loop,F,I):
    n = len(loop)
    dummy = loop[:]
    dummy.append(loop[0])
    res = []
    for i in range(n):
        res.append( is_monotonic(F[dummy[i+1]],True)[1][list(I[dummy[i+1]]).index(dummy[i])] )
    return res

def get_loop_type_number(types_vector):
    if not set(types_vector).issubset(set(['decreasing', 'increasing'])):
        return -1 if 'not essential' not in types_vector else -2
    else:
        nr_type = int(np.sum([1 if el=='decreasing' else 0 for el in types_vector]))
    return nr_type

def is_pos_loop(types_vector):
    if not set(types_vector).issubset(set(['decreasing', 'increasing'])):
        return np.nan
    else:
        POSITIVE = int(np.sum([1 if el=='decreasing' else 0 for el in types_vector])) % 2 == 0
    return POSITIVE







## 6) Interesting plots
def generate_networkx_graph(I,constants,variables):
    names = list(variables)+list(constants)
    G=nx.DiGraph()
    G.add_nodes_from(names)
    G.add_edges_from([(names[I[i][j]],names[i]) for i in range(len(variables)) for j in range(len(I[i]))])
    return G #plot this using nx.draw(G)

def generate_networkx_graph_from_edges(I,n_variables):
    edges = []
    for j,regulators in enumerate(I):
        if j>=n_variables: #exclude constant self-loops
            break
        for i in regulators:
            edges.append((i,j))
    return nx.DiGraph(edges) #plot this using nx.draw(generate_networkx_graph_from_edges(I,n_variables))

def bias_vs_sens(n,nsim_per_k=1000,avg_sens_sim=1000): 
    assert n>1
    f,ax = plt.subplots()
    colorlist = np.array(list(set(mcolors.BASE_COLORS.keys())-set('w')))
    biass,senss = [],[]
    ks = (list(range(n-1)) + [n])
    for k in ks:
        biass.append([])
        senss.append([])
        for i in range(nsim_per_k):
            F = random_k_canalizing(n,k,True)#random_non_degenerated_function(n)
            biass[-1].append(absolute_bias(F,n))
            senss[-1].append(average_sensitivity(F,avg_sens_sim))
        ax.plot(biass[-1],senss[-1],'o',color=colorlist[k],label=str(k),alpha=0.4)
    ax.legend(loc='best',title='canalizing depth')
    ax.set_xlabel('absolute bias')
    ax.set_ylabel('average sensitivity')
    ax.set_title('Boolean functions in %i essential variables' % n)
    plt.savefig('bias_vs_sens_v01_n%i_nsim%i_avgsenssim%i.pdf' % (n,nsim_per_k,avg_sens_sim))

def bias_vs_sens_v2(n,nsim_per_k=1000,avg_sens_sim=1000): 
    assert n>1
    f,ax = plt.subplots()
    colorlist = np.array(list(set(mcolors.BASE_COLORS.keys())-set('w'))+['orange'])
    biass,senss = [],[]
    ks = (list(range(n-1)) + [n])
    for k in ks:
        biass.append([])
        senss.append([])
        for i in range(nsim_per_k):
            F = random_k_canalizing(n,k,True)#random_non_degenerated_function(n)
            biass[-1].append(absolute_bias(F,n))
            senss[-1].append(average_sensitivity(F,avg_sens_sim))
        which_biases = list(set(biass[-1]))
        which_biases.sort()
        median_sens = [np.median(np.array(senss[-1])[np.array(biass[-1]) == el]) for el in which_biases]
        ax.plot(which_biases,median_sens,'x-',color=colorlist[k],label=str(k),alpha=0.4)
    ax.legend(loc='best',title='canalizing depth')
    ax.set_xlabel('absolute bias')
    ax.set_ylabel('median average sensitivity')
    ax.set_title('Boolean functions in %i essential variables' % n)
    f.savefig('bias_vs_sens_v02_n%i_nsim%i_avgsenssim%i.pdf' % (n,nsim_per_k,avg_sens_sim))
    
    g,ax = plt.subplots()
    ax.boxplot(senss)
    ax.set_xticklabels(ks)
    ax.set_xlabel('canalizing depth')
    ax.set_ylabel('average sensitivity')
    ax.set_title('Boolean functions in %i essential variables' % n)
    g.savefig('average_sensitivity_v02_n%i_nsim%i_avgsenssim%i.pdf' % (n,nsim_per_k,avg_sens_sim))
    
