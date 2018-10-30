#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
from collections import Counter
from tabulate import tabulate
 
 
def Huckel(H):
    print()
    print('Energies and Degeneracies of ' + u'\u03A0'.lower() + '-System from Huckel Theory')
    print()
 
    evals,evacs = np.linalg.eig(H)
    Decimals = Counter(np.round(np.array(evals),3)).items()
    print(tabulate(sorted(Decimals), headers=('Energy','Degeneracy')))
 
     
def linear_huckel(n):
    H = np.zeros((n,n))
    for i in range(0,n-1):
        H[i,i+1] = -1
        for j in range(0,n):
            H[j,i] = H[i,j]
    Huckel(H)
 
 
def cyclic_huckel(n):
    H = np.zeros((n,n))
    for i in range(0,n-1):
        H[i,i+1] = -1
        H[0,n-1] = -1
        for j in range(0,n):
            H[j,i] = H[i,j]
    Huckel(H)
 
 
def platonic_huckel(n):
    if n == 4:
        H = -np.ones((n,n))
        for i in range(0,n):
            H[i,i] = 0
        Huckel(H)
 
    elif n == 8:
        H = np.mat("0,-1,-1,0,-1,0,0,0;-1,0,0,-1,0,-1,0,0;-1,0,0,-1,0,0,1,0;0,-1,-1,0,0,0,0,-1;-1,0,0,0,0,-1,-1,0;0,-1,0,0,-1,0,0,-1;0,0,-1,0,-1,0,0,-1;0,0,0,-1,0,-1,-1,0")
        Huckel(H)
         
    elif n == 20:
        H = np.mat("0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,-1,0,0,0,0;0,0,0,0,-1,-1,0,0,0,0,0,0,-1,0,0,0,0,0,0,0;0,0,0,0,0,0,-1,0,0,0,0,0,0,-1,0,0,0,0,-1,0;0,0,0,0,0,0,0,-1,0,0,0,0,0,0,-1,0,0,0,0,-1;0,-1,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,-1,0;0,-1,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,-1;0,0,-1,0,0,0,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0;0,0,0,-1,0,0,0,0,0,0,0,-1,0,0,0,-1,0,0,0,0;0,0,0,0,0,0,0,0,0,-1,0,0,0,-1,0,0,-1,0,0,0;0,0,0,0,0,0,0,0,-1,0,0,0,0,0,-1,0,0,-1,0,0;0,0,0,0,-1,0,-1,0,0,0,0,-1,0,0,0,0,0,0,0,0;0,0,0,0,0,-1,0,-1,0,0,-1,0,0,0,0,0,0,0,0,0;0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,0,0;-1,0,-1,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0;-1,0,0,-1,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0;-1,0,0,0,0,0,-1,-1,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,-1,0,0,0,-1,0,0,0,0,0,-1,0;0,0,0,0,0,0,0,0,0,-1,0,0,-1,0,0,0,0,0,0,-1;0,0,-1,0,-1,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0;0,0,0,-1,0,-1,0,0,0,0,0,0,0,0,0,0,0,-1,0,0")
        Huckel(H)
         
    else:
        print('This does not define a Platonic solid (which exists as a hydrocarbon anyway).')
        
def buckminster_huckel(n): 
        H = np.zeros((60,60))
        H[0,4] = -1
        H[0,8] = -1
        H[1,11] = -1
        H[2,14] = -1
        H[3,17] = -1
        H[5,19] = -1
        H[6,21] = -1
        H[7,24] = -1
        H[9,25] = -1
        H[10,28] = -1
        H[12,29] = -1
        H[13,32] = -1
        H[15,33] = -1
        H[16,36] = -1
        H[18,37] = -1
        H[20,39] = -1
        H[22,41] = -1
        H[23,43] = -1
        H[26,44] = -1
        H[27,46] = -1
        H[30,47] = -1
        H[31,49] = -1
        H[34,50] = -1
        H[35,52] = -1
        H[38,53] = -1
        H[40,54] = -1
        H[42,56] = -1
        H[45,57] = -1
        H[48,58] = -1
        H[51,59] = -1
        H[55,59] = -1
        for i in range(0,n-1):
            H[i,i+1] = -1
            for j in range(0,n):
                H[j,i] = H[i,j]
        Huckel(H)
        
         
def solver():
 
    print()
    print('Is the molecule a linear polyene, a cyclic polyene, a Platonic solid or the Buckmister Fullerene?')
    t = input('Please enter linear, cyclic, Platonic or Buckmister Fullerene: ').lower()
 
    if t == 'linear':
        print()
        print('How many C atoms are there?')
        n = int(input('Please enter the number of C atoms: '))
        linear_huckel(n)
 
    elif t == 'cyclic':
        print()
        print('How many C atoms are there?')
        n = int(input('Please enter the number of C atoms: '))
        cyclic_huckel(n)
 
    elif t == 'platonic':
        print()
        print('How many C atoms are there?')
        n = int(input('Please enter the number of C atoms: '))
        platonic_huckel(n)
        
    elif t == 'buckminster fullerene':
        print()
        buckminster_huckel(60)
 
    else:
        print()
        print('Oops! Look\'s like that isn\'t a valid system for a hydrocarbon. Maybe try again?')
        print()
        solver()
 
    print()
    print('Would you like to solve for another molecule?')
 
    r = input('Please enter either yah or nah: ').lower()
 
    if r == 'yah':
        solver()
    else:
        print('Thank you for using my Marvellous Huckel Calculator.')
 
 
print('Ever had trouble calculating the energy levels of a ' + u'\u03A0'.lower() + 'system and wanted an easier method? \n Well you came to the right place my friend! This is the Marvellous Huckel Calculator')
solver()

