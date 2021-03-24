#!/usr/bin/env python
"""
2021-02-17
Karn Parmar
an organic chemist

This script allows you to take a compound containing an aromatic ring and manipulate the xyz coordinates 
such that it is placed in the XZ plane. A ghost atom is then placed at a user chosen position. The
script attempts to find all aromatic rings automatically and only fails at locating central aromatic
rings or ones containing trans double bonds.

Note: Since the script begins with a random atom in the input coordinates it is not necessarily the case
that the script will output the same XYZ coordinates when the same file is input.

This script is beneficial since it can be used with the command line, allowing for a high throughput
of structures. Rings may be chosen manually in cases where rings are not located automatically.

Note: Atom numbers start at "0"

The script is not 100% complete and has a lot of untrimmed/useless code. However, it is usable in its current form.

Have fun

=)


"""



import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import math
from operator import itemgetter
import random

x = [0, 0, 0]
y = [0, 1, 0]
z = [0, 0 ,1]

"""
Things left to do:
[X]manual input for file
[X]automatic input for atomcount
[X]automatic search for rings
[X]request user input for which ring is placed in XY plane
[X]various types of NICS calculations
[X]place Bq's
[ ]give option to place Bq's on BOTH sides of ring
[ ]create .txt savefile with Bq #'s that correspond with a ring and its coordinates
[ ]create a logfile or save coordinates to an .xyz file

""" 
 
import numpy as np
import math
import sys
import os

# TEMP VARIABLE ASSIGNMENTS
TEMPX = 3
TEMPY = 26


test = "testing"

coord = []

specifyring = False
globalatom1 = 0
globalatom2 = 1
globalatom3 = 2

filetext = []
#filename = '//home//ksp325//UTILS//PYTHON//testing//ANNULENE.xyz'
#filename = os.getcwd()+ '/thiophene.xyz'

npycoord = np.array([])
npycentered = np.array([])
atomlist = []
globfinalarr = []
glob_rings = []
rings_for_file = []
bqs_for_file = []

def rotate_mol(angle, center, nullaxis):
    #test rotation on xy axis
    if(nullaxis == 1):
        for i in range(0, len(npycoord)):
            newx = math.cos(float(math.radians(angle)))*(float(npycoord[i][0])) - math.sin(float(math.radians(angle)))*(float(npycoord[i][2]))
            newz = math.cos(float(math.radians(angle)))*(float(npycoord[i][2])) + math.sin(float(math.radians(angle)))*(float(npycoord[i][0]))
            npycoord[i][0] = float(newx)
            npycoord[i][2] = float(newz)
    if(nullaxis == 2):
        for i in range(0, len(npycoord)):
            newy = math.cos(float(math.radians(angle)))*(float(npycoord[i][1])) - math.sin(float(math.radians(angle)))*(float(npycoord[i][0]))
            newx = math.cos(float(math.radians(angle)))*(float(npycoord[i][0])) + math.sin(float(math.radians(angle)))*(float(npycoord[i][1]))
            npycoord[i][1] = float(newy)
            npycoord[i][0] = float(newx)
    if(nullaxis == 0):
        for i in range(0, len(npycoord)):
            newy = math.cos(float(math.radians(angle)))*(float(npycoord[i][1])) - math.sin(float(math.radians(angle)))*(float(npycoord[i][2]))
            newz = math.cos(float(math.radians(angle)))*(float(npycoord[i][2])) + math.sin(float(math.radians(angle)))*(float(npycoord[i][1]))
            npycoord[i][1] = float(newy)
            npycoord[i][2] = float(newz)               

def move_center(central_point):
    deltaX = 0
    deltaY = 0
    deltaZ = 0
    if (float(npycoord[central_point][0]) > 0):
        deltaX = -float(npycoord[central_point][0])
    elif (float(npycoord[central_point][0]) < 0):
        deltaX = abs(float(npycoord[central_point][0]))
    else: 
        deltaX = 0
    if (float(npycoord[central_point][1]) > 0):
        deltaY = -float(npycoord[central_point][1])
    elif (float(npycoord[central_point][1]) < 0):
        deltaY = abs(float(npycoord[central_point][1]))
    else: 
        deltaY = 0
    if (float(npycoord[central_point][2]) > 0):
        deltaZ = -float(npycoord[central_point][2])
    elif (float(npycoord[central_point][2]) < 0):
        deltaZ = abs(float(npycoord[central_point][2]))
    else: 
        deltaZ = 0    
    for i in range(0, len(npycoord)):
        npycoord[i][0] = float(npycoord[i][0]) + deltaX
        npycoord[i][1] = float(npycoord[i][1]) + float(deltaY)
        npycoord[i][2] = float(npycoord[i][2]) + deltaZ
        
    newarr = "\n".join(map(str, npycoord))


    return npycoord


#Extract all information from a given string in a file, place all text between spaces in an array
def slice_string(string1):
    tempstr = '' 	 	
    finalarr = []
    whitespace = []
	
# Find locations of all whitespaces and record values
    whitespace.append(0) # add a white space for the start
    for i in range(0,len(string1)):
        if(string1[i] == ' '):
            whitespace.append(i+1)  ##MODIFICATION (+1)
    whitespace.append(len(string1))
    print(whitespace)

# Append all values inbetween whitespaces into an array
# the number of whitespaces determines the number of iterations
    for j in range(0,len(whitespace)-1):
        for k in range(whitespace[j],whitespace[j+1]):
            tempstr += string1[k]
        finalarr.append(tempstr) # Append the string to the array
        tempstr = "" # Clear variable after use
        
# Delete all the extra white spaces in the array
    while (' ' in finalarr):
        finalarr.remove(' ')
    if('' in finalarr):
        finalarr.remove('')

    counter=0
    for i in range(0, len(finalarr) - 1):
        finalarr[i] = finalarr[i].replace(" ", "")
        finalarr[i] = finalarr[i].replace('\n', "")
    atomlist.append(finalarr[0])
    return finalarr



def getAngle(a, b, c):
    ang = math.degrees(math.atan2(c[1]-b[1], c[0]-b[0]) - math.atan2(a[1]-b[1], a[0]-b[0]))
    return float(ang)




def locate_ring():

    global specifyring
    global globalatom1
    global globalatom2
    global globalatom3
    global glob_rings
    global param_choosering_now
    global param_choosering_atoms
    global param_choosering
    #### DEBUGGING
    #print(glob_rings)
    #print(glob_rings[0][0][0])
    atom1 = glob_rings[0][0][0]
    atom2 = glob_rings[0][1][0]
    atom3 = glob_rings[0][2][0]
    chooseringatomfound = 0
    
    if(param_choosering_now == 1):
        for i in range(0, len(glob_rings)):
            for j in range(0, len(glob_rings[i])):
                if(glob_rings[i][j][0] == int(param_choosering_atoms[0])):
                    print("......Chose a ring!")
                    atom1 = int(param_choosering_atoms[0]) + 1 
                    for k in range(0, len(glob_rings[i])):
                        if (atom1 is not glob_rings[i][k][0]):
                            atom2 = glob_rings[i][k][0]
                    for k in range(0, len(glob_rings[i])):
                        if (atom1 is not glob_rings[i][k][0] and atom2 is not glob_rings[i][k][0]):
                            atom3 = glob_rings[i][k][0]                    
                    #atom2 = glob_rings[i][0][0]
                    #atom3 = glob_rings[i][1][0]
                    #print(atom1, atom2, atom3)
                    chooseringatomfound = 1
        if(chooseringatomfound != 1):
            print(".........ERROR: ATOM NOT IN RING")

    elif(specifyring == True):
        atom1 = globalatom1
        atom2 = globalatom2
        atom3 = globalatom3
    else:    
        atom1 = glob_rings[0][0][0]
        atom2 = glob_rings[0][1][0]
        atom3 = glob_rings[0][2][0]
    
    specifyring = False

    sliceinput = [atom1, atom2, atom3]
    move_center(int(sliceinput[0]))
    anglezx = getAngle((float(npycoord[int(sliceinput[1])][0]), float(npycoord[int(sliceinput[1])][2])), (float(0.0), float(0.0)), (10.0, 0.0))
    rotate_mol(anglezx, npycoord[0], 1)
    anglexy = getAngle((float(npycoord[int(sliceinput[1])][0]), float(npycoord[int(sliceinput[1])][1])), (float(0.0), float(0.0)), (10.0, 0.0))
    rotate_mol(-anglexy, npycoord[0], 2)
    angleyz = getAngle((float(npycoord[int(sliceinput[2])][1]), float(npycoord[int(sliceinput[2])][2])), (float(0.0), float(0.0)), (0.0, 10.0))
    rotate_mol(angleyz, npycoord[0], 0)
    
    tempringforfile = []
    for i in range(0, len(glob_rings)):
        tempringforfile = []
        for j in range(0, len(glob_rings[i])):
            tempringforfile.append(glob_rings[i][j][0])
        rings_for_file.append(tempringforfile)

    return sliceinput

def displaygraph():
    
    global atomlist
    xcoordinates = []
    ycoordinates = []
    zcoordinates = []
    # x axis values 
    for i in range(0, len(npycoord)):
        xcoordinates.append(float(npycoord[i][0]))
    for i in range(0, len(npycoord)):
        ycoordinates.append(float(npycoord[i][1]))
    for i in range(0, len(npycoord)):
        zcoordinates.append(float(npycoord[i][2]))
    #TESTING FEATURE:    
    #matplotlib.use('Agg')
    x = xcoordinates
    # corresponding y axis values 
    y = ycoordinates 
    z = zcoordinates
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    n = 100
    counter = 0
    xs = []
    ys = []
    zs = []
    for i in range(0, len(xcoordinates)):
        xs = x[i]
    for i in range(0, len(ycoordinates)):
        ys = y[i]
    for i in range(0, len(zcoordinates)):
        zs = z[i]

    foundatom = 0
    foundbq = 0
    carbonatoms = []
    bqatoms = []
    for i in range(0, len(atomlist)):
        if (atomlist[i] != 'H' and atomlist[i] != 'Bq'):
            carbonatoms.append(i)
        if (atomlist[i] == 'Bq'):
            bqatoms.append(i)
    
    for i in range(0, len(xcoordinates)):
        foundatom = 0

        for j in range(0, len(carbonatoms)):

            if (i == int(carbonatoms[j])):
                #ax.scatter(x[i], y[i], z[i], c='r', marker='^')
                foundatom = 1
                break;
        for j in range(0, len(bqatoms)):
            if (i == int(bqatoms[j])):
                foundbq = 1
                break;
                
        if(foundatom == 1):
            ax.scatter(x[i], y[i], z[i], c='r', marker='o')
        elif(foundatom == 0 and foundbq == 1):
            ax.scatter(x[i], y[i], z[i], c='g', marker='^')
        else:
            ax.scatter(x[i], y[i], z[i], c='b', marker='.')
            

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    sortedx = sorted(xcoordinates)
    axessort = []
    difference = float(sortedx[0])-float(sortedx[-1])
    axessort.append(difference)
    sortedy = sorted(ycoordinates)
    difference = float(sortedy[0])-float(sortedy[-1])
    axessort.append(difference)
    sortedz = sorted(zcoordinates)
    difference = float(sortedz[0])-float(sortedz[-1])
    axessort.append(difference)
    
    sorted(axessort)

    
    ax.set_xbound(float(-axessort[-1]), float(axessort[-1]))
    ax.set_ybound(float(-axessort[-1]), float(axessort[-1]))
    ax.set_zbound(float(-axessort[-1]), float(axessort[-1]))
    
    plt.show()
    
def autosearchring():

    carbonlist = []
    trimmedcarbonlist = [] #contains only carbon atoms that are bonded to two others
    trimmedcarbonbonds = []
    tempbondlist = []
    tempring = [] # holds carbon numbers of potential rings
    counter2 = 0
    counter = 0
    #find first carbon atom from beginning
    for i in range(0, len(atomlist)):
        counter = 0
        if(atomlist[i] == 'C' or atomlist[i] == 'S' or atomlist[i] == 'O' or atomlist[i] == 'N' or atomlist[i] == 'B' or atomlist[i] == 'P'):
            for j in range(0, len(atomlist)): ## REMOVED -1 <<len(atomlist-1>> on 2020-11-16
                if(atomlist[i] == 'S' or atomlist[j] == 'S' or atomlist[i] == 'P' or atomlist[j] == 'P' or atomlist[i] == 'B' or atomlist[j] == 'B'):
                    #print("FOUND_SULFUR_ETC")
                    if(distance_finder(i, j) < 1.8 ):
                        counter = counter +  1
                else:
                    if(distance_finder(i, j) < 1.6):
                        counter = counter +  1                
            if(counter < 5):
                carbonlist.append(i)
    #print("CARBONLIST", carbonlist)            # DEBUGGING
    #trim carbon list to include only carbons with two others bonded
    for i in range(0, len(carbonlist)):
        tempval = carbonlist[i]
        tempbondlist = []
        counter = 0
        for j in range(0, len(carbonlist)):
            tempval2 = carbonlist[j]
            if(atomlist[tempval] == 'S' or atomlist[tempval] == 'N' or atomlist[tempval] == 'O' or atomlist[tempval] == 'P' or atomlist[tempval] == 'B'):

                if (distance_finder(carbonlist[i], carbonlist[j]) < 1.8):

                    if(i is not j):
                        counter = counter + 1
                        tempbondlist.append(carbonlist[j])        #make an array containing all bonded atoms  

            else:
                if(atomlist[tempval2] == 'S' or atomlist[tempval2] == 'P'):
                    if (distance_finder(carbonlist[i], carbonlist[j]) < 1.8):
                        if(i is not j):
                            counter = counter + 1
                            tempbondlist.append(carbonlist[j]) 

                else:             
                    if (distance_finder(carbonlist[i], carbonlist[j]) < 1.6):
                        if(i is not j):
                            counter = counter + 1
                            tempbondlist.append(carbonlist[j])                
                                    
        if (counter>1 and counter<4):
            trimmedcarbonlist.append(carbonlist[i])
            trimmedcarbonbonds.append(tempbondlist)
    
    
    todelete = []
    noadd = True
    duplicatetrimmedcarbonlist = trimmedcarbonlist
    duplicatetrimmedcarbonbonds = trimmedcarbonbonds
    #print("TRIMMEDCARBONLIST", trimmedcarbonlist) #DEBUGGING
    #print("TRIMMEDCARBONBONDS", trimmedcarbonbonds) #DEBUGGING
    #[ ] Items in trimmedcarbonlist must have two atoms/heteroatoms bonded
    # WHILE (INITIAL != FINAL):
        #[x] GET INITIAL ARRAY
        #[x] Delete items in trimmedcarbonbonds that are NOT in trimmed carbonlist
            #[x] Clear array elements that are blank
        #[ ] Delete items in trimmedcarbonlist that are not bonded to two others in the list
        #[ ] GET FINAL ARRAY
    
    trimmedcarbonbondsloop = []
    trimmedcarbonlistloop = []
    inlist = False
    
  
    for z in range(0, 100):
        for m in range(0, len(trimmedcarbonbonds)):
            counter = 0
            for i in range(0, len(trimmedcarbonbonds[m])):
                counter = 0
                for j in range(0, len(trimmedcarbonlist)):
                    if(trimmedcarbonbonds[m][i] == trimmedcarbonlist[j]):
                            inlist = True
                if (inlist == False):
                    duplicatetrimmedcarbonbonds[m][i] = '-1'      
                inlist = False
        for i in range(0, len(trimmedcarbonbonds)):        
            try:       
                duplicatetrimmedcarbonbonds[i].remove('-1')
            except ValueError:
                pass
        trimmedcarbonlist = duplicatetrimmedcarbonlist
        trimmedcarbonbonds = duplicatetrimmedcarbonbonds  
        for m in range(0, len(trimmedcarbonbonds)): # clear array elements that are blank/have less than 2 bonds
            if (len(trimmedcarbonbonds[m])< 2):
                duplicatetrimmedcarbonbonds[m] = '-1'
                duplicatetrimmedcarbonlist[m] = '-1'
        for i in range(0, len(trimmedcarbonbonds)):        
            try:       
                duplicatetrimmedcarbonbonds.remove('-1')

            except ValueError:
                pass      
        for i in range(0, len(trimmedcarbonlist)):        
            try:       
                duplicatetrimmedcarbonlist.remove('-1')
            except ValueError:
                pass  
        duplicatetrimmedcarbonlist = trimmedcarbonlist
        duplicatetrimmedcarbonbonds = trimmedcarbonbonds       
 
        counter = 0
        #print(len(trimmedcarbonlist))
        #print(len(trimmedcarbonbonds))
        #print(trimmedcarbonbonds)
        #print(trimmedcarbonlist)
        for m in range(0, len(trimmedcarbonlist)-1): # delete all items containing less than 2 bonds to others in the list
            for i in range(0, len(trimmedcarbonbonds)):
                for j in range(0, len(trimmedcarbonbonds[i])):
                     #print(m, i, j, "[", trimmedcarbonlist[m], "] [", trimmedcarbonbonds[i][j], "]")
                     #print(trimmedcarbonlist.index(trimmedcarbonlist[m]), trimmedcarbonbonds.index(trimmedcarbonbonds[i][j]))
                     #print("trimmedcarbonbonds[i]", trimmedcarbonbonds[i]) #DEBUGGING
                     #print("trimmedcarbonlist[m]", trimmedcarbonlist[m])   #DEBUGGING    
                     #print(trimmedcarbonlist[m], trimmedcarbonbonds[i][j])     
                     try:         
                         if (trimmedcarbonlist[m] == trimmedcarbonbonds[i][j]):
                             counter = counter + 1
                     except:
                         pass 
            if (counter <2):
                try:
                    duplicatetrimmedcarbonlist.pop(m)   
                    duplicatetrimmedcarbonbonds.pop(m)
                except:
                    pass
            counter = 0

    trimmedcarbonbonds = duplicatetrimmedcarbonbonds
    trimmedcarbonlist = duplicatetrimmedcarbonlist
    
    isinring = 0
    potentialring = []
    minangle = []
    ring = []
    rings = []
    ringadd = []
    foundnew = False
    ## NEW RING FINDING ALGORITHM ##
    #print("\n\n\nRING SEARCH:\n\n")
    ringfound = False
    ringsearchlist = list(trimmedcarbonlist)
    ringsearchlistcounter = 0
    while(len(ringsearchlist) > 0):
        #print("START_WHILE_LOOP")
        whileloopcount = 0
        ringsearchlistcounter += 1
        i_rand = random.choice(ringsearchlist)
        i_counter = 0   
        for i in range(0, len(trimmedcarbonlist)): #changed from ringsearchlist 
            if(i_rand == trimmedcarbonlist[i]):
                i_counter = i
        #add the first atom and bonds
        ringadd = []
        ring = []
        ringadd = [trimmedcarbonlist[i_counter],trimmedcarbonbonds[i_counter][0]]  # Add the first two carbon atoms (the bond)
        ring.append(ringadd)
        for i in range(0, len(trimmedcarbonlist)):
            if(ringadd[1] == trimmedcarbonlist[i]):
                i_counter = i

        ringfound = False
        ringfoundcounter = 0
        while(ringfound == False):
            ringfoundcounter += 1
            if(ringfoundcounter > 500):
                ringfound = True

            for i in range(0, len(trimmedcarbonlist)):                             # Search for the second carbon atom in the bond in the TCL
                try:
                    if(trimmedcarbonlist[i] == ring[-1][1]):

                        i_counter = i
                        foundnew = True
                except ValueError:
                    foundnew = True
                    ringfound = True

            j_counter = 0
            halting = False
            whileloopcount = 0
            while (halting == False):
                whileloopcount += 1
                j = random.choice(trimmedcarbonbonds[i_counter])
                if(j != ring[-1][0] and trimmedcarbonbonds[i_counter] != ring[-1]):
                    if(len(ring)>3):
                        if(measureangle(ring[-4][0], ring[-3][0], ring[-2][0])> (measureangle(ring[-4][0], ring[-3][0], ring[-1][0]) + 1)):
                            #print(trimmedcarbonbonds[i_counter][j_counter])
                            ringadd = [trimmedcarbonlist[i_counter],j]  # Add this new bond to a temporary array 
                            j_counter = len(trimmedcarbonbonds[i_counter])   
                            ring.append(ringadd)   # Add the temporary array to the ring array
                            halting = True
                        else:
                            #print("ring before del", ring)
                            del ring[-1]
                            del ring[-1]
                            if(len(ring) < 4):
                                del ring[-1]
                            if(len(ring) == 2):
                                del ring[-1]
                            #print("ring after del", ring)
                            #print("Deleting ring[-1]")
                            if (len(ring) < 1):
                                ringfound = True
                               
                            j_counter = len(trimmedcarbonbonds[i_counter])
                            halting = True
                    else:   
                        #print("------------else: ringadd/j_counter")
                        ringadd = [trimmedcarbonlist[i_counter],j]  # Add this new bond to a temporary array   
                        j_counter = len(trimmedcarbonbonds[i_counter])
                        ring.append(ringadd)     # Add the temporary array to the ring array
                        halting = True
                j_counter = j_counter + 1     
                if(whileloopcount > 20): # clear the entire ring array
                    #print("\n------------whileloopcount>20")
                    rings = []
                    ring = []
                    halting = True
                    ringfound = True
                    break;

            ringadd = []
            # Evaluate whether or not a ring has been located:
            ringpositions = -1
            tempring = []
            trimrings = False
            for i in range(0, len(ring)):
                if (ring[-1] == ring[i]):
                    if ((len(ring) - i) > 3):
                        trimrings = True
                        ringpositions = i
                        #print("FOUND A RING")
                        #print("STRING OF BONDS BEFORE TRIMMING", ring)
            if(trimrings == True):
                tempring = ring[ringpositions:-1]
                ring = tempring
                ringfound = True
            # Delete all items in list not found in the ring
        isring = False
        for i in range(0, len(ring)):
            if(ring[-1][1] == ring[i][0]):
                isring = True
        if(isring == True):    
            for i in range(0, len(ring)):
                for j in range(0, len(ringsearchlist)):
                    if (ringsearchlist[j] == ring[i][0]):
                        ringsearchlist[j] = '-1'
        """for i in range(0, len(trimmedcarbonlist)):
            try: 
                ringsearchlist.remove('-1')  
            except ValueError:
                pass """
        rings.append(ring)
        if(ringsearchlistcounter > atomcount*2):
            break;
            ringsearchlist[i_rand] = '-1'
        for i in range(0, len(trimmedcarbonlist)):
            try: 
                ringsearchlist.remove('-1')  
            except ValueError:
                pass         
        #print("ringsearchlist", ringsearchlist)
        #print("RINGS", rings)
        #print("---------\n")
    global glob_rings
    glob_rings = rings
    temp_glob_rings = glob_rings
    
    rings_to_delete = []
    ## The following code prevents non-rings from being placed into glob_rings)
    ## There is a current error in the ring finding algorithm and strings of four bonds at ring junctions 
    ##    are sometimes identified as rings 
    ## This ensures all items in glob_rings are actually rings!
    ## This does not affect rings that are manually input
    for i in range(0, len(glob_rings)):
        templength = len(glob_rings[i])
        if(glob_rings[i][0][0] is not glob_rings[i][templength-1][1]):
            temp_glob_rings = glob_rings
            rings_to_delete.append(i) 
    rings_to_delete.sort(reverse=True)
    for i in range(0, len(rings_to_delete)):
        del temp_glob_rings[rings_to_delete[i]]
    glob_rings = temp_glob_rings
    ## End glob_rings clean up
    
    ## OPTIONS THAT MODIFY THE rings/glob_rings VARIABLES
    if(param_manual == 1):
        ring = []
        rings = []
        glob_rings = []
        temp_mra = []
        temp = []
        for j in range(0, len(manual_ring_atoms)):
            for i in range(0, len(manual_ring_atoms[j])):
                temp_mra = manual_ring_atoms[j][i]
                temp = [temp_mra, 0] # add the zero to make the array fit with the rest of the program!
                rings.append(temp)             
            glob_rings.append(rings)
            temp_mra = []
            temp = []
            rings = []
            print("-------------PARAM_MANUAL-----------")
            print(glob_rings)
    if(param_automanual == 1):
        rings = []
        for i in range(0, len(automanual_ring_atoms)):
            temp = [int(automanual_ring_atoms[i]), 0] # add the zero to make the array fit with the rest of the program!
            rings.append(temp)
        glob_rings.append(rings)          

    
def measureangle(atom1, atom2, atom3):

    a = [float(npycoord[atom1][0]),float(npycoord[atom1][1]),float(npycoord[atom1][2])]
    b = [float(npycoord[atom2][0]),float(npycoord[atom2][1]),float(npycoord[atom2][2])]
    c = [float(npycoord[atom3][0]),float(npycoord[atom3][1]),float(npycoord[atom3][2])]

    ba = [float((a[0])-float(b[0])), (float(a[1])-float(b[1])), (float(a[2])-float(b[2]))]
    ac = [float((c[0])-float(b[0])), (float(c[1])-float(b[1])), (float(c[2])-float(b[2]))]

    cosine_angle = np.dot(ba, ac) / (np.linalg.norm(ba) * np.linalg.norm(ac))
    angle = np.arccos(cosine_angle)
    
    return(np.degrees(angle))

def distance_finder(atom1,atom2) :

    return (((float(npycoord[atom2][0])-float(npycoord[atom1][0]))**2)+((float(npycoord[atom2][1])-float(npycoord[atom1][1]))**2)+((float(npycoord[atom2][2])-float(npycoord[atom1][2]))**2))**(1/2)

def prettyprint():
    print("...Writing file")
    global npyatomcoord
    npyatomcoord = np.insert(npycoord, 0, "", axis = 1)
    for i in range(0, len(npyatomcoord)):
        npyatomcoord[i][0] = atomlist[i]
        for j in range(1, 4):
            npyatomcoord[i][j] = round(float(npyatomcoord[i][j]), 5)
            if ((float(npyatomcoord[i][j]) < 0.00005) and (float(npyatomcoord[i][j]) > -0.00005)):
                npyatomcoord[i][j] = 0.0
    print("\n\ntransformed XYZ coordinates with Bq:")
    for i in range(0, len(npyatomcoord)):
        print(npyatomcoord[i][0], " ", npyatomcoord[i][1], " ", npyatomcoord[i][2], " ", npyatomcoord[i][3])
    tempxyzfilename = sys.argv[1]
    for i in range(0, int(len(tempxyzfilename)/3)):
        tempxyzfilename = tempxyzfilename.strip('.xyz')
        tempxyzfilename = tempxyzfilename.strip('.log')    
    new_file_path = os.path.join(os.getcwd(), tempxyzfilename + "_transformed.xyz")
    comfile = open(new_file_path, 'w')    
    for i in range(0, len(npyatomcoord)):
        temp = str(npyatomcoord[i][0]) + "    " + str(npyatomcoord[i][1]) +  "     " + str(npyatomcoord[i][2]) + "     " + str(npyatomcoord[i][3]) + "\n"
        comfile.write(temp)



def find_centroid():
    bq_atom_count = atomcount
    global bqs_for_file
    global double_bq

    #add bq's to the atomlist

    for i in range(0, len(glob_rings)):
        atomlist.append('Bq')
    if(double_bq == 1):
        print("-d is currently unavailable!")
    #find average x
    sum_x = 0
    average_x = 0
    sum_y = 0
    average_y = 0
    sum_z = 0
    average_z = 0
    global npyatomcoord
    global npycoord
    global op_x
    global specifyring
    global globalatom1 
    global globalatom2 
    global globalatom3 
    
    npyatomcoord = np.insert(npycoord, 0, "", axis = 1) # necessary to include otherwise prettyprint() must be run before find_centroid()
    coords_all = [[]]
    coords = []
    tempnpycoord = np.array([])
    for i in range(0, len(glob_rings)): 
        sum_x = 0
        average_x = 0
        sum_y = 0
        average_y = 0
        sum_z = 0
        average_z = 0
        coords = []
        tempnpycoord = npycoord

        # Rotate the molecule to place each ring in plane to ensure out of plane Bq's are placed properly
        specifyring = True
        globalatom1 = glob_rings[i][0][0]
        globalatom2 = glob_rings[i][1][0]
        globalatom3 = glob_rings[i][2][0]
        locate_ring()
        specifyring = False

        for j in range(0, len(glob_rings[i])):
            location = glob_rings[i][j][0]
            sum_x = float(sum_x) + float(npycoord[location][0]) 
            sum_y += float(npycoord[location][1]) 
            sum_z += float(npycoord[location][2]) 
            average_x = sum_x/float(len(glob_rings[i]))
            average_y = sum_y/float(len(glob_rings[i]))
            average_z = sum_z/float(len(glob_rings[i]))
        #bqs_for_file[i] = bq_atom_count
        #bq_atom_count += 1
        NICS_ZERO_XYZ = []
        NICS = []
        NICS_ZERO_XYZ.append(average_x)
        if(double_bq == 1):
            if (op_zero == 0 and op_vert == 0 and op_x == 0):
                NICS_ZERO_XYZ.append(average_y - 1)
    
            elif (op_x == 1):
                NICS_ZERO_XYZ.append(average_y + (-1*float(op_x_value[0])))
        else:        
            if (op_zero == 1 or op_vert == 1):
                    NICS_ZERO_XYZ.append(average_y)
            elif (op_zero == 0 and op_vert == 0 and op_x == 0):
                       NICS_ZERO_XYZ.append(average_y + 1)
        
            elif (op_x == 1):
                NICS_ZERO_XYZ.append(average_y + float(op_x_value[0]))
                            
        NICS_ZERO_XYZ.append(average_z)
        NICS.append(NICS_ZERO_XYZ)
        npycoord = np.append(npycoord, NICS, axis = 0) 
        NICS[0].insert(0, 'Bq')
        npyatomcoord = np.append(npyatomcoord, NICS, axis = 0)        
        bq_npycoord = np.array([average_x,average_y,average_z])        
        tempnpycoord = np.append(tempnpycoord, bq_npycoord)

        
        ### TESTING TESTING
    bq_loc_sum = 0
    if(double_bq == 1):
        multiple = -1
    if(double_bq == 0):
        multiple = 1
        
    if(op_vert == 1):
        Bq_locations = multiple*(float(op_vert_locations[0])/float(op_vert_locations[1]))
        for i in range(0, int(Bq_locations)):
            for j in range(0, len(glob_rings)):
                atomlist.append('Bq')
        #find average x
        
        for k in range(0, int(Bq_locations)): 
            coords_all = [[]]
            coords = []
            tempnpycoord = npycoord
            bq_loc_sum = bq_loc_sum + float(op_vert_locations[1])
            for i in range(0, len(glob_rings)):
                sum_x = 0
                average_x = 0
                sum_y = 0
                average_y = 0
                sum_z = 0
                average_z = 0
                coords = []
                specifyring = True
                globalatom1 = glob_rings[i][0][0]
                globalatom2 = glob_rings[i][1][0]
                globalatom3 = glob_rings[i][2][0]
                locate_ring()
                specifyring = False
                for j in range(0, len(glob_rings[i])):
                    location = glob_rings[i][j][0]
                    sum_x = float(sum_x) + float(npycoord[location][0]) 
                    sum_y += float(npycoord[location][1]) 
                    sum_z += float(npycoord[location][2]) 
                    average_x = sum_x/float(len(glob_rings[i]))
                    average_y = sum_y/float(len(glob_rings[i]))
                    average_z = sum_z/float(len(glob_rings[i]))
                NICS_ZERO_XYZ = []
                NICS = []
                NICS_ZERO_XYZ.append(average_x)                
                NICS_ZERO_XYZ.append(average_y + multiple*bq_loc_sum) 
                NICS_ZERO_XYZ.append(average_z)
                NICS.append(NICS_ZERO_XYZ)
                npycoord = np.append(npycoord, NICS, axis = 0) 
                NICS[0].insert(0, 'Bq')
                npyatomcoord = np.append(npyatomcoord, NICS, axis = 0)                
                bq_npycoord = np.array([average_x,average_y,average_z])                
                tempnpycoord = np.append(tempnpycoord, bq_npycoord)
                coords.append(average_x)
                coords.append(average_y)
                coords.append(average_z)                
                coords_all.append(coords)



def menu():
    print("\n\n") 
    print("This script allows you to generate XYZ coordinates for different NICS \ncalculations by placing Bq's (ghost atoms) at various locations\n\n")
    """How this script works:")
    1. Rings are located either automatically or randomly
    2. A ring is chosen to be placed in the ZY plane
    3. A ghost atom (Bq) is placed at a specified location
    """
    print("USAGE:")
    print("----------------------------------------")
    print("nics.py [filename] [NICS options] [parameters] [other options]")
    print("----------------------------------------\n")
    print(">>>Parameters:\n")
    print("1. -m : [M]anual ring input. Must be followed by a listing of ALL ring atoms (integer input only)")
    print("\t\t\t example > nics.py [inputfile] [NICS option] -m atom1 atom2 atom3 ... \n")
    print("2. -am : [A]utomatic ring search with an additional [M]anually input ring (integer input only)\n")
    print("\t\t\t example > nics.py [inputfile] [NICS option] -am atom1 atom2 atom3 ... \n")
    print("3. -cr : [C]hoose a [R]ing to place in the ZY plane by picking one or more atoms in that ring (integer input only)")
    print("\t\t\t example > nics.py [inputfile] [NICS option] <parameters> -cr atom#")
    print("4. No parameters results in automatic ring search with a random ring selected to place in the ZY plane\n")

    print("Current known issues: \n 1. Central aromatic rings of polycyclic aromatic hydrocarbons are not found")
    print(" <<This issue can be bypassed by manually inputting rings>>\n\n")
    print(">>>NICS options: \n")
    print("1. -N0 : NICS(0),  Places Bq at centroid of all rings (including manually specified rings) ")
    print("\t\t\t example > nics.py [inputfile] -N0 [parameters]")    
    print("2. -N1 : NICS(1)/NICSzz(1), Places Bq one angstrom above ring centroid, transforms XYZ coordinates")
    print("\t\t\t example > nics.py [inputfile] -N1 [parameters]")
    print("3. -Nx : NICS(x)/NICSzz(x), Places Bq x angstrom(s) above ring centroid, transforms XYZ coordinates")
    print("\t\t\t example > nics.py [inputfile] -N1 [parameters]")
    print("4. -Nvert : Places Bq's from the centroid to a specified distance above the ring with specified interval (angstroms)")
    print("\t\t\t example > nics.py [inputfile] -Nvert <distance> <interval> [parameters]")
    print("\t\t\t example > nics.py [inputfile] -Nvert 5.0 1.0 [parameters]")
    print("\t\t\t - NOTE: The above example places 5 Bq's at 0, 1, 2, 3, 4, 5 \n\t\t\t\t angstroms above the ring.  ")
    print("\t\t\t - NOTE: Default is to place a Bqs at the centroid and up until the largest integer")
    print("\t\t\t - NOTE: if distance = 5 and interval = 3 Bq's will only be \n\t\t\t\t placed at 0 and 3 angstroms.")
    print(">>>Other options:\n")
    print("1. -g : Displays a 3D graph of the XYZ coordinates after transformation/Bq placement) (Not compatible on all terminals!) ")
    print("\t\t\t example > nics.py [inputfile] [NICS option] [parameters] -g ")
    print("2. -d <<CURRENTLY UNAVAILABLE>>: Add Bq's to both sides of a ring (useful when both sides of the ring are different)")
    print("\t\t\t example > nics.py [inputfile] [NICS options] [parameters] -d ")
    
def shortmenu():
    print("short menu")


    
## reads the user inputs (see menu() function)
def inputselect():
    global param_manual
    global manual_ring_atoms 
    global param_automanual 
    global automanual_ring_atoms
    global param_choosering
    global param_choosering_atoms
    global full_auto
    global choosering_atom
    global op_zero
    global op_one
    global op_x
    global op_x_value
    global op_vert
    global op_vert_locations
    global op_displaygraph
    global runprogram
    global double_bq

    for i in range(1, len(sys.argv)):
        if(sys.argv[i] == '-d'):
            #double_bq = 1 ## CURRENTLY UNAVAILABLE!
            print("...Double up Bq's")
    for i in range(1, len(sys.argv)):
        # search for parameters and NICS options

        # manual ring input
        if (sys.argv[i] == '-m'):
            param_manual = 1
            temp_manual = []
            for j in range(i+1, len(sys.argv)):
                temp = sys.argv[j]
                if(temp.startswith('-')): # stop searching the inputs after the next string starting with '-' is found
                    break;
                if temp[0].isdigit():
                    temp_manual.append(int(temp))
            manual_ring_atoms.append(temp_manual)
            print(temp_manual)
            print(manual_ring_atoms)
    for i in range(1, len(sys.argv)):
        if (sys.argv[i] == '-am'):
            param_automanual = 1
            for j in range(i+1, len(sys.argv)):
                temp = sys.argv[j]
                if(temp[0].startswith('-')): # stop searching the inputs after the next string starting with '-' is found
                    break;
                if temp.isdigit():
                    automanual_ring_atoms.append(temp)
    for i in range(1, len(sys.argv)):
        if (sys.argv[i] == '-cr'):
            param_choosering = 1
            for j in range(i+1, len(sys.argv)):
                temp = sys.argv[j]
                if(temp.startswith('-')): # stop searching the inputs after the next string starting with '-' is found
                    break;
                if temp.isdigit():
                    param_choosering_atoms.append(temp)

    for i in range(1, len(sys.argv)):
        if (sys.argv[i] == '-N0'):
            op_zero = 1

            break;
            
        if (sys.argv[i] == '-N1'):
            op_one = 1
            break;
        if(sys.argv[i] == '-Nx'):
            op_x = 1
            for j in range(i+1, len(sys.argv)):
                temp = sys.argv[j]
                if(temp.startswith('-')): # stop searching the inputs after the next string starting with '-' is found
                    break;
                elif temp[0].isdigit():
                    op_x_value.append(temp)
            if(len(op_x_value) < 1):
                print("\n>>Error: Did not specify a distance for -Nx !\n")
            if(len(op_x_value) > 1):
                print("Too many inputs for -Nx. Selecting the first value for distance", op_x_value[0])     
            break;
                  
        if (sys.argv[i] == '-Nvert'):
            op_vert = 1
            for j in range(i+1, len(sys.argv)):
                temp = sys.argv[j]
                if(temp.startswith('-')): # stop searching the inputs after the next string starting with '-' is found
                    break;
                elif temp[0].isdigit():
                    op_vert_locations.append(temp)
            if(len(op_vert_locations) < 2):
                print("\n>>Error: Did not specify a distance and interval!\n")
                break;
            if(len(op_vert_locations) > 2):
                print("Too many inputs for -Nvert. Selecting the first two values for distance", op_vert_locations[0], " and interval", op_vert_locations[1] ) 
                temp_op_vert_locations = [op_vert_locations[0],op_vert_locations[1]]
                op_vert_locations = temp_op_vert_locations
            print("...Placing Bqs at range of positions above centroid")
            break;
                  
    for i in range(1, len(sys.argv)):    
        if (sys.argv[i] == '-g'):
            op_displaygraph = 1
            print("...Displaying graph (blue=H, green=Bq, red=all other atoms (C,N,O,etc))")
            
     
     
          
    if (param_automanual == 1 and param_manual == 1):
        print("\n>>ERROR: Cannot select -m and -am simultaneously !\n")
        #quit()
        
    if (param_automanual == 0 and param_manual == 0):
        print("...Proceeding with automatic ring search")
    if(param_automanual == 1):
        print("...Including manually selected ring with automatically found rings")
    if(param_manual == 1):
        print("...Including ONLY manually selected ring")
    if(op_zero == 1):
        print("...Placing Bq at ring centroid (middle)")
    if(op_one == 1):
        print("...Placing Bq at one angstrom above ring centroid (middle)")
    if(op_x == 1):
        print("...Placing Bq at specified position above ring centroid (middle)")
    if(param_choosering == 1):
        print("...Placing ring containing specified atom in plane")



        #Default choice in effect, proceeding with manual ring search

#to avoid showing the menu multiple times
showmenu = 0
#in case the program should not proceed (1 runs the program, 0 does not)
runprogram = 1

try:
    filename = sys.argv[1]
except:
    menu() 
    showmenu = 1
try:
    with open(filename, 'r') as f:
        atomcount = 0
        start_at_2 = False
        for (i, line) in enumerate(f): 
            templine = line[0]
            if(templine.isdigit()):
                start_at_2 = True
            if(start_at_2 == False):
                filetext.append(line)
            else:
                if(i > 1):
                    filetext.append(line)
        #print("LINE[0]", lines[5])
        for i in range(0, len(filetext)):  
            currentline = slice_string(filetext[i])  
            currentline.remove(currentline[0])
            npycoord = np.append(npycoord, currentline, axis = 0)
        for i in range(0, len(npycoord) - 1):
            npycoord[i] = npycoord[i].replace(" ", "")
            npycoord[i] = npycoord[i].replace('\n', "")  
        atomcount=len(filetext) #2019-11-12
except:
#    print("\n>>ERROR: File not found!")
    runprogram = 0

# option variables

param_manual = 0
manual_ring_atoms = []
param_automanual = 0
automanual_ring_atoms = []
param_choosering = 0
param_choosering_now = 0
param_choosering_atoms = []
full_auto = 0
choosering_atom = []
op_zero = 0
op_one = 0
op_x = 0
op_x_value = []
op_vert = 0
op_vert_locations = []
op_displaygraph  = 0
double_bq = 0



if (len(sys.argv) < 2 and showmenu == 0):
    menu()

##UNCOMMENT THE FOLLOWING CODE FOR RUNNING PROGRAM (WONT SHOW SPECIFIC PYTHON ERRORS)
"""
try:
    if (runprogram == 1):
        inputselect()
        npycoord.shape = (atomcount,3)
        chosenatoms = []
        autosearchring()
        locate_ring()
        find_centroid()
        if(double_bq == 1):
            find_centroid()
            print("FINDING CENTROID AGAIN")
        if(param_choosering == 1):   
            param_choosering_now == 1
            locate_ring()      
        prettyprint()
        if(op_displaygraph == 1):
            displaygraph()
except:
    print("\n\n>> ERROR !")

###UNCOMMENT THE FOLLOWING CODE FOR DEBUGGING (WILL SHOW ERRORS FROM PYTHON INTERPRETER)
"""
if (runprogram == 1):
    inputselect()
    npycoord.shape = (atomcount,3)
    chosenatoms = []
    autosearchring()
    locate_ring()
    find_centroid()
    if(double_bq == 1):
        find_centroid()
    if(param_choosering == 1):   
        param_choosering_now = 1
        print("...CHOOSING RING")
        locate_ring()      
    prettyprint()
    if(op_displaygraph == 1):
        displaygraph()
###"""    


