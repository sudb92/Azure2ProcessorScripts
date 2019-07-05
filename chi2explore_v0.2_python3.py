'''
chi2explore.py
v0.2

Based on
parameters2azr.py
version 0.4

Python script to read an .azr file, and do a chi2 map around the current values using calculations, as follows:

1. Print a ID'd list of all parameters free to vary in a fit. These are the parameter-space used.
2. Ask the user for how many parameters are to be varied in the map - 1 or 2.
3. Upon getting the response, accept the IDs of parameters to vary.
4. Prompt using the ID and the present value, for the user to enter start value, end value and step size for all the varied params.
5. Use pexpect to run Azure in text mode and perform calculations at each value. Write progres to a log file, and results to a .dat file

update v0.2
6. Takes care of bookkeeping when the energy and width of the same state are varied.
When this is the case, the updated program changes energies for all the involved sublevels and changes widths appropriately.

B.Sudarsan
6 Feb 2019
sbalak2@lsu.edu

Ported to python3 on 4 July 2019
'''

#Prologue: Library imports, and function declarations
import numpy as np
import lxml.etree as ET
import os
import pexpect as px
import sys
import time
import matplotlib.pyplot as plt


#Filenames used:
input_azr_file = 'F17-Dec-Pratt-test4-out4.azr'##hu0junk-out.azr'  #Specify the name of 18the input .azr file
working_azr_file = input_azr_file[:-4]+'-chi2test.azr' #Specify the name of the output .azr file
param_out_path_file = "./output/parameters.out" #Path to parameters.out
normalization_out_path_file = './output/normalizations.out' #Path to normalizations.out
chi2_out_path_file = './output/chiSquared.out' #Path to chiSquared.out

AZURE_EXECUTABLE_FULL_PATH = "/home/sud/Desktop/NuclearPhysics/Azure/AZURE2/build/AZURE2"


#Boolean switches to set
save_copy_of_azr_files = True
save_chiSquared_out_files = False

'''
Function and Dictionary definitions:
'''
def read_proper_units(valuestr, unitstr):
    '''
    A dictionary to convert energies in multiple units to eV
    '''
    units_dict = {'meV':1e-3,'eV':1,'keV':1e3,'MeV':1e6,'GeV':1e9}
    for key in list(units_dict.keys()):
        if key in unitstr:
            multiplier = np.float(units_dict[key])
            return np.float(valuestr)*multiplier #value in eV!

def xml_maker(infile,outfile):
    '''
    xml_maker(string infile, string outfile):

    Function to convert .azr files to proper readable .xml files by prefixing and suffixing the appropriate XML tag
    '''
    f = open(infile,"r")
    fo = open(outfile,"w")

    lines = f.readlines()
    head1 = "<firstElement>\n"
    foot1 = "</firstElement>"

    fo.write(head1)

    for line in lines:
        fo.write(line)

    fo.write(foot1)

    f.close()
    fo.close()

def azr_maker(infile,outfile):
    '''
    azr_maker(string infile, string outfile):

    Function to convert .xml files to proper readable .azr files by getting rid of prefix/suffix tag in XML file
    '''
    f = open(infile,"r")
    fo = open(outfile,"w")

    lines = f.readlines()
    for line in lines[1:-1]:
        fo.write(line)

    f.close()
    fo.close()


#Dictionary mapping the contents of one line of 'levels' by category
levelDict = {'J-channel':0,
             'Pi-channel':1,
            'ExcEnergyChannelMeV':2,
             'FixE?':3,
             'UnknownFlag':4,
             'ParticlePair#':5,
             '2S':6,
             '2L':7,
             'ChannelIndex':8,
             'IncludeLevel?':9,
             'FixWidth?':10,
             'WidthChanneleV':11,
             'J-light':12,
             'Pi-light':13,
             'J-heavy':14,
             'Pi-heavy':15,
             'ExcEnergyInputMeV':16,
             'A-Light':17,
             'A-Heavy':18,
             'Z-light':19,
             'Z-Heavy':20,
             'UnknownSeparationEnergyMeV':21,
             'ParticlePair#SeparationEnergyMeV':22,
             'Unknown#1':23,
             'Unknown#1':24,
             'Unknown#1':25,
             'Unknown#1':26,
             'ChannelRadiusfm':27,
             'Unknown#1':28,
             'Unknown#1':29,
             'Unknown#1':30
             }

#Two dictionaries to use with segment data
#Dictionary to use when segments have angle-integrated,differential,or angle-integrated-total-capture data
segmentDict1 = {'Include?':0,
                'EntrancePair':1,
                'ExitPair':2,#Becomes -1 when DataType=3. 
                'LowLabEnergyMeV':3,
                'HighLabEnergyMeV':4,
                'LowLabAngleDeg':5,#Automatically goes to 0 for DataType=0,3
                'HighLabAngleDeg':6,#Automatically goes to 180 for DataType=0,3
                'DataType':7,#0 - Angle Integrated, 1 - Differential, 3- Angle Integrated Total Capture
                'Normalization':8,
                'VaryNorm?':9,
                'NormError%':10,
                'DataFilePath':11
                }

#Dictionary to use when segments have phase-shift data, 
segmentDict2 = {'Include?':0,
                'EntrancePair':1,
                'ExitPair':2,#Becomes -1 when DataType=3. 
                'LowLabEnergyMeV':3,
                'HighLabEnergyMeV':4,
                'LowLabAngleDeg':5,#Automatically goes to 0 for DataType=0,3
                'HighLabAngleDeg':6,#Automatically goes to 180 for DataType=0,3
                'DataType':7,#2 - Phase-shift
                'J':8,
                'l':9,
                'Normalization':10,
                'VaryNorm?':11,
                'NormError%':12,
                'DataFilePath':13
                }


'''
Act 1: Read through the input azr file, find the energy, non-zero width parameters that have been allowed to vary according to the 'tick' marks. 
'''
print('Converting input .azr file to .xml and parsing level data..', end=' ')
#Make XML file from azure by suffixing and prefixing a category
xml_maker(input_azr_file,"temp-in.xml")

#Parse the XML file, get root tree, levels, and segmentDetails(for norm)
doc = ET.parse("temp-in.xml")
root = doc.getroot()
memoryElem = doc.find('levels') #Levels
SegmentDetails = doc.find('segmentsData') #Segments for Normalization

#Copy all levels first, from parameters.out --> alllevellist_param
levels =  memoryElem.text
testlevels = levels.split('\n') #Separate each level details into an array of level-details
outlevels = ''
print('')
Evarycount = 0
Evarylist = []

Widthvarycount = 0
Widthvarylist = []

NumChannelsAtEvalue = 1
Engy_prev = np.nan
J_prev = np.nan
Pi_prev = np.nan
FixE_prev = np.nan

i = 0

print('List of all found parameters:')
print("FixE?\tE(MeV)\tFixW?\tW(eV)\tJ\tPi\tL\tS")



paramoutcounter = 0
for testlevel in testlevels:
    #Study one level at a time
    if len(testlevel)>0:
        #If level details are present
        levelarray = testlevel.split() #Split level details into an array at each space separator. This ordered array will contain all the 0-30 elements numbered in levelDict

        #From the parameters.out, the level dictionary is ('J','pi','Energy','Width','s,'l'), remember!
        #Extract j-pi, L, S values from both databases
        J_azr = np.float(levelarray[levelDict['J-channel']])
        #J_pout = np.float((alllevellist_param[paramoutcounter])[0])
        Pi_azr = np.float(levelarray[levelDict['Pi-channel']])
        #Pi_pout = np.float((alllevellist_param[paramoutcounter])[1])
        Ell_azr = np.float(levelarray[levelDict['2L']])/2.0
        #Ell_pout = np.float((alllevellist_param[paramoutcounter])[5])
        Ess_azr = np.float(levelarray[levelDict['2S']])/2.0
        #Ess_pout = np.float((alllevellist_param[paramoutcounter])[4])

        Includelevel = np.int(levelarray[levelDict['IncludeLevel?']])
        FixE = np.int(levelarray[levelDict['FixE?']])
        FixW = np.int(levelarray[levelDict['FixWidth?']])
        Engy = np.float(levelarray[levelDict['ExcEnergyChannelMeV']])
        Widthu = np.float(levelarray[levelDict['WidthChanneleV']])

        #print (Engy_prev,Engy),
        if Engy_prev == Engy and J_prev == J_azr and Pi_prev == Pi_azr and FixE_prev == FixE:
            NumChannelsAtEvalue = NumChannelsAtEvalue + 1
        elif Engy_prev != np.nan:
            NumE = NumChannelsAtEvalue + 1
            #print 'Group',NumE,' levels above'
            NumChannelsAtEvalue = 0
            if(FixE_prev==0) and ((Engy_prev,J_prev,Pi_prev,NumE) not in Evarylist) and (Includelevel_prev==1):
                Evarycount = Evarycount + 1
                Evarylist.append((Engy_prev,J_prev,Pi_prev,NumE))
                print('Level: ',(Engy_prev,J_prev,Pi_prev,NumE))

        Includelevel_prev = Includelevel
        FixE_prev = FixE
        Engy_prev = Engy
        J_prev = J_azr
        Pi_prev = Pi_azr
            
        #Zero width states are not going to be touched by azure whether or not they're varied. Still include them in the list for completeness
        if (FixW==0) and ((Engy,J_azr,Pi_azr,Ell_azr,Ess_azr,Widthu) not in Widthvarylist) and (Includelevel==1):
            Widthvarycount = Widthvarycount + 1
            Widthvarylist.append((Engy,J_azr,Pi_azr,Ell_azr,Ess_azr,Widthu))

        print(FixE, '\t', Engy, '\t', FixW, '\t', Widthu, '\t', J_azr,'\t', Pi_azr,'\t', Ell_azr,'\t', Ess_azr)
        
    else:
        continue

print(Evarycount, ' energies varied.')
print(Widthvarycount, ' widths are varied.')

ctr = 1
print('Energies varied\n(ID,(Energy, J, pi)):')
for E1 in Evarylist:
    #print ctr,'\t',E1
    Evarylist[ctr-1] = (ctr,Evarylist[ctr-1])
    print(Evarylist[ctr-1])
    ctr = ctr + 1

print('Widths varied\n(ID,(Energy,J,Pi,L,S,Width)):')
for W1 in Widthvarylist:
    #print ctr,'\t',W1
    Widthvarylist[ctr-1-Evarycount] = (ctr,Widthvarylist[ctr-1-Evarycount])
    print(Widthvarylist[ctr-1-Evarycount])
    ctr = ctr + 1

#print Evarylist
#print Widthvarylist

#stringtomyself1 = 'Just identify the ID of the parameter to vary. Leave it to the user what he varies. It could be energy/width, it could be width/width. Just accept two IDs'
#stringtomyself2 = ' systematically. For energies, make a list of matching segments for each index, all of which will be replaced during each iteration. Widths its just one segment'
#print stringtomyself1+stringtomyself2


idstovary = []

#while 1:
#print 'How many parameters to vary? (1 or 2):'
#numparam = raw_input('How many parameters to vary? (1 or 2):')
#numparam = int(numparam)
#if numparam not in [1,2]:
#    numparam = 1

numparam = 2
print('I can vary upto two parameters at a time..')
ctr = 0
while ctr < int(numparam):
    id1 = input('Enter ID of param number '+str(ctr+1)+' :')
    idstovary.append(int(id1))
    ctr = ctr + 1

print(idstovary)
thingstovary = []

for ID in idstovary:
    if(ID-1 < len(Evarylist)):
        a = Evarylist[ID-1][0]
        if a == ID:
            thingstovary.append(Evarylist[ID-1])
            print('Energy found:',Evarylist[ID-1])

    print(ID)
    print(len(Evarylist))
    print(ID-Evarycount)
    print(len(Widthvarylist))
    if abs(ID-1-Evarycount)<len(Widthvarylist):
        b = Widthvarylist[ID-Evarycount-1][0]
        if b== ID:
            thingstovary.append(Widthvarylist[ID-Evarycount-1])
            print('Width found:',Widthvarylist[ID-Evarycount-1])

if len(thingstovary)==0:
    print('Enter the right indices and try again, exiting..')
    exit()

thingstovary_withrange = []

for thing in thingstovary:
    if len(thing[1]) == 4:
        lowE = input('Varying energy at '+str(thing[1])+', enter low value:')
        lowE = float(lowE)
        highE = input('Varying energy at '+str(thing[1])+', enter high value:')
        highE = float(highE)
        NstepsE = input('Varying energy at '+str(thing[1])+', enter # of steps:')
        NstepsE = int(NstepsE)
        if(lowE>highE) or (NstepsE<=0):
            print('Erroneous range.. choosing default values..', end=' ')
            lowE = thing[1][1] - 0.1*thing[1][1]
            highE = thing[1][1] + 0.1*thing[1][1]
            NstepsE = 10
            print(' lowE:',lowE,' highE:',highE,' NstepsE:',NstepsE)
        
        thing = (thing[0],thing[1],(lowE,highE,NstepsE),"Energy")
        thingstovary_withrange.append(thing)
        
    elif len(thing[1]) == 6:
        lowW = input('Varying width at '+str(thing[1])+', enter low value:')
        lowW = float(lowW)
        highW = float(input('Varying width at '+str(thing[1])+', enter high value:'))
        highW = float(highW)
        NstepsW = int(input('Varying width at '+str(thing[1])+', enter # of steps:'))
        NstepsW = int(NstepsW)
        if(lowW>highW) or (NstepsW<=0):
            print('Erroneous range.. choosing default values..', end=' ')
            lowW = thing[1][5] - 0.1*thing[1][5]
            highW = thing[1][5] + 0.1*thing[1][5]
            NstepsW = 10
            print(' lowW:',lowW,' highW:',highW,' NstepsW:',NstepsW)
        
        
        thing = (thing[0],thing[1],(lowW,highW,NstepsW),"Width")
        thingstovary_withrange.append(thing)

print(thingstovary_withrange[0],thingstovary_withrange[1])

'''
#Find out if we're varying energy and width of the same level
HAVE_COMMON_ENERGIES = False
if thingstovary_withrange[0][1][0] == thingstovary_withrange[1][1][0] :
    HAVE_COMMON_ENERGIES = True
'''
#exit()
#print thingstovary
print(thingstovary_withrange)

del thingstovary

print('About to vary '+str(len(thingstovary_withrange))+' parameters to study chi2 dependence..')
for thing in thingstovary_withrange:
    if len(thing[1])==4:
        print('Energy ',thing[1][0],' MeV will be varied in ',thing[2][2],' steps from ',thing[2][0],' to ',thing[2][1],' ..')
    elif len(thing[1])==6:
        print('Width ',thing[1][5],' keV will be varied in ',thing[2][2],' steps from ',thing[2][0],' to ',thing[2][1],' ..')
#exit()

'''
Act 2 : Generate the working .azr file by copying input file.
'''
print('Preparing working file at ',working_azr_file,' ...')
os.system("cp "+input_azr_file+" "+working_azr_file)
print('done.')

os.system("mkdir chi2search_folder")

'''
Act 3 : Take the workig .azr file and run calculations in a while loop
'''
#Hold the tuple and a description like 'Width' or 'Energy'
working_param1 = list(thingstovary_withrange[0][1])
working_param2 = list(thingstovary_withrange[1][1])

energies_are_equal = False

if(working_param2[0] == working_param1[0]):
    energies_are_equal = True
    print('Same E for both parameters:', end=' ')
else:
    print('Not same E for both parameters:', end=' ')

energy_width_combo=False
energy_energy_combo = False
width_width_combo = False
if(thingstovary_withrange[0][3]!=thingstovary_withrange[1][3]):
    print('energy and width.')
    energy_width_combo = True
elif (thingstovary_withrange[0][3]=='Width'):
    print('width and width.')
    width_width_combo = True
elif (thingstovary_withrange[0][3]=='Energy'):
    print('energy and energy.')
    energy_energy_combo = True

               

print(working_param1)
print(working_param2)

#exit()

param1array = np.linspace(thingstovary_withrange[0][2][0],thingstovary_withrange[0][2][1],thingstovary_withrange[0][2][2])
param2array = np.linspace(thingstovary_withrange[1][2][0],thingstovary_withrange[1][2][1],thingstovary_withrange[1][2][2])

print(param1array, ' ', end=' ')
print(param2array)

input("press key to continue..")

#Output list holding (p1,p2,chisquared)
chisqlist = []
index = 0
for p1 in param1array:
    for p2 in param2array:


        #Make XML file from working azure file by suffixing and prefixing a category
        xml_maker(working_azr_file,"temp-in.xml")

        #Parse the XML file, get root tree, levels, and segmentDetails(for norm)
        doc = ET.parse("temp-in.xml")
        root = doc.getroot()
        memoryElem = doc.find('levels') #Levels
        SegmentDetails = doc.find('segmentsData') #Segments for Normalization
        #Copy all levels first, from parameters.out --> alllevellist_param
        levels =  memoryElem.text
        testlevels = levels.split('\n') #Separate each level details into an array of level-details
        outlevels = ''

        print('Parameters in present iteration:',p1,' ',p2)
        Ectr1 = 0
        Ectr2 = 0
        #print working_param1,' ',working_param2
        
        i = 0
        paramoutcounter = 0
        for testlevel in testlevels:
            #Study one level at a time
            if len(testlevel)>0:
                #If level details are present
                levelarray = testlevel.split() #Split level details into an array at each space separator. This ordered array will contain all the 0-30 elements numbered in levelDict

                #From the parameters.out, the level dictionary is ('J','pi','Energy','Width','s,'l'), remember!
                #Extract E,W, j-pi, L, S values from both databases
                E_azr = np.float(levelarray[levelDict['ExcEnergyChannelMeV']])
                W_azr = np.float(levelarray[levelDict['WidthChanneleV']])
                
                J_azr = np.float(levelarray[levelDict['J-channel']])
                #J_pout = np.float((alllevellist_param[paramoutcounter])[0])
                Pi_azr = np.float(levelarray[levelDict['Pi-channel']])
                #Pi_pout = np.float((alllevellist_param[paramoutcounter])[1])
                Ell_azr = np.float(levelarray[levelDict['2L']])/2.0
                #Ell_pout = np.float((alllevellist_param[paramoutcounter])[5])
                Ess_azr = np.float(levelarray[levelDict['2S']])/2.0
                #Ess_pout = np.float((alllevellist_param[paramoutcounter])[4])

                FixE = np.int(levelarray[levelDict['FixE?']])
                FixW = np.int(levelarray[levelDict['FixWidth?']])
                    
             
                if(len(working_param2)==4):
                    if (E_azr == working_param2[0]) and (J_azr == working_param2[1]) and (Pi_azr == working_param2[2]):#(FixE==0) and
                        levelarray[levelDict['ExcEnergyChannelMeV']]=str(p2)
                        Ectr2 = Ectr2+1 #Count the number of sublevels edited

                        #print 'Ctr:',Ectr2
                        #print 'Lim:',working_param2[3]

                        #Replace energy only the right number of times, count just the right number of sublevels needed.
                        if Ectr2 == working_param2[3]: #working_paramx[3] holds the number of sublevels at the same energy
                            working_param2[0] = p2
                            Ectr2 = 0.0                        
                
                
                if(len(working_param1)==4):
                    if (E_azr == working_param1[0]) and (J_azr == working_param1[1]) and (Pi_azr == working_param1[2]): #(FixE==0) and 
                        levelarray[levelDict['ExcEnergyChannelMeV']]=str(p1)
                        Ectr1 = Ectr1+1 #Count the number of sublevels edited
                        #print 'Ctr:',Ectr1
                        #print 'Lim:',working_param1[3]


                        #Replace energy only the right number of times, count just the right number of sublevels needed. 
                        if Ectr1 == working_param1[3]: #working_paramx[3] holds the number of sublevels at the same energy
                            working_param1[0] = p1
                            Ectr1 = 0.0

                        
                if(len(working_param1)==6):
                    #'''
                    if(FixW==0) and (E_azr == working_param1[0]) and (J_azr == working_param1[1]) and (Pi_azr == working_param1[2]) and (Ell_azr == working_param1[3]) and (Ess_azr == working_param1[4]) and (W_azr == working_param1[5]):
                        levelarray[levelDict['WidthChanneleV']]=str(p1)
                        working_param1[5] = p1
                    
                        if energies_are_equal and energy_width_combo: #Then working_param2 would have the energy we need
                            levelarray[levelDict['ExcEnergyChannelMeV']]=str(p2)
                            working_param1[0] = p2
                    '''#Update version 0.2 below
                    if(E_azr == working_param1[0]) and (J_azr == working_param1[1]) and (Pi_azr == working_param1[2]):
                        if energies_are_equal: #Then working_param2 would have the energy we need
                            levelarray[levelDict['ExcEnergyChannelMeV']]=str(p2)
                            working_param1[0] = p2
                        if (FixW==0) and (Ell_azr == working_param1[3]) and (Ess_azr == working_param1[4]) and (W_azr == working_param1[5]):
                            levelarray[levelDict['WidthChanneleV']]=str(p1)
                            working_param1[5] = p1
                    #'''    
                    
                if(len(working_param2)==6):
                    #'''
                    if(FixW==0) and (E_azr == working_param2[0]) and (J_azr == working_param2[1]) and (Pi_azr == working_param2[2]) and (Ell_azr == working_param2[3]) and (Ess_azr == working_param2[4]) and (W_azr == working_param2[5]):
                        levelarray[levelDict['WidthChanneleV']]=str(p2)
                        working_param2[5] = p2

                        if energies_are_equal and energy_width_combo: #Then working_param1 would have the energy we need
                            levelarray[levelDict['ExcEnergyChannelMeV']]=str(p1)
                            working_param2[0] = p1
                    '''#Update version 0.2 below
                    if (E_azr == working_param2[0]) and (J_azr == working_param2[1]) and (Pi_azr == working_param2[2]):
                        if energies_are_equal: #Then working_param2 would have the energy we need
                            levelarray[levelDict['ExcEnergyChannelMeV']]=str(p1)
                            working_param2[0] = p1
                            #print p2, working_param2[5], W_azr
                        if(FixW==0) and (Ell_azr == working_param2[3]) and (Ess_azr == working_param2[4]) and (W_azr == working_param2[5]):
                            levelarray[levelDict['WidthChanneleV']]=str(p2)
                            working_param2[5] = p2
                    '''
                        
                #print working_param1,' ',working_param2
                #For energies, match E, J and pi.
                #For widths, match E, J, pi, L, S, W
                
                '''
                #When the parameters match, copy E, G values from .out file to .azr file
                if((J_azr == J_pout) and (Pi_azr==Pi_pout) and (Ell_azr == Ell_pout) and (Ess_azr == Ess_pout)):
                    levelarray[levelDict['ExcEnergyChannelMeV']]=str(np.float((alllevellist_param[paramoutcounter])[2]))
                    levelarray[levelDict['WidthChanneleV']]=str(np.float((alllevellist_param[paramoutcounter])[3]))
                    paramoutcounter = paramoutcounter+1
                '''
                #Convert levelarray back to string
                outlevel = '\n'
                for string in levelarray:
                    outlevel += '   '
                    outlevel += string

                #Stich things together to a large string separated by newlines
                outlevels += outlevel
                
            else:
                continue


        outlevels += '\n'

        #Locate the address of the levels string within the root XML file
        #This would be an array of one element.
        levelloc = root.xpath("//firstElement/levels")
        levelloc[0].text = outlevels #Replace the element by the new string we've synthesized


        ET.ElementTree(root).write('temp-out.xml') #Write the edited tree to a new XML file.

        '''
        Act 3 : Read in the new XML file, and generate the updated working .azr file by stripping its ends.
        '''
        print('Writing output azr file..', end=' ')
        azr_maker("temp-out.xml",working_azr_file)
        print('done.')

        
        #Act 4: Run Azure with updated parameters, generate chiSquared.out
        child = px.spawn(str(AZURE_EXECUTABLE_FULL_PATH)+" "+working_azr_file+" --no-gui --use-brune",timeout=None)
        #child.logfile = sys.stdout
        child.expect(".*azure2:")
        child.sendline("1")
        child.expect(".*new file")
        child.sendline("")
        child.expect("Thanks for using AZURE2.")
        #child.expect(px.EOF)
        time.sleep(1.5)
        child.close()
        #if child.isalive(): 
            


        #Read-in chiSquared.out and extract the chiSquared value
        print('Reading chiSquared.out to find chi2 data..')
        
        try:
            f = open(chi2_out_path_file,"r+")
            lines = f.readlines()
            #print lines
            for line in lines[1:-1]:
                array = line.split()
                #print array
                #Array output would look like
                #['Segment', '#2', 'Chi-Squared/N:', '76.383']
                #['Total', 'Chi-Squared:', '6339.79']
                #We need the second element in every second line
                chisqlist.append((index,p1,p2,float(array[2])))
                print("p1:",p1,"  p2:",p2,"  chi2:",float(array[2]))
        finally:
            f.close()
        print('done.\n')
        
        if save_chiSquared_out_files:
            os.system("cp ./output/chiSquared.out ./chi2search_folder/chiSquared-"+str(index)+".out")
        if save_copy_of_azr_files:
            os.system("cp "+working_azr_file+" ./chi2search_folder/"+working_azr_file[:-4]+"-"+str(index)+".azr")
        index = index + 1
        #exit()

np.savetxt("chisquared-output.dat",chisqlist,fmt="%1.4f")

'''
Epilogue:

The levels and segments in both .azr files and .out files are sorted. This is what makes the loops doing the
replacements funny. In case a future update of azure decides to not order the levels, nested loops may be needed.

'''

