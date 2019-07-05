'''
parameters2azr.py
version 0.3

Python script that populates level, width data in parameters.out into a .azr file.

v0.1 did not transfer normalization forward.
v0.2 reads normalizations.out file, and reassigns the starting normalization values to
     only those segments that have the normalization varied during a fit.
v0.3 bugfix about counting normalization levels, added an extra break condition to exit if allnormlevels has been scanned fully.

B.Sudarsan
19 April 2018
sbalak2@lsu.edu

(Ported to python3 on July 04 2019)
'''

#Prologue: Library imports, and function declarations
import numpy as np
import lxml.etree as ET

#Filenames used:
input_azr_file = 'F17-Dec-Pratt-test4-out4.azr'##hu0junk-out.azr'  #Specify the name of the input .azr file
output_azr_file = 'F17-Dec-Pratt-test4-out5.azr' #Specify the name of the output .azr file
param_out_path_file = "./output/parameters.out" #Path to parameters.out
normalization_out_path_file = './output/normalizations.out' #Path to normalizations.out


'''
Function definitions:
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







'''
Act 1 : Get all the levels in parameters.out into a list of levels called alllevelslist_param
'''

print('Looking for parameters.out..', end=' ')

#Path to parameters.out goes here
try:
    f = open(param_out_path_file,"r+")
    lines = f.readlines()[2:] #Skip 2-line header in parameters.out and get to level data

    alllevellist_param = [] # Stores all levels in parameters.out

    n = 0 #linecount
    nmax = len(lines)
    while(True):
        line = lines[n]  #Read nth line as string
        array = line.split() #Split string at every whitespace
        #print array
        if 'J' in line:
            #If the line has J-pi, E values
            if '-' in array[2]:
                parity = -1
                J = np.float(array[2].replace('-',''))
            elif '+' in array[2]:
                parity = +1
                J = np.float(array[2].replace('+',''))
            else:
                print('Error reading parity! Check input file!')
            #Energy = np.float(array[5])
            
            Energy = read_proper_units(array[5],array[6])/1.0e6 #Convert energy to MeV no matter what units it comes in
            #Extracted j-pi, E values
            
            flag = 0
            while (len(line)>=4): #Keep reading until the file encounters a line with only '\n' that separates sections
                flag = 1
                line = lines[n]
                if 's =' in line:
                    array2 = line.split()
                    #Index:                                         5               8               11          12
                    #array2 dictionary = ['R', '=', '1', 'l', '=', '3', 's', '=', '2.0', 'G', '=', '0.000000', 'meV', 'g_int', '=', '0.000000', 'MeV^(1/2)', 'g_ext', '=', '(0.000000,0.000000)', 'MeV^(1/2)']

                    ell = np.float(array2[5])
                    ess = np.float(array2[8])
                    width = read_proper_units(array2[11],array2[12]) #Convert width to eV so we can enter it to .azr file safely
                                  
                    #Remember this order - ('J','pi','Energy','Width','s,'l')
                    alllevellist_param.append([J,parity,Energy,width,ess,ell])
                
                n = n + 1

        #Avoid infinite loop
        if flag == 0:
            n = n + 1

        if n>=nmax:
            break
finally:
    f.close()

print('done.')

print('Reading normalizations.out to find norm data..', end=' ')
allnormlist = []
try:
    f = open(normalization_out_path_file,"r+")
    lines = f.readlines()
    for line in lines:
        array = line.split()
        #Array output would look like 'Segment','Key','#n','1.043'. We store segment# 'n' and the normalization.
        allnormlist.append([array[2].replace('#',''),array[3]])
    print('done.')
finally:
    f.close()

'''
Act 2: Read through the azr file, replace the energies and widths everytime J-pi and ell, ess values match
'''
print('Converting .azr file to .xml and parsing level data..', end=' ')
#Make XML file from azure by suffixing and prefixing a category
xml_maker(input_azr_file,"temp-in.xml")

#Parse the XML file, get root tree, levels, and segmentDetails(for norm)
doc = ET.parse("temp-in.xml")
root = doc.getroot()
memoryElem = doc.find('levels') #Levels
SegmentDetails = doc.find('segmentsData') #Segments for Normalization

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

#Copy all levels first, from parameters.out --> alllevellist_param
levels =  memoryElem.text
testlevels = levels.split('\n') #Separate each level details into an array of level-details
outlevels = ''

i = 0
paramoutcounter = 0
for testlevel in testlevels:
    #Study one level at a time
    if len(testlevel)>0:
        #If level details are present
        levelarray = testlevel.split() #Split level details into an array at each space separator. This ordered array will contain all the 0-30 elements numbered in levelDict

        #From the parameters.out, the level dictionary is ('J','pi','Energy','Width','s,'l'), remember!
        #Extract j-pi, L, S values from both databases
        J_azr = np.float(levelarray[levelDict['J-channel']])
        J_pout = np.float((alllevellist_param[paramoutcounter])[0])
        Pi_azr = np.float(levelarray[levelDict['Pi-channel']])
        Pi_pout = np.float((alllevellist_param[paramoutcounter])[1])
        Ell_azr = np.float(levelarray[levelDict['2L']])/2.0
        Ell_pout = np.float((alllevellist_param[paramoutcounter])[5])
        Ess_azr = np.float(levelarray[levelDict['2S']])/2.0
        Ess_pout = np.float((alllevellist_param[paramoutcounter])[4])
        
        #When the parameters match, copy E, G values from .out file to .azr file
        if((J_azr == J_pout) and (Pi_azr==Pi_pout) and (Ell_azr == Ell_pout) and (Ess_azr == Ess_pout)):
            levelarray[levelDict['ExcEnergyChannelMeV']]=str(np.float((alllevellist_param[paramoutcounter])[2]))
            levelarray[levelDict['WidthChanneleV']]=str(np.float((alllevellist_param[paramoutcounter])[3]))
            paramoutcounter = paramoutcounter+1
        
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

print('done.')

print(allnormlist)

print('Copying segment data.. ', end=' ')
segments = SegmentDetails.text
testsegments = segments.split('\n')
outsegments = ''
counter_azr = 1 #counts the line number for the segment in the .azr file
counter_norm_out = 0
for testsegment in testsegments[1:-1]:
    segmentarray = testsegment.split()
    #print segmentarray
    datatype = np.int(segmentarray[segmentDict1['DataType']])
    if counter_norm_out <= len(allnormlist) and counter_azr == np.int((allnormlist[counter_norm_out])[0]): #Update v0.3, v0.4 to have <= len(allnormlist)
            if datatype==2:
                segmentarray[segmentDict2['Normalization']]=allnormlist[counter_norm_out][1]
            elif datatype==0 or datatype==1 or datatype==3:
                segmentarray[segmentDict1['Normalization']]=allnormlist[counter_norm_out][1]
            else:
                print('Error classifying segments.')
            counter_norm_out += 1
            
    counter_azr +=1

    outsegment = '\n'
    for string in segmentarray:
        outsegment += '   '
        outsegment += string

    #Stich things together to a large string separated by newlines
    outsegments += outsegment

outsegments+='\n'

#Locate the address of the segmentData string within the root XML file
#This would be an array of one element.
segmentloc = root.xpath("//firstElement/segmentsData")
segmentloc[0].text = outsegments #Replace the element by the new string we have synthesized from normalizations.out

print('done.')
ET.ElementTree(root).write('temp-out.xml') #Write the edited tree to a new XML file.

'''
Act 3 : Read in the new XML file, and generate the new .azr file by stripping its ends.
'''
print('Writing output file..', end=' ')
azr_maker("temp-out.xml",output_azr_file)
print('done.')

'''
Epilogue:

The levels and segments in both .azr files and .out files are sorted. This is what makes the loops doing the
replacements funny. In case a future update of azure decides to not order the levels, nested loops may be needed.

'''

