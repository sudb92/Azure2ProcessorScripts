'''
pretty_printer_python3.py
(formerly spreadsheet_parser.py)
version 0.3

Python program that parses a parameters.out file to provide a readable output file

B.Sudarsan
24 June 2018
sbalak2@lsu.edu

Ported to python3 on 4 July 2019
'''

#Prologue: Library imports, and function declarations
import numpy as np
import lxml.etree as ET

#Provide directories
param_out_path_file = "parameters.out" #Path to parameters.out
output_path_file = 'parsed-level-width.txt' #Output filename
normalization_out_path_file = 'normalizations.out' #Path to normalizations.out
chi2_out_path_file = 'chiSquared.out' #path to chiSquared.out

WRITE_DUMMY_LEVELS = False

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
            print(str(multiplier)+'\t'+valuestr+'\t'+unitstr+'\n')
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
f = open(param_out_path_file,"r+")

#Path to parameters.out goes here
try:
    f = open(param_out_path_file,"r+")
    lines = f.readlines()[2:] #Skip 2 line header

    alllevellist_param = []

    n = 0 #linecount
    nmax = len(lines)
    while(True):
        line = lines[n]
        array = line.split()
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


print('Reading chiSquared.out to find chi2 data..', end=' ')
allchi2list = []
try:
    f = open(chi2_out_path_file,"r+")
    lines = f.readlines()
    x = 1 #Variable to keep track of lines in pairs of 2
    ctr = 1 #counter for no of segments
    for line in lines:
        if len(line) > 1:
            if x == 2 :
                array = line.split()
                #print array,array[2]
                totchi2 = np.float(array[2])
            if x == 1: 
                array = line.split()
                print(array, array[3])
                chi2overN = np.float(array[3])
            if x == 2 :
                x = 1
                allchi2list.append([ctr, chi2overN,totchi2])
                ctr = ctr + 1
            else:
                x = 2
            
                
        #Array output would look like 'Segment','Key','#n','1.043'. We store segment# 'n' and the normalization.
        #allchi2list.append([array[2].replace('#',''),array[3]])
        
    print('done.')
finally:
    f.close()


f = open("parsed-level-width.txt","w+")
f2 = open("parsed-level-width2.txt","w+")
f3 = open("parsed-level-width3.txt","w+")

f.write("Energy(MeV)\tJ\tpi\tL\tS\tWidth(keV)\n")
f2.write("Energy(MeV)\tJ\tpi\tL\tS\tWidth(keV)\n")
f3.write("Energy(MeV)\tJ\tpi\tL\tS\tWidth(keV)\n")

for level in alllevellist_param:
    #Remember this order - ('J','pi','Energy','Width','s,'l')
    
    f.write(str(level[2])+'\t'+str(int(level[0]))+'\t'+str(level[1])+'\t'+str(level[5])+'\t'+str(level[4])+'\t'+str(level[3]/1.0e3)+'\n')
    if(int(level[2]) == 20) and int(level[3])==0 :
        print('skip\n')
    else :
        f2.write(str(level[2])+'\t'+str(int(level[0]))+'\t'+str(level[1])+'\t'+str(level[5])+'\t'+str(level[4])+'\t'+str(level[3]/1.0e3)+'\n')
        if np.abs(level[3])>0:
            f3.write(str(level[2])+'\t'+str(int(level[0]))+'\t'+str(level[1])+'\t'+str(level[5])+'\t'+str(level[4])+'\t'+str(level[3]/1.0e3)+'\n')

f.write("Norm:\n")
f2.write("Norm:\n")
f3.write("Norm:\n")
for norm in allnormlist:
    f.write(str(norm[0])+'\t'+str(norm[1])+'\n')
    f2.write(str(norm[0])+'\t'+str(norm[1])+'\n')
    f3.write(str(norm[0])+'\t'+str(norm[1])+'\n')


f.write("Chi2:\n")
f2.write("Chi2:\n")
f3.write("Chi2:\n")
for chi2 in allchi2list:
    f.write(str(chi2[0])+'\t'+str(chi2[1])+'\t'+str(chi2[2])+'\n')
    f2.write(str(chi2[0])+'\t'+str(chi2[1])+'\t'+str(chi2[2])+'\n')
    f3.write(str(chi2[0])+'\t'+str(chi2[1])+'\t'+str(chi2[2])+'\n')

f.close()
f2.close()
f3.close()

