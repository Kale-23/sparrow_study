import matplotlib.pyplot as plt  
import numpy as np
import sys
from os import walk

def filter_names( fulln , target ):
    ret = []
    for f in fulln:
        if( (f.find(target)) > (-1) ) :
            ret.append(f)
    return ret

def weighted_sum( lengthD, coverageD ) :
    S = 0.00
    for chromosome in lengthD.keys() :
        chrlength = float( lengthD[chromosome] ) 
        cover = coverageD[chromosome] 
        # print( "(chrlength,cover)=" , chrlength , cover ) 
        wlen = cover * chrlength 
        S=S+wlen 
    return S   


def create_prefix( fname ) :
   return( fname.split('_')[0] ) 

def sex_acsension_numbers() :
    san = ["NC_080464.1" , "NC_080465.1" ,  "NC_080632.1" , "NC_080668.1", "NC_080669.1"]   
    return (san)  

def coverage_dictionary( prefix, fntbls , directory ) :
    sexChromos = sex_acsension_numbers() 
    D = {}
    thefile = filter_names( fntbls, prefix )
    fullpath = directory + thefile[0]  
    with open( fullpath, 'r') as fhandle:
        for line in fhandle:
            toks = line.split()
            if( (toks[0])[0] == "N" ) : 
                if(  not (toks[0] in sexChromos)  ) : 
                    chromosome = toks[0].strip()
                    depthX = float( toks[6].strip() )
                    D[chromosome] = depthX 
    return D 

def split_scaffolds( wholeD ) :
    scaffolds = {}
    chromosomes = {}
    for key in wholeD.keys() :
        if(  key.find("NC_") > (-1) ) :
            chromosomes[key] = wholeD[key] 
        if(  key.find("NW_") > (-1) ) :
            scaffolds[key] = wholeD[key]
    return chromosomes,scaffolds 

def individual_indeces( specfiles ) :
    retii = []
    for sf in specfiles :
        retii.append( sf.replace(".tbl", ''  )  )  
    return retii 
        
#####
#
#####
indirTBLS = "C:\\MCBS913\\speciestrend\\tbls\\"
species = "Acaudacuta"
taxonomy= "Ammospiza caudacuta"
maxcoverageclip = 185.0
filenamesTBLS = next(walk(indirTBLS), (None, None, []))[2]


specnamesuns = filter_names( filenamesTBLS , species )
specnames = sorted(specnamesuns)

facecolors = []
labels = []
all_data = []
all_individuals = individual_indeces( specnames )
  
for individual in all_individuals : 
    depD = coverage_dictionary( individual , specnames , indirTBLS ) 
    chromosomes,scaffolds = split_scaffolds( depD )
    chromoX = np.array( list(chromosomes.values()) , dtype=np.float32) 
    scaffX = np.array(  list(scaffolds.values()) , dtype=np.float32) 
    all_data.append( np.clip(chromoX,0,maxcoverageclip) )
    all_data.append( np.clip(scaffX,0,maxcoverageclip )  )
    labels.append(individual)
    labels.append(individual)
    facecolors.append( '#aa5522' ) 
    facecolors.append( 'white'   )  
    #
    meandepth = np.median(   np.clip(chromoX,0,maxcoverageclip)  )
    print( individual , " , " , meandepth ) 
    #


plt.figure(figsize=(8,11)) 
plt.style.use('seaborn-v0_8-whitegrid')
bplothandle = plt.boxplot( all_data , vert=False, patch_artist=True, labels=labels , showfliers=False, whis=99.0 ) 
#plt.title("{} , {} individuals".format( taxonomy,   len(specnames))  )
plt.title("") 
#plt.title("single individual , no truncated data".format( taxonomy,   len(specnames))  ) 
plt.xlabel('X depth of coverage' ) 
plt.xlim(0.0, 90.0 ) 
plt.xticks( np.arange( 0.0, 100.0, 10.0 ) )   
plt.yticks( fontsize=7  ) 
plt.subplots_adjust(left=0.35,top= 0.96, bottom=0.15 )
for patch,c in zip(bplothandle['boxes'],facecolors):
    patch.set_facecolor(c)  
plt.show()

