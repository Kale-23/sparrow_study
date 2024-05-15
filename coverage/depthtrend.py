import matplotlib.pyplot as plt  
import numpy as np
import sys
from os import walk
import pandas as pd

"""
#
"""
def reads_from_filename( filename, dframe ) :
    basen = filename.replace(".tbl", '' )
    return( dframe.loc[basen,"reads"]   )
"""
#
"""
def filter_names( fulln , target ):
    ret = []
    for f in fulln:
        if( (f.find(target)) > (-1) ) :
            ret.append(f)
    return ret

"""
#
"""
def sex_acsension_numbers() :
    san = ["NC_080464.1" , "NC_080465.1" ,  "NC_080632.1" , "NC_080668.1", "NC_080669.1"]   
    return (san)  

"""
#
"""
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

"""
#
"""
def split_scaffolds( wholeD ) :
    scaffolds = {}
    chromosomes = {}
    for key in wholeD.keys() :
        if(  key.find("NC_") > (-1) ) :
            chromosomes[key] = wholeD[key] 
        if(  key.find("NW_") > (-1) ) :
            scaffolds[key] = wholeD[key]
    return chromosomes,scaffolds 

"""
#
"""
def species_to_sdata( indir, species, taxonomy, dframe ):
    filenamesTBLS = next(walk(indir), (None, None, []))[2]
    specnamesuns = filter_names( filenamesTBLS , species )
    specnames = sorted(specnamesuns)

    reads_D = { }
    NWavg_D = { } 
    NCavg_D = { } 

    for individual in specnames : 
        ilreads = reads_from_filename(individual, dframe ) 
        print( "individual ", individual )
        print( "    reads " , ilreads ) 
        depD = coverage_dictionary( individual , specnames , indirTBLS ) 
        chromosomes,scaffolds = split_scaffolds( depD )
        chromoX = np.array( list(chromosomes.values()) , dtype=np.float32) 
        scaffX = np.array(  list(scaffolds.values()) , dtype=np.float32)
        reads_D[individual] = float(ilreads) / float(100_000_000) 
        NWavg_D[individual] = np.median( scaffX )
        NCavg_D[individual] = np.median( chromoX ) 

    xNCl = [ ]
    yNCl = [ ]
    xNWl = [ ]
    yNWl = [ ] 

    for individual in specnames :  
        xNCl.append( reads_D[individual] )
        yNCl.append( NCavg_D[individual] ) 
        xNWl.append( reads_D[individual] )
        yNWl.append( NWavg_D[individual] )  

    print( "xNCl  ", xNCl , " ", len(xNCl) )
    print( "yNCl  ", yNCl , " ", len(yNCl) )
    print( "xNWl  ", xNWl , " ", len(xNWl) )
    print( "yNWl  ", yNWl , " ", len(yNWl) )
    indivs = len(xNCl) 

    return xNCl,yNCl,xNWl,yNWl,indivs  
    
"""
# Return only a slope and a R^2. 
# This assumes the constraint (0,0) must be in the estimated line.
"""
def linear_regression( xlist , ylist ):
    xreg = np.array(  xlist , dtype=np.float32 ) 
    yreg = np.array(  ylist , dtype=np.float32 )
    xreg = xreg[:,np.newaxis]
    slope,resid,_,_ = np.linalg.lstsq(xreg,yreg) 
    R2 = 1 - resid / (yreg.size * yreg.var())
    return  slope,R2

"""
# Return a Dataframe for the requested X coverage depths , X_reqs
# Assume a model y = Ax + 0.0 
# multf is applied to reach integer reads, e.g. "in units of 100 million reads"    
"""
def read_estimates( A , X_reqs, multf ):
    df = pd.DataFrame(columns = ['Depth req', 'Reads'] ) 
    depths = [ ]
    for elt in X_reqs :
        depths.append( "{}X".format(elt)  ) 

    ysu = np.array(  X_reqs , dtype=np.float64 ) 
    ysu = ysu * (1.0/A)
    ysu = ysu * multf 
    Reads = list(  ysu.astype(np.int64) ) 
    for elt in (zip(depths , Reads)):
        df.loc[len(df.index)] = [ elt[0] ,elt[1]  ] 

    return  ( df )

# #############
#<- left margin
# #############
df = pd.read_csv( "flagsarchive.csv" , index_col="basename" ) 
print( df.head() ) 

indirTBLS = "C:\\MCBS913\\speciestrend\\tbls\\"
batch = [ "Acaudacuta" , "Anelsoni" , "Mgeorgiana" ]
taxonomy= "Melospiza georgiana"

xNC = [ ]
yNC = [ ] 
xNW = [ ]
yNW = [ ]
gtotal = 0 
for species in batch:
    xnc,ync,xnw,ynw,total = species_to_sdata( indirTBLS , species, taxonomy , df )  
    xNC += xnc
    yNC += ync
    xNW += xnw
    yNW += ynw 
    gtotal =gtotal+total 

Anw,R2nw = linear_regression(xNW,yNW) 
Anc,R2nc = linear_regression(xNC,yNC) 
print( "Anw=", Anw, "  R2nw=",R2nw )
print( "Anc=", Anc, "  R2nc=",R2nc )
depthRequestTable = read_estimates(Anc,[1,2,5,10], 1.00e8 )
print(depthRequestTable.to_string(index=False))

fig = plt.figure()
plt.style.use('seaborn-v0_8-whitegrid')
ax1 = fig.add_subplot(111)

# Plot points
ax1.scatter(xNC,yNC, s=10, c='b', marker="s", label='NC_')
ax1.scatter(xNW,yNW, s=10, c='#bbbb00', marker="o", label='NW_')

# Plot regression lines.  
pts=np.linspace(0,2.6,26) 
ax1.plot( pts, Anc * pts, 'b-' )
ax1.plot( pts, Anw * pts, color='#bbbb00', linestyle='solid') 

# Plot anchors between NC_ and NW_ pairs. 
anchors = np.array( (xNC,yNC,yNW) ).T 
for row in anchors :
    top = max( [row[1] , row[2]] )
    bott = min(  [row[1] , row[2]] )
    plt.vlines(x=row[0], ymin=bott, ymax=top, color="k", linestyle='dashed',linewidth=0.5)


plt.xlim(0.0, 2.5 )
plt.ylim(0.0, 70.0 )
leg = plt.legend()
leg.get_frame().set_edgecolor('b')
plt.legend(loc='lower right', frameon=True)
species=""
plt.title("{} , {} individuals".format(species,gtotal) )  
plt.ylabel('median X depth' )
plt.xlabel("Read Count (100 millions)")  
plt.show()