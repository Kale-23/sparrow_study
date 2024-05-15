"""
File is called mleaflet.py ,  actually uses alternative library, folium. 

python .\mleaflet.py geodata\pcangsd_cov_10_output.cov  .\geodata\bam_list_cov_10.txt 5.0
"""

import sys
from io import StringIO 
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt 
import folium
import webbrowser 
from kmeans import run_Kmeans 


def usage():
    print("python .\mleaflet.py geodata\pcangsd_cov_10_output.cov  .\geodata\\bam_list_cov_10.txt 5.0") 

def cases_up(sqt):
    ss = list(sqt)
    stop = len(ss)-1
    for i in range(stop):
        if ss[i].islower():
            if ss[i+1].isupper():
                return(i)
    return (-1)
    
def un_squash(sqt) :
    lsqt = "{}".format(sqt)     
    cup = cases_up(lsqt)
    while( cup > 0 ) :
        sprev = "{}{}".format( lsqt[cup], lsqt[cup+1] )
        snew = "{} {}".format( lsqt[cup], lsqt[cup+1] )
        lsqt = lsqt.replace(sprev,snew)
        cup = cases_up(lsqt)
    return( lsqt )



def extract_city( basen ) :
    locus = basen.find("Mmelodia_")
    crop = basen[locus:] 
    trunc = crop.replace("Mmelodia_",'')
    minus = trunc.find("-")
    if (minus < 0 ) :
        toks = trunc.split("_") 
        citystate = toks[0]
        squash = citystate[:-2] 
    else :
        toks = trunc.split("-") 
        squash = toks[0] 
    city = squash
    while( city[-1].isupper()  ) :
        city = city[:-1] 
    
    return(  un_squash(city)  )

def PCA_mirroring( pc1 , pc2 ) :
    pc1abs = np.abs(pc1)
    pc2abs = np.abs(pc2)
    pc1sgn = np.sign(pc1)
    pc2sgn = np.sign(pc2) 
    index = np.argmax(  np.add(pc1abs , pc2abs ) )
    xflip = 1 if (pc1sgn[index] >= 0) else -1
    yflip = 1 if (pc2sgn[index] >= 0) else -1
    return xflip,yflip   

"""
2551-86994_Mmelodia_MedouieCreekNantucketMA_20160712
L1_S01_SOSP-Yellow-13_Mmelodia_ClaytorLake-VA_2011
L1_S15_2831-61091_Mmelodia_PeaIslandNWR-NC_2022
"""
geo_coords = {"Aroostook":(46.812,-68.477),
            "Assateague":(38.093,-75.207),
            "Blacksburg":(37.228,-80.413), 
            "Claytor Lake":(37.056,-80.613), 
            "Deal Island":(38.159,-75.948),
            "Delaware Seashore":(38.634,-75.067),
            "Fishing Bay":(38.260,-76.005),
            "Fortescue":(39.238,-75.172),
            "Frostburg":(39.658,-78.928), 
            "Great Bay Blvd":(39.552,-74.341),
            "Guard Shore":(37.846,-75.678),
            "Hobart Stream":(44.842,-67.251), 
            "Hog Island":(37.418352,-75.69354),
            "John HChafee":(41.439,-71.471),
            "Lieutenant Island":(41.896,-70.016),
            "Little Creek":(39.167,-75.448),
            "Loring Airfield":(46.946,-67.909),  
            "Medouie Creek Nantucket":(41.311,-70.015),
            "Pea Island":(35.718,-75.497),
            "Ponquogue Beach":(40.837,-72.492),
            "Primehook":(38.856,-75.248),
            "Rochester":(43.155,-77.609),  
            "Saxis":(37.893,-75.683), 
            "Scarborough Marsh":(43.551,-70.333),
            "Swan Bay":(39.562,-74.485),
            "Tuckahoe":(39.290,-74.754),
            "Wertheim":(40.794,-72.880) }
colors = [
    'red',
    'orange',
    'gray',
    'pink',
    'blue',
    'purple',
    'darkgreen',
    'darkred',
    'lightred',
    'beige',
    'green',
    'lightgreen',
    'darkblue',
    'lightblue',
    'darkpurple',
    'cadetblue',
    'lightgray',
    'black'
]

""" ###################################
<--- left margin 
""" #################################
if (len(sys.argv) < 4) :
    usage()
    sys.exit(-1) 

with open( sys.argv[1], 'r') as fhandle:
    lines = fhandle.readlines()
    Cfile = ''.join(lines)  
    C = np.genfromtxt( StringIO(Cfile), dtype=np.float32 ) 

print(C)

citiesOrd = [ ]
indivOrd = [ ]  
 
with open( sys.argv[2] , 'r' ) as fhandle:
    lines = fhandle.readlines() 
    for elt in lines :
      tokens = elt.strip().split('/') 
      if( len(tokens) > 3 ) : 
        basen = tokens[-1]
        basen = basen.replace("_aligned_reads_sorted.bam", "") 
        indivOrd.append ( basen ) 
        citiesOrd.append( extract_city(basen) ) 

latitudesOrd = [ ] 
longitudesOrd = [ ]
for n in range(len(citiesOrd)) :
    print( citiesOrd[n] , " , " , indivOrd[n] )
    gcoo = geo_coords[ citiesOrd[n] ]
    latitudesOrd.append( gcoo[0] )
    longitudesOrd.append( gcoo[1] )  

rowsC,colsC = C.shape  
print( "Found covariance matrix C ({} by {})".format(rowsC,colsC)  )  
print( "Found {} individuals , {} cities".format( len(indivOrd) , len(citiesOrd) )  )

groupNumD = { } 
groupDict = { } 
with open( "geosamples.tsv" , 'r' ) as fhandle :
    lines = fhandle.readlines()
    for elt in lines:
        tokens = elt.strip().split() 
        tabtoks = elt.strip().split('\t')
        if ( tokens[1].find("_Mmelodia_")  > (-1)  ) :
            groupDict[tokens[1]] = tabtoks[5]
            groupNumD[tokens[1]] = tabtoks[6]   

groupOrd = [ ]
groupNumOrd = [ ]  
for n in range( len(indivOrd) ) :
    group = groupDict[ indivOrd[n] ] 
    groupOrd.append( group )
    groupNumber = groupNumD[ indivOrd[n] ] 
    groupNumOrd.append( groupNumber ) 

print(groupOrd)
print(groupNumOrd) 


eigenvalues, eigenvectors = np.linalg.eig(C) 
explained_variance = 100.0 * ( eigenvalues / np.sum(eigenvalues) )  
order_of_importance = np.argsort(eigenvalues)[::-1] 
sorted_eigenvalues = eigenvalues[order_of_importance] 
sorted_eigenvectors = eigenvectors[:,order_of_importance] 
sorted_variances = explained_variance[order_of_importance] 

#print(sorted_eigenvectors) 

## Use  4 dimensions for K means clustering. 
uPC1 = sorted_eigenvectors[:,0]
uPC2 = sorted_eigenvectors[:,1]
uPC3 = sorted_eigenvectors[:,2]
uPC4 = sorted_eigenvectors[:,3] 

##  but only plot the first 2.
pc1mir,pc2mir = PCA_mirroring( uPC1,uPC2 ) 
PC1 = float(pc1mir) * uPC1
PC2 = float(pc2mir) * uPC2 

print( "PC1 " , PC1, sorted_variances[0] )
print( "PC2 ",  PC2, sorted_variances[1] ) 

#  Perform K means clustering with 5 assumed clusters.
#   useRand=False  attempts to have the same colors roughly in the same locations.  
rowpc = PC1.shape[0] 
Smatrix = np.hstack( ( PC1.reshape(rowpc,1) , PC2.reshape(rowpc,1) , uPC3.reshape(rowpc,1), uPC4.reshape(rowpc,1)  ) )
clusters = run_Kmeans( 5, Smatrix, useRand=False ) 
print( "clusters ", clusters ) 


# Convert cluster assignments to plot colors. 
iclusters = np.array( [int(i) for i in clusters] ) 
lclust = iclusters.tolist() 
print( lclust , len(lclust) , type(lclust) )
assigncol = [ colors[i] for i in lclust ]

# X,Y,color,group,coverage 
coverageOrd = [ (sys.argv[3]) for i in range(len(assigncol)) ]
pc1_pcntOrd = [ (sorted_variances[0]) for i in range(len(assigncol)) ]
pc2_pcntOrd = [ (sorted_variances[1]) for i in range(len(assigncol)) ]
df = pd.DataFrame(  {'X' : PC1 ,'Y':PC2 ,'color':assigncol,'group':groupOrd, 'group_num':groupNumOrd,  'coverage':coverageOrd, 'pc1_pcnt':pc1_pcntOrd, 'pc2_pcnt':pc2_pcntOrd  } ) 
print( df.head() )
outfncsv = "pca_clusters_{}.csv".format( sys.argv[3] )      
df.to_csv( outfncsv, index=False ) 

# Plot PCA first 2 principle components
# Scatter plot of PCA in 2 dimensions. 
fig = plt.figure()
plt.style.use('seaborn-v0_8-whitegrid')
ax1 = fig.add_subplot(111)
ax1.grid(False) 
ax1.scatter(PC1,PC2, s=26, c=assigncol , marker="o", edgecolors='black' )
print("ax1.scatter()") 


# Plotting attributes and options 
plt.xlim(-0.2, 0.6 )
plt.ylim(-0.4, 0.5 )
plt.title("{}X".format(sys.argv[3] ) , weight='bold' )  
plt.ylabel('PC2 ({:.2f}%)'.format( sorted_variances[1] ) , weight='bold'  )
plt.xlabel('PC1 ({:.2f}%)'.format( sorted_variances[0] ), weight='bold' )  
plt.show()

# Make a data frame with dots to show on the map
data = pd.DataFrame({
   'lon'   :longitudesOrd,
   'lat'   :latitudesOrd ,
   'name'  :indivOrd ,
   'value' :assigncol
}, dtype=str)


mapfilename = sys.argv[1].replace(".cov",".html").replace("pcangsd","foliummap")  

# Make an empty map
m = folium.Map(location=[20,0], tiles="OpenStreetMap", zoom_start=4)
feature_group = folium.FeatureGroup(name='FeatureGroup')

# add marker one by one on the map
for i in range(0,len(data)):
    fmarker = folium.Marker(
      location=[data.iloc[i]['lat'], data.iloc[i]['lon']],
      popup=data.iloc[i]['name'],
      icon=folium.Icon(  color=data.iloc[i]['value'])  )  
    fmarker.add_to( feature_group).add_to(m)

#  folium.Marker([lat, lon], popup=str(name)+': '+color+'-'+str(clname), icon=folium.Icon(color=color)).add_to(feature_group)

# Zoom and focus to our data points on the map.
southwest = data[['lat','lon']].min().values.tolist()
northeast = data[['lat','lon']].max().values.tolist()
m.fit_bounds([southwest,northeast]) 
             
# Save interactive HTML 
m.save( mapfilename )

# Open that file in the local browser. 
webbrowser.open( mapfilename ) 



