import sys

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


def joseph_coords( Dspacey ) :
    ret = { }
    for k in (Dspacey.keys()) :
        rk = k.replace(' ','')
        ret[rk] = (Dspacey[k])
    return (ret)  
        

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

def extract_state( basen ) : 
    locus = basen.find("Mmelodia_")
    crop = basen[locus:] 
    trunc = crop.replace("Mmelodia_",'')
    minus = trunc.find("-")
    if (minus < 0 ) :
        toks = trunc.split("_") 
        citystate = toks[0]
        state = citystate[-2:] 
    else :
        state = "{}{}".format( trunc[minus+1] , trunc[minus+2] ) 
    return (state) 

def extract_year(basen) :
    locus = basen.find("Mmelodia_")
    crop = basen[locus:] 
    trunc = crop.replace("Mmelodia_",'')
    toks = trunc.split("_") 
    fulldate = toks[-1]
    return( fulldate[0:4] ) 


"""  ##########################
<-- left margin
"""  

joD = joseph_coords( geo_coords )
print( joD) 

fullpath = "samples.tsv"
outlines = [ ] 
with open( fullpath, 'r') as fhandle:
    for line in fhandle:
        basen = line.strip() 
        if(  (basen.find('subspecies') ) > (-1) ) :
            linemod = "{}\tlat\tlong\n".format(basen)
        else :
            match = ''
            for k in (joD.keys()) : 
                if(  (basen.find(k) ) > (-1) ) :
                    match=k
            
            assert len(match)>0, "Cannot find that subspecies. Quiting."
            geoc = joD[match] 
            Dlat = geoc[0]
            Dlong = geoc[1] 
            linemod = "{}\t{}\t{}\n".format(basen,Dlat,Dlong) 

        outlines.append( linemod ) 

fulloutpath = "geosamples.tsv"
with open( fulloutpath , 'w' ) as fhandle :
    fhandle.writelines(  outlines ) 

sys.exit(-1)

Dcity = {}
Dstate = {}
Dyear = {}
Dlat = {}
Dlong = {} 
fullpath = "baseMmelodias.txt"
with open( fullpath, 'r') as fhandle:
    for line in fhandle:
        basen = line.strip() 
        srch = basen.find("Mmelodia_")
        if( srch > (-1) ) :
            Dcity[basen] = extract_city(basen) 
            Dstate[basen] = extract_state(basen)
            Dyear[basen] = extract_year(basen)
            geoc = geo_coords[ Dcity[basen] ]   
            Dlat[basen] = geoc[0]
            Dlong[basen] = geoc[1] 

fulloutf= "geo_Mmelodia.csv"
with open( fulloutf, 'w') as fhandle:
    fhandle.write( "name,city,state,year,latitude,longitude\n" ) 
    for k in (Dcity.keys()) :
        row = "{},{},{},{},{},{}\n".format( k, Dcity[k] , Dstate[k], Dyear[k], Dlat[k], Dlong[k] )  
        fhandle.write(row) 

    
"""
cities = list( Dcity.values() )
setcities = set(cities) 
uniquecities = list( setcities )
uniquecities.sort() 
for c in uniquecities :
    outc = "\"{}\":(38.260,-76.005),".format(c)
    print(outc) 
"""  
