# Import python modules
import os
import sys

# data handling libraries
import pandas as pd
import numpy as np
import collections as col
import csv
from datetime import datetime, timedelta

# graphical control libraries
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt

# shape and layer libraries
import fiona
from shapely.geometry import MultiPolygon, shape, point, box

# data wrangling libraries
import urllib2
import ftplib
import wget
import bz2
from bs4 import BeautifulSoup as bs

print 'Version 9/8/17 5:36  cb'


def getFullShape(shapefile):
    """
    Generate a MultiPolygon to represent each shape/polygon within the shapefile
    """
    shp = fiona.open(shapefile)
    mp = [shape(pol['geometry']) for pol in shp]
    mp = MultiPolygon(mp)
    shp.close()
    return(mp)
    
    
def getShapeBbox(polygon):
    """
    Generate a geometric box to represent the bounding box for the polygon, shapefile connection, or MultiPolygon
    """
    # identify the cardinal bounds
    minx, miny, maxx, maxy = polygon.bounds
    bbox = box(minx, miny, maxx, maxy, ccw=True)
    return(bbox)


def readShapefileTable(shapefile):
    """
    read in the datatable captured within the shapefile properties
    """
    shp = fiona.open(shapefile)
    centroid = [eachpol['properties'] for eachpol in shp]
    cent_df = pd.DataFrame.from_dict(centroid, orient='columns')
    shp.close()
    return(cent_df)

def filterPointsinShape(shape, 
                        points_lat, points_lon, points_elev=None, 
                        buffer_distance=0.06, buffer_resolution=16, labels=['LAT', 'LONG_', 'ELEV']):
    """
    filter for datafiles that can be used
    """
    # add buffer region
    region = shape.buffer(buffer_distance, resolution=buffer_resolution)
    
    # Intersection each coordinate with the region
    limited_list = []
    for lon, lat, elev in zip(points_lon, points_lat, points_elev):
        gpoint = point.Point(lon, lat)
        if gpoint.intersects(region):
            limited_list.append([lat, lon, elev])
    maptable = pd.DataFrame.from_records(limited_list, columns=labels)
    print('Number of gridded points/files: '+ str(len(maptable)))
    return(maptable)


# scrape the datafiles from a url of interest
def scrapeftp(url):
    # grab the html of the ftp url
    page = urllib2.urlopen(url).readlines()
        
    # grapple each line for the filename segment on the right side
    filenames = [eachline.replace('\r\n','').rsplit(' ',1)[1] for eachline in page]

    # separate the lon and lat points from the filename
    lonlat = [eachfile.replace('.bz2','').rsplit('_',2) for eachfile in filenames]
    filename_df = pd.DataFrame.from_records(lonlat, columns = ['prefix','LONG_','LAT'])
    return(filename_df)

# scrape the datafiles from a url of interest
def scrapeurl(url, startswith=None, hasKeyword=None):
    # grab the html of the url, and prettify the html structure
    page = urllib2.urlopen(url).read()
    page_soup = bs(page, 'lxml')
    page_soup.prettify()

    # loop through and filter the hyperlinked lines
    if startswith is not None:
        temp = [anchor['href'] for anchor in page_soup.findAll('a', href=True) if anchor['href'].startswith(startswith)]
    else:
        temp = [anchor['href'] for anchor in page_soup.findAll('a', href=True) if hasKeyword in anchor['href']]

    # convert to dataframe then separate the lon and lat coordinate values
    temp = pd.DataFrame(temp, columns = ['filenames'])
    temp[['filetype','Lat','Lon']] = temp['filenames'].apply(lambda x: pd.Series(str(x).rsplit('_', 2)))
    temp['Lon'] = temp.Lon.astype(float)
    temp['Lat'] = temp.Lat.astype(float)
    return(temp)


def treatgeoself(shapefile, NAmer, folder_path=os.getcwd(), outfilename='monkeysonatree.csv', buffer_distance=0.06):
    """
    TreatGeoSelf to some [data] lovin'!
    
    :param shapefile:
    :param Namer:
    :param Webservice:
    :param file_prefix:
    :param folder_path:
    :param outfilename:
    :param buffer_distance:
    """
    
    # read shapefile into a multipolygon shape-object
    shape_mp = getFullShape(shapefile)

    # read in the North American continental DEM points for the station elevations
    NAmer_datapoints = readShapefileTable(NAmer).rename(columns={'Lat':'LAT','Long':'LONG_','Elev':'ELEV'})
    
    # generate maptable
    maptable = filterPointsinShape(shape_mp,
                                   points_lat=NAmer_datapoints.LAT,
                                   points_lon=NAmer_datapoints.LONG_,
                                   points_elev=NAmer_datapoints.ELEV,
                                   buffer_distance=buffer_distance,
                                   buffer_resolution=16,
                                   labels=['LAT', 'LONG_', 'ELEV'])
    maptable.reset_index(inplace=True)
    maptable = maptable.rename(columns={"index":"FID"})
    print(maptable.shape)
    print(maptable.tail())
    
    # print the mappingfile
    mappingfile=os.path.join(folder_path, outfilename)
    maptable.to_csv(mappingfile, sep=',', header=True, index=False)
    return mappingfile

# ## Define internal operations for downloading data

def mapContentFolder(resid):
    path = os.path.join('/home/jovyan/work/notebooks/data', str(resid), str(resid), 'data/contents')
    return path


# read in reference locations
def read_in_longlats(mappingfile):
    maptable=[]
    with open(mappingfile, 'rb') as csvfile:
        longlat = csv.reader(csvfile, delimiter=',')
        for row in longlat:
            maptable.append(row)
    csvfile.close()
    return(maptable)

# ### CIG (DHSVM)-oriented functions

# index and extract longitude and latitude points
def compile_bc_Livneh2013_locations(maptable):
    # index the lat long fields
    latitude=maptable[0].index('LAT')
    longitude=maptable[0].index('LONG_')

    # compile a list of file urls for Livneh et al. 2013 (CIG)
    locations2013=[]
    for row in maptable:
        if maptable.index(row)!=0:
            basename='_'.join(['data',row[latitude], row[longitude]])
            url=['http://cses.washington.edu/rocinante/Livneh/bcLivneh_WWA_2013/forcings_ascii/',basename]
            locations2013.append(''.join(url))
    return(locations2013)

# index and extract longitude and latitude points
def compile_Livneh2013_locations(maptable):
    # index the lat long fields
    latitude=maptable[0].index('LAT')
    longitude=maptable[0].index('LONG_')

    # compile a list of lats and longs for Livneh 2013
    locations2013=[]
    for row in maptable:
        if maptable.index(row)!=0:
            basename='_'.join(['data',row[latitude], row[longitude]])
            url=['http://www.cses.washington.edu/rocinante/Livneh/Livneh_WWA_2013/forcs_dhsvm/',basename]
            locations2013.append(''.join(url))
    return(locations2013)


# ### VIC-oriented functions

# In[3]:

# compile file URLs
def compile_VICASCII_Livneh2015_locations(maptable):
    # index the lat long fields
    latitude=maptable[0].index('LAT')
    longitude=maptable[0].index('LONG_')
    
    # compile the VIC.ASCII data from Livneh 2016
    locations2015=[]
    for row in maptable:
        if maptable.index(row)!=0:
            loci='_'.join(['Fluxes_Livneh_NAmerExt_15Oct2014',row[latitude], row[longitude]])
            url=["ftp://192.12.137.7/pub/dcp/archive/OBS/livneh2014.1_16deg/VIC.ASCII/latitude.",row[latitude],'/',loci,'.bz2']
            locations2015.append(''.join(url))
    return(locations2015)

# compile file URLs
def compile_VICASCII_Livneh2013_USA_locations(maptable):
    # index the lat long fields
    latitude=maptable[0].index('LAT')
    longitude=maptable[0].index('LONG_')
    
    # compile the VIC.ASCII data from Livneh 2013
    locations2013=[]
    for row in maptable:
        if maptable.index(row)!=0:

            loci='_'.join(['VIC_fluxes_Livneh_CONUSExt_v.1.2_2013',row[latitude], row[longitude]])
            url=["ftp://ftp.hydro.washington.edu/pub/blivneh/CONUS/Fluxes.asc.v.1.2.1915.2011.bz2/fluxes.125.120.37.49/",loci,".bz2"]
            locations2013.append(''.join(url))
    return(locations2013)


def compile_VICASCII_Livneh2013_CAN_locations(maptable):
    # index the lat long fields
    latitude=maptable[0].index('LAT')
    longitude=maptable[0].index('LONG_')
    
    # compile the VIC.ASCII data from Livneh 2013
    locations2013=[]
    for row in maptable:
        if maptable.index(row)!=0:

            loci='_'.join(['VIC_fluxes_Livneh_CONUSExt_v.1.2_2013',row[latitude], row[longitude]])
            url=["ftp://ftp.hydro.washington.edu/pub/blivneh/CONUS/Fluxes.asc.v.1.2.1915.2011.bz2/fluxes.canada.columbia/",loci,".bz2"]
            locations2013.append(''.join(url))
    return(locations2013)


# ### Climate (Meteorological observations)-oriented functions

# In[ ]:

# index and extract longitude and latitude points for Livneh 2013
def compile_dailyMET_Livneh2013_locations(maptable):
    # index the lat long fields
    latitude=maptable[0].index('LAT')
    longitude=maptable[0].index('LONG_')
    
    # compile the daily MET data from Livneh 2013
    locations2013=[]
    for row in maptable:
        if maptable.index(row)!=0:
            loci='_'.join(['Meteorology_Livneh_CONUSExt_v.1.2_2013',row[latitude], row[longitude]])
            url=["ftp://ftp.hydro.washington.edu/pub/blivneh/CONUS/Meteorology.asc.v.1.2.1915.2011.bz2/data.125.120.37.49/",loci,".bz2"]
            locations2013.append(''.join(url))
    return(locations2013)

# index and extract longitude and latitude points for Livneh 2016
def compile_dailyMET_Livneh2015_locations(maptable):
    # index the lat long fields
    latitude=maptable[0].index('LAT')
    longitude=maptable[0].index('LONG_')

    # compile the daily data for Livneh 2016
    locations2015=[]
    for row in maptable:
        if maptable.index(row)!=0:
            loci='_'.join(['Meteorology_Livneh_NAmerExt_15Oct2014',row[latitude], row[longitude]])
            url=["ftp://192.12.137.7/pub/dcp/archive/OBS/livneh2014.1_16deg/ascii/daily/latitude.",row[latitude],"/",loci,".bz2"]
            locations2015.append(''.join(url))
    return(locations2015)


# ### WRF-oriented functions

# compile file URLs  

# index and extract longitude and latitude points for raw WRF NNRP
def compile_wrfnnrp_raw_Salathe2014_locations(maptable):
    # index the lat long fields
    latitude=maptable[0].index('LAT')
    longitude=maptable[0].index('LONG_')

    # compile a list of file urls for Livneh et al. 2013 (CIG)
    locations2014=[]
    for row in maptable:
        if maptable.index(row)!=0:
            basename='_'.join(['data',row[latitude], row[longitude]])
            url=['http://cses.washington.edu/rocinante/WRF/NNRP/vic_16d/WWA_1950_2010/raw/forcings_ascii/',basename]
            locations2014.append(''.join(url))
    return(locations2014)

# index and extract longitude and latitude points for bc WRF NNRP
def compile_wrfnnrp_bc_Salathe2014_locations(maptable):
    # index the lat long fields
    latitude=maptable[0].index('LAT')
    longitude=maptable[0].index('LONG_')

    # compile a list of file urls for Livneh et al. 2013 (CIG)
    locations2014=[]
    for row in maptable:
        if maptable.index(row)!=0:
            basename='_'.join(['data',row[latitude], row[longitude]])
            url=['http://cses.washington.edu/rocinante/WRF/NNRP/vic_16d/WWA_1950_2010/bc/forcings_ascii/',basename]
            locations2014.append(''.join(url))
    return(locations2014)


# ## Data file migration functions

# In[15]:

# check if the destination folder directory exists; if not, create it
def ensure_dir(f):
    if not os.path.exists(f):
        os.makedirs(f)
    os.chdir(f)

# Download the livneh 2013 files to the livneh2013 subdirectory
def wget_download(listofinterest):
           
    # check and download each location point, if it doesn't already exist in the download directory
    for fileurl in listofinterest:
        basename=os.path.basename(fileurl)
        filename=os.path.join(os.getcwd(), basename) # file location on the local directory
        if os.path.isfile(filename):
            print('file already exists: ' + basename)
            continue
        try:
            wget.download(fileurl)
            wget.close()
            print('downloaded: ' + basename)
        except:
            print('File does not exist at this URL: ' + basename)
        
# Download the livneh 2013 files to the livneh2013 subdirectory
def wget_download_one(fileurl):
    # check and download each location point, if it doesn't already exist in the download directory
    basename=os.path.basename(fileurl)
    filename=os.path.join(os.getcwd(), basename) # file location on the local directory
    if os.path.isfile(filename):
        print('file already exists: ' + basename)
    else:
        try:
            wget.download(fileurl)
            print('downloaded: ' + basename)
        except:
            print('File does not exist at this URL: ' + basename)
    
def wget_download_p(listofinterest):
    from multiprocessing import Pool
    pool = Pool(10) # submit 10 at once
    pool.map(wget_download_one, listofinterest)
    pool.close()
    pool.terminate()

# Download and decompress the livneh 2016 files
def ftp_download(listofinterest):
    for loci in listofinterest:
        
        # establish path info
        fileurl=loci.replace('ftp://','') # loci is already the url with the domain already appended
        ipaddress=fileurl.split('/',1)[0] # ip address
        path=os.path.dirname(fileurl.split('/',1)[1]) # folder path
        filename=os.path.basename(fileurl) # filename
        
        if os.path.isfile(filename):
            print('file already exists')
            continue
        
        # download the file from the ftp server
        ftp=ftplib.FTP(ipaddress)
        ftp.login()
        ftp.cwd(path)
        try:
            ftp.retrbinary("RETR " + filename ,open(filename, 'wb').write)
            ftp.close()
            
            # decompress the file
            decompbz2(filename)
        except:
            print('File does not exist at this URL: '+fileurl)
        
        

# Download and decompress the livneh 2016 files
def ftp_download_one(loci):
    # establish path info
    fileurl=loci.replace('ftp://','') # loci is already the url with the domain already appended
    ipaddress=fileurl.split('/',1)[0] # ip address
    path=os.path.dirname(fileurl.split('/',1)[1]) # folder path
    filename=os.path.basename(fileurl) # filename
        
    if os.path.isfile(filename):
        print('file already exists')
    else:        
        # download the file from the ftp server
        ftp=ftplib.FTP(ipaddress)
        ftp.login()
        ftp.cwd(path)
        try:
            ftp.retrbinary("RETR " + filename ,open(filename, 'wb').write)
            ftp.close()
            
            # decompress the file
            decompbz2(filename)
        except:
            print('File does not exist at this URL: '+fileurl)

def ftp_download_p(listofinterest):
    from multiprocessing import Pool
    pool = Pool(10) # submit 10 at once
    pool.map(ftp_download_one, listofinterest)
    pool.close()
    pool.terminate()
    
# unzip the file
def decompbz2(filename):
    with open(filename.split(".bz2",1)[0], 'wb') as new_file, open(filename, 'rb') as zipfile:
        decompressor = bz2.BZ2Decompressor()
        for data in iter(lambda : zipfile.read(100 * 1024), b''):
            new_file.write(decompressor.decompress(data))
    os.remove(filename)
    zipfile.close()
    new_file.close()
    print(os.path.splitext(filename)[0] + ' unzipped')


# ## Wrapper scripts

# ### Get Daily Meteorological data from Livneh 2013

# read in the longitude and latitude points from the reference mapping file
def getClimateData_DailyVIC_USA_livneh2013(homedir, mappingfile):
    
    # generate table of lats and long coordinates
    maptable = read_in_longlats(mappingfile)
    
    # compile the longitude and latitude points
    dailyVIClocations2013 = compile_VICASCII_Livneh2013_USA_locations(maptable)

    # check and generate VIC_ASCII Flux model livneh 2016 data directory
    filedir=homedir+'/livneh2013/Daily_VIC_1915_2011/'
    ensure_dir(filedir)

    # Download the livneh 2016 VIC_ASCII Flux model data files
    ftp_download_p(dailyVIClocations2013)
    os.chdir(homedir)
    return(filedir)

def getClimateData_DailyVIC_CAN_livneh2013(homedir, mappingfile):
    
    # generate table of lats and long coordinates
    maptable = read_in_longlats(mappingfile)
    
    # compile the longitude and latitude points
    dailyVIClocations2013 = compile_VICASCII_Livneh2013_CAN_locations(maptable)

    # check and generate VIC_ASCII Flux model livneh 2016 data directory
    filedir=homedir+'/livneh2013/Daily_VIC_1915_2011/'
    ensure_dir(filedir)

    # Download the livneh 2016 VIC_ASCII Flux model data files
    ftp_download_p(dailyVIClocations2013)
    os.chdir(homedir)
    return(filedir)

# read in the longitude and latitude points from the reference mapping file
def getClimateData_DailyMET_livneh2013(homedir, mappingfile):
    
    # generate table of lats and long coordinates
    maptable = read_in_longlats(mappingfile)
    
    # compile the longitude and latitude points
    dailyMETlocations2013 = compile_dailyMET_Livneh2013_locations(maptable)

    # check and generate VIC_ASCII Flux model livneh 2016 data directory
    filedir=homedir+'/livneh2013/Daily_MET_1915_2011/'
    ensure_dir(filedir)

    # Download the livneh 2016 VIC_ASCII Flux model data files
    ftp_download_p(dailyMETlocations2013)
    os.chdir(homedir)
    return(filedir)

def getClimateData_DailyMET_livneh2015(homedir, mappingfile):
    
    # generate table of lats and long coordinates
    maptable = read_in_longlats(mappingfile)
    
    # compile the longitude and latitude points
    dailyMETlocations2015 = compile_dailyMET_Livneh2015_locations(maptable)

    # check and generate VIC_ASCII Flux model livneh 2016 data directory
    filedir=homedir+'/livneh2015/Daily_MET_1950_2013/'
    ensure_dir(filedir)

    # Download the livneh 2016 VIC_ASCII Flux model data files
    ftp_download_p(dailyMETlocations2015)
    os.chdir(homedir)
    return(filedir)


def getClimateData_DailyMET_bcLivneh2013(homedir, mappingfile):
    
    # read in the longitude and latitude points from the reference mapping file
    maptable = read_in_longlats(mappingfile)
    
    # compile the longitude and latitude points
    dailyMET_bcLivneh_locations2013 = compile_bc_Livneh2013_locations(maptable)

    # check and generate baseline_corrected livneh2013 data directory
    filedir=homedir+'/livneh2013/Daily_MET_1915_2011/bc/'
    ensure_dir(filedir)

    # download the data
    wget_download_p(dailyMET_bcLivneh_locations2013)
    os.chdir(homedir)
    return(filedir)

# read in the longitude and latitude points from the reference mapping file
def getClimateData_DailyVIC_livneh2015(homedir, mappingfile):
    
    # generate table of lats and long coordinates
    maptable = read_in_longlats(mappingfile)
    
    # compile the longitude and latitude points
    dailyVIClocations2015 = compile_VICASCII_Livneh2015_locations(maptable)

    # check and generate VIC_ASCII Flux model livneh 2016 data directory
    filedir=homedir+'/livneh2015/Daily_VIC_1950_2013/'
    ensure_dir(filedir)

    # Download the livneh 2016 VIC_ASCII Flux model data files
    ftp_download_p(dailyVIClocations2015)
    os.chdir(homedir)
    return(filedir)


# ### Get Daily Meteorological data from Livneh 2015
# formerly 'getClimateData_DailyMET_Livneh2016.py'

# In[ ]:

# read in the longitude and latitude points from the reference mapping file
def getClimateData_DailyMET_rawWRF(homedir, mappingfile):
    
    # read in the longitude and latitude points from the reference mapping file
    maptable = read_in_longlats(mappingfile)
    
    # compile the longitude and latitude points
    dailyMET_WRF_locations2014 = compile_wrfnnrp_raw_Salathe2014_locations(maptable)

    # check and generate baseline_corrected livneh2013 data directory
    filedir=homedir+'/Salathe2014/WWA_1950_2010/raw/'
    ensure_dir(filedir)

    # download the data
    wget_download_p(dailyMET_WRF_locations2014)
    os.chdir(homedir)
    return(filedir)


def getClimateData_DailyMET_bcWRF(homedir, mappingfile):
    
    # read in the longitude and latitude points from the reference mapping file
    maptable = read_in_longlats(mappingfile)
    
    # compile the longitude and latitude points
    dailyMET_WRF_locations2014 = compile_wrfnnrp_bc_Salathe2014_locations(maptable)

    # check and generate baseline_corrected livneh2013 data directory
    filedir=homedir+'/Salathe2014/WWA_1950_2010/bc/'
    ensure_dir(filedir)

    # download the data
    wget_download_p(dailyMET_WRF_locations2014)
    os.chdir(homedir)
    return(filedir)

# # Data Processing libraries

# In[ ]:


# ## Define functions for reading-in downloaded files

# In[18]:

# create a list of files with their paths to be added to the HydroShare resource.
def compileContentfiles(directory):
    files = []
    filelist = os.listdir(directory) # the list of filenames
    for eachfile in filelist:
        files.append(directory + eachfile) # filepaths list
    return(files)


# Read in the mappingfile as a data frame
def mappingfileToDF(mappingfile):
    map_df= pd.read_csv(mappingfile)
    map_df.columns = map_df.columns.map(lambda x:x.upper())

    # Select for station, latitude, longitude, and elevation columns
    if 'ELEV' in map_df.columns:
        map_df = map_df[['FID','LAT','LONG_','ELEV']]
    elif 'RASTERVALU' in map_df.columns:
        map_df = map_df[['FID','LAT','LONG_','RASTERVALU']]

    # Rename columns in climate locations dataframe
    map_df.columns = ['station','latitude','longitude','elevation']
    
    # compile summaries
    print map_df[0:4]
    print
    print 'Number of climate stations:', len(map_df)
    print 'Minimum elevation of climate stations:', np.min(map_df.elevation), 'm'
    print 'Mean elevation of climate stations:', int(np.mean(map_df.elevation)), 'm'
    print 'Maximum elevation of climate stations:', np.max(map_df.elevation), 'm'
    
    return(map_df, len(map_df))


# appends the folder path with the list of files within it
def filesWithPath(filedir):
    file_names=[]
    listOfFiles = os.listdir(filedir)
    for eachfile in listOfFiles:
        # logical rule to exclude hidden files
        if not eachfile.startswith('.'):
            file_names.extend([os.path.join(filedir,eachfile)])
    return(file_names)

# define a function to read in the met files. The start and end date of the data must be input into the function since it is not 
# included in the raw data download.
def read_in_all_met_files(file_names, file_start_date, file_end_date, subset_start_date, subset_end_date, n_stations):
    
    #initialize matrix and time sequence
    met_daily=[]

    met_daily_dates=pd.date_range(file_start_date, file_end_date, freq='D') # daily
    
    # Identify separator for files in the dataset
    sample = pd.read_table(file_names[0], header=None, sep='\t', nrows=1)
    if len(sample) != 1:
        separator = '\t'
    else:
        separator = '\s+'
    
    # import data for all climate stations
    for j in range(0, n_stations):
        met_daily.append(pd.read_table(file_names[j], header=None, sep=separator))
        met_daily[j].columns=['precip_mm','tmax_c', 'tmin_c', 'wind_m_s']
        
        # set the index as the date
        met_daily[j].set_index(met_daily_dates, inplace=True)
        met_daily[j] = met_daily[j].ix[subset_start_date:subset_end_date]
        
    # display first 10 lines of the first station's data frame
    print met_daily[0][0:10]
    return(met_daily)
    



# read in a daily streamflow data set
def read_daily_streamflow(file_name, drainage_area_m2, file_colnames=None, delimiter='\t', header='infer'):
    
    # if file_colnames are supplied, use header=None
    if file_colnames is not None:
        header=None
    
    # read in the data
    daily_data=pd.read_table(file_name, delimiter=delimiter, header=header) 
    
    # set columns, if header=None
    if file_colnames is not None:
        daily_data.columns=file_colnames
    else:
        file_colnames=list(daily_data.columns)
        
    # calculate cfs to cms conversion, or vice versa
    if 'flow_cfs' in daily_data.columns:
        flow_cfs=daily_data['flow_cfs']
        flow_cms=flow_cfs/(3.28084**3)
        flow_mmday=flow_cms*1000*3600*24/drainage_area_m2
        
    elif 'flow_cms' in daily_data.columns:
        flow_cms=daily_data['flow_cms']
        flow_cfs=flow_cms*(3.28084**3)
        flow_mmday=flow_cms*1000*3600*24/drainage_area_m2
            
    # determine the datetime
    date_index=[file_colnames.index(each) for each in ['year','month','day']]
    row_dates=pd.to_datetime(daily_data[date_index])
    
    # generate the daily_flow and set the datetime as row indices
    daily_flow=pd.concat([flow_cfs, flow_cms, flow_mmday],axis=1)
    daily_flow.set_index(row_dates, inplace=True)
    daily_flow.columns=['flow_cfs', 'flow_cms', 'flow_mmday']
    return(daily_flow)

# read in a daily precipitation data set
def read_daily_precip(file_name, file_colnames=None, header='infer', delimiter='\s+'):
    
    # if file_colnames are supplied, use header=None
    if file_colnames is not None:
        header=None
    
    # read in the data
    daily_data=pd.read_table(file_name, delimiter=delimiter, header=header) 
    
    # set columns, if header=None
    if file_colnames is not None:
        daily_data.columns=file_colnames
    else:
        file_colnames=list(daily_data.columns)
    
    # calculate cfs to cms conversion, or vice versa
    if 'precip_m' in daily_data.columns:
        precip_m=daily_data['precip_m']
        precip_mm=precip_m*1000
    
    # determine the datetime
    date_index=[file_colnames.index(each) for each in ['year','month','day']]
    row_dates=pd.to_datetime(daily_data[date_index])
    
    # generate the daily_flow and set the datetime as row indices
    daily_precip=pd.concat([precip_m, precip_mm],axis=1)
    daily_precip.set_index(row_dates, inplace=True)
    daily_precip.columns=['precip_m', 'precip_mm']
    return(daily_precip)

# read in a daily SNOTEL observation data set
def read_daily_snotel(file_name, file_colnames=None, usecols=None, delimiter=',', header='infer'):
    
    # if file_colnames are supplied, use header=None
    if file_colnames is not None:
        header=None
    
    # read in the data
    daily_data=pd.read_table(file_name, usecols=usecols, header=header, delimiter=delimiter)
    
    # reset the colnames
    daily_data.columns=['Date', 'Tmax_C', 'Tmin_C', 'Tavg_C', 'Precip_mm']
    
    # transform the data
    daily_data['Tmax_C']=(daily_data['Tmax_C'] -32)/1.8
    daily_data['Tmin_C']=(daily_data['Tmin_C'] -32)/1.8
    daily_data['Tavg_C']=(daily_data['Tavg_C'] -32)/1.8
    daily_data['Precip_mm']=daily_data['Precip_mm'] *25.4
    
    # determine the datetime
    row_dates=pd.to_datetime(daily_data.Date)
    
    # generate the daily_flow and set the datetime as row indices
    daily_snotel=daily_data[['Tmax_C', 'Tmin_C', 'Tavg_C', 'Precip_mm']]
    daily_snotel.set_index(row_dates, inplace=True)
    return(daily_snotel)


def read_daily_coop(file_name, file_colnames=None, usecols=None, delimiter=',', header='infer'):
    
    # if file_colnames are supplied, use header=None
    if file_colnames is not None:
        header=None
    
    # read in the data
    daily_data=pd.read_table(file_name, usecols=usecols, header=header, delimiter=delimiter, 
                                 date_parser=lambda x: pd.datetime.strptime(x, '%Y%m%d'), 
                                 parse_dates=[0],
                                 na_values=-9999)
    
    # reset the colnames
    daily_data.columns=['Date', 'Precip_mm','Tmax_C', 'Tmin_C', 'Tavg_C']
    
    # transform the data
    daily_data['Tmax_C']=(daily_data['Tmax_C'] -32)/1.8
    daily_data['Tmin_C']=(daily_data['Tmin_C'] -32)/1.8
    daily_data['Tavg_C']=(daily_data['Tavg_C'] -32)/1.8
    daily_data['Precip_mm']=daily_data['Precip_mm'] *25.4
    
    # determine the datetime
    row_dates=pd.to_datetime(daily_data.Date)
    
    # generate the daily_flow and set the datetime as row indices
    daily_coop=daily_data[['Precip_mm','Tmax_C', 'Tmin_C', 'Tavg_C']]
    daily_coop.set_index(row_dates, inplace=True)
    return(daily_coop)


def aggregate_space_time_sum(VarTable, n_stations, start_date, end_date):
    Var_daily = VarTable.loc[start_date:end_date, range(0,n_stations)]
    
    # Average precipitation per month at each station
    permonth_daily=Var_daily.groupby(pd.TimeGrouper("M")).sum()
    
    # Average precipitation per month averaged at all stations
    meanpermonth_daily=permonth_daily.mean(axis=1)
    
    # Average monthly precipitation averaged at all stations
    meanmonth_daily= meanpermonth_daily.groupby(meanpermonth_daily.index.month).mean()
    
    return(Var_daily,
          permonth_daily,
          meanpermonth_daily,
          meanmonth_daily)

# ### Data Processing functions

# In[ ]:

# Determine dates (rows) and stations (columns). Number of stations 
# is the same for each dataset but number of dates may be different
def generateVarTables (listOfDates, listOfTables, n_stations):
    # NOTE: listOfTable must contain:
    # tmin_c
    # tmax_c
    # precip_mm
    # wind_m_s
    
    len_listOfDates=len(listOfDates) # number of dates
    station_list=range(0, n_stations) # list of stations number
    
    # Create arrays of for each variable of interest (Tmin, Tmax, Precip).
    # Rows are dates of analysis and columns are the station number
    temp_min_np=np.empty([len_listOfDates,n_stations])
    temp_max_np=np.empty([len_listOfDates,n_stations])
    precip_np=np.empty([len_listOfDates,n_stations])
    wind_np=np.empty([len_listOfDates,n_stations])
    
    # fill in each array with values from each station
    for i in station_list:
        temp_min_np[:,i]=listOfTables[i].tmin_c.values.astype(float)
        temp_max_np[:,i]=listOfTables[i].tmax_c.values.astype(float)
        precip_np[:,i]=listOfTables[i].precip_mm.values.astype(float)
        wind_np[:,i]=listOfTables[i].wind_m_s.values.astype(float)
        
    # generate each variable dataframe with rows as dates and columns as stations
    temp_min_df=pd.DataFrame(temp_min_np, columns=station_list, index=listOfDates)    
    temp_max_df=pd.DataFrame(temp_max_np, columns=station_list, index=listOfDates)    
    precip_df=pd.DataFrame(precip_np, columns=station_list, index=listOfDates)    
    wind_df=pd.DataFrame(wind_np, columns=station_list, index=listOfDates)
    
    # Create average temperature data frame as the average of Tmin and Tmax
    temp_avg_df=pd.DataFrame((temp_min_np+temp_max_np)/2, columns=station_list, index=listOfDates)
    
    # generate each variable dataframe with rows as dates and columns as stations
    
    
    return(temp_min_df, temp_max_df, precip_df, wind_df, temp_avg_df)

# compare two date sets for the start and end of the overlapping dates
def overlappingDates(date_set1, date_set2):
    # find recent date
    if date_set1[0] > date_set2[0]:
        start_date = date_set1[0]
    else:
        start_date = date_set2[0]
    
    # find older date
    if date_set1[-1] < date_set2[-1]:
        end_date = date_set1[-1]
    else:
        end_date = date_set2[-1]
    return(start_date, end_date)

# Calculate means by 8 different methods
def multigroupMeans(VarTable, n_stations, start_date, end_date):
    Var_daily = VarTable.loc[start_date:end_date, range(0,n_stations)]
    
    # e.g., Mean monthly temperature at each station
    month_daily=Var_daily.groupby(Var_daily.index.month).mean() # average monthly minimum temperature at each station
    
    # e.g., Mean monthly temperature averaged for all stations in analysis
    meanmonth_daily=month_daily.mean(axis=1)
    
    # e.g., Mean monthly temperature for minimum and maximum elevation stations
    meanmonth_min_maxelev_daily=Var_daily.loc[:,analysis_elev_max_station].groupby(Var_daily.index.month).mean()
    meanmonth_min_minelev_daily=Var_daily.loc[:,analysis_elev_min_station].groupby(Var_daily.index.month).mean()
    
    # e.g., Mean annual temperature
    year_daily=Var_daily.groupby(Var_daily.index.year).mean()
    
    # e.g., mean annual temperature each year for all stations
    meanyear_daily=year_daily.mean(axis=1)
    
    # e.g., mean annual min temperature for all years, for all stations
    meanallyear_daily=np.nanmean(meanyear_daily)
    
    # e.g., anomoly per year compared to average
    anom_year_daily=meanyear_daily-meanallyear_daily
    
    return(month_daily, 
           meanmonth_daily, 
           meanmonth_min_maxelev_daily, 
           meanmonth_min_minelev_daily, 
           year_daily, 
           meanyear_daily, 
           meanallyear_daily,
           anom_year_daily)

def specialTavgMeans(VarTable):
    Var_daily = VarTable.loc[start_date:end_date, range(0,n_stations)]
    
    # Average temperature for each month at each station
    permonth_daily=Var_daily.groupby(pd.TimeGrouper("M")).mean()
    
    # Average temperature each month averaged at all stations
    meanpermonth_daily=permonth_daily.mean(axis=1)
    
    # Average monthly temperature for all stations
    meanallpermonth_daily=meanpermonth_daily.mean(axis=0)
    
    # anomoly per year compared to average
    anom_month_daily=(meanpermonth_daily-meanallpermonth_daily)/1000
    
    return(permonth_daily,
          meanpermonth_daily,
          meanallpermonth_daily,
          anom_month_daily)

#HIDE THIS CODE IN OGH

# ADD NOTES TO EVERY CELL ON WHAT HAPPENED - DID IT RUN SUCCCESSFULLY


# Calculate means by 8 different methods
def aggregate_space_time_average(VarTable, n_stations, elev_min_station, elev_mid_station,elev_max_station, start_date, end_date):
    Var_daily = VarTable.loc[start_date:end_date, range(0,n_stations)]
    
    # e.g., Mean monthly temperature at each station
    month_daily=Var_daily.groupby(Var_daily.index.month).mean() 
    
    # e.g., Mean monthly temperature averaged for all stations in analysis
    meanmonth_daily=Var_daily.loc[:,elev_max_station+elev_mid_station+elev_min_station].groupby(Var_daily.index.month).mean().mean(axis=1)
                
    # e.g., Mean monthly temperature for minimum and maximum elevation stations
    meanmonth_maxelev_daily=Var_daily.loc[:,elev_max_station].groupby(Var_daily.index.month).mean().mean(axis=1)
    meanmonth_midelev_daily=Var_daily.loc[:,elev_mid_station].groupby(Var_daily.index.month).mean().mean(axis=1)
    meanmonth_minelev_daily=Var_daily.loc[:,elev_min_station].groupby(Var_daily.index.month).mean().mean(axis=1)   
         
                     
    # e.g., Mean annual temperature
    year_daily=Var_daily.groupby(Var_daily.index.year).mean()
    
    # e.g., mean annual temperature each year for all stations
    meanyear_daily=year_daily.mean(axis=1)
    
    # e.g., mean annual min temperature for all years, for all stations
    meanallyear_daily=np.nanmean(meanyear_daily)
    
    # e.g., anomoly per year compared to average
    anom_year_daily=meanyear_daily-meanallyear_daily
    
    return(month_daily, 
           meanmonth_daily,
           meanmonth_maxelev_daily, 
           meanmonth_midelev_daily,
           meanmonth_minelev_daily, 
           year_daily, 
           meanyear_daily, 
           meanallyear_daily,
           anom_year_daily)

def aggregate_space_time_sum(VarTable, n_stations, start_date, end_date):
    Var_daily = VarTable.loc[start_date:end_date, range(0,n_stations)]
    
    # Average precipitation per month at each station
    permonth_daily=Var_daily.groupby(pd.TimeGrouper("M")).sum()
    
    # Average precipitation per month averaged at all stations
    meanpermonth_daily=permonth_daily.mean(axis=1)
    
    # Average monthly precipitation averaged at all stations
    meanmonth_daily= meanpermonth_daily.groupby(meanpermonth_daily.index.month).mean()
    
    return(Var_daily,
          permonth_daily,
          meanpermonth_daily,
          meanmonth_daily)

def specialTavgMeans(VarTable):
    Var_daily = VarTable.loc[start_date:end_date, range(0,n_stations)]
    
    # Average temperature for each month at each station
    permonth_daily=Var_daily.groupby(pd.TimeGrouper("M")).mean()
    
    # Average temperature each month averaged at all stations
    meanpermonth_daily=permonth_daily.mean(axis=1)
    
    # Average monthly temperature for all stations
    meanallpermonth_daily=meanpermonth_daily.mean(axis=0)
    
    # anomoly per year compared to average
    anom_month_daily=(meanpermonth_daily-meanallpermonth_daily)/1000
    
    return(permonth_daily,
          meanpermonth_daily,
          meanallpermonth_daily,
          anom_month_daily)

def plotTavg(dictionary, loc_name, start_date, end_date):
    # Plot 1: Monthly temperature analysis of Livneh data
    if 'meanmonth_temp_avg_liv2013_met_daily' and 'meanmonth_temp_avg_wrf2014_met_daily' not in dictionary.keys():
        pass

    # generate month indices
    wy_index=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    wy_numbers=[10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    month_strings=[ 'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sept']
    
    # initiate the plot object
    fig, ax=plt.subplots(1,1,figsize=(10, 6))

    if 'meanmonth_temp_avg_liv2013_met_daily' in dictionary.keys():
        # Liv2013
        plt.plot(wy_index, dictionary['meanmonth_maxelev_temp_avg_liv2013_met_daily'][wy_numbers],'r*--',linewidth=1, label='Liv Tavg- Max Elev='+str(dictionary['analysis_elev_max_cutoff'])+"-"+str(dictionary['analysis_elev_max'])+'m')
        
        plt.plot(wy_index, dictionary['meanmonth_midelev_temp_avg_liv2013_met_daily'][wy_numbers],'r-', linewidth=1, label='Liv Tavg- Mid Elev='+str(dictionary['analysis_elev_min_cutoff'])+"-"+str(dictionary['analysis_elev_max_cutoff'])+'m')
        
        plt.plot(wy_index, dictionary['meanmonth_minelev_temp_avg_liv2013_met_daily'][wy_numbers],'rX--',linewidth=1, label='Liv Tavg- Min Elev='+str(dictionary['analysis_elev_min'])+"-"+str(dictionary['analysis_elev_min_cutoff'])+'m')
    
    
    if 'meanmonth_temp_avg_wrf2014_met_daily' in dictionary.keys():
        # WRF2014
        plt.plot(wy_index, dictionary['meanmonth_maxelev_temp_avg_wrf2014_met_daily'][wy_numbers],'b^--',linewidth=1, label='WRF Tavg- Max Elev='+str(dictionary['analysis_elev_max_cutoff'])+"-"+str(dictionary['analysis_elev_max'])+'m')
        
        plt.plot(wy_index, dictionary['meanmonth_midelev_temp_avg_wrf2014_met_daily'][wy_numbers],'b-',linewidth=1, label='WRF Tavg- Mid Elev='+str(dictionary['analysis_elev_min_cutoff'])+"-"+str(dictionary['analysis_elev_max_cutoff'])+'m')
        
        plt.plot(wy_index, dictionary['meanmonth_minelev_temp_avg_wrf2014_met_daily'][wy_numbers],'bo--',linewidth=1, label='WRF Tavg- Min Elev='+str(dictionary['analysis_elev_min'])+"-"+str(dictionary['analysis_elev_min_cutoff'])+'m')

    if 'meanmonth_temp_avg_livneh2013_wrf2014bc_met_daily' in dictionary.keys():
        # WRF2014
        plt.plot(wy_index, dictionary['meanmonth_maxelev_temp_avg_livneh2013_wrf2014bc_met_daily'][wy_numbers],'g^--',linewidth=1, label='WRFbc Tavg- Max Elev='+str(dictionary['analysis_elev_max_cutoff'])+"-"+str(dictionary['analysis_elev_max'])+'m')
        
        plt.plot(wy_index, dictionary['meanmonth_midelev_temp_avg_livneh2013_wrf2014bc_met_daily'][wy_numbers],'g-',linewidth=1, label='WRFbc Tavg- Mid Elev='+str(dictionary['analysis_elev_min_cutoff'])+"-"+str(dictionary['analysis_elev_max_cutoff'])+'m')
        
        plt.plot(wy_index, dictionary['meanmonth_minelev_temp_avg_livneh2013_wrf2014bc_met_daily'][wy_numbers],'go--',linewidth=1, label='WRFbc Tavg- Min Elev='+str(dictionary['analysis_elev_min'])+"-"+str(dictionary['analysis_elev_min_cutoff'])+'m')

        
    # add reference line at y=0
    plt.plot([1, 12],[0, 0], 'k-',linewidth=1)

    plt.ylabel('Temperature (deg C)',fontsize=14)
    plt.xlabel('Month',fontsize=14)
    plt.xlim(1,12);
    plt.xticks(wy_index, month_strings);
        
    plt.tick_params(labelsize=12)
    plt.legend(loc='best')
    plt.grid(which='both')
    plt.title(str(loc_name)+'\nAverage Temperature\n Years: '+str(start_date.year)+'-'+str(end_date.year)+'; Elevation: '+str(dictionary['analysis_elev_min'])+'-'+str(dictionary['analysis_elev_max'])+'m', fontsize=16)
    
    plt.savefig('avg_monthly_temp'+str(loc_name)+'.png')
    plt.show()
    
def plotPavg(dictionary, loc_name, start_date, end_date):
    # Plot 1: Monthly temperature analysis of Livneh data
    if 'meanmonth_precip_liv2013_met_daily' and 'meanmonth_precip_wrf2014_met_daily' not in dictionary.keys():
        pass

    # generate month indices
    wy_index=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    wy_numbers=[10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    month_strings=[ 'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sept']
    
    # initiate the plot object
    fig, ax=plt.subplots(1,1,figsize=(10, 6))

    if 'meanmonth_precip_liv2013_met_daily' in dictionary.keys():
        # Liv2013

        plt.plot(wy_index, dictionary['meanmonth_maxelev_precip_liv2013_met_daily'][wy_numbers],'r^--',linewidth=1, label='Liv Precip- Max Elev='+str(dictionary['analysis_elev_max_cutoff'])+"-"+str(dictionary['analysis_elev_max'])+'m')
        plt.plot(wy_index, dictionary['meanmonth_midelev_precip_liv2013_met_daily'][wy_numbers],'r-', linewidth=1, label='Liv Precip- Mid Elev='+str(dictionary['analysis_elev_min_cutoff'])+"-"+str(dictionary['analysis_elev_max_cutoff'])+'m')  
        plt.plot(wy_index, dictionary['meanmonth_minelev_precip_liv2013_met_daily'][wy_numbers],'ro--',linewidth=1, label='Liv Precip- Min Elev='+str(dictionary['analysis_elev_min'])+"-"+str(dictionary['analysis_elev_min_cutoff'])+'m')
    
    if 'meanmonth_temp_avg_wrf2014_met_daily' in dictionary.keys():
        # WRF2014

        plt.plot(wy_index, dictionary['meanmonth_maxelev_precip_wrf2014_met_daily'][wy_numbers],'b^--',linewidth=1, label='WRF Precip- Max Elev='+str(dictionary['analysis_elev_max_cutoff'])+"-"+str(dictionary['analysis_elev_max'])+'m')
        plt.plot(wy_index, dictionary['meanmonth_midelev_precip_wrf2014_met_daily'][wy_numbers],'b-',linewidth=1, label='WRF Precip- Mid Elev='+str(dictionary['analysis_elev_min_cutoff'])+"-"+str(dictionary['analysis_elev_max_cutoff'])+'m')
        plt.plot(wy_index, dictionary['meanmonth_minelev_precip_wrf2014_met_daily'][wy_numbers],'bo--',linewidth=1, label='WRF Precip- Min Elev='+str(dictionary['analysis_elev_min'])+"-"+str(dictionary['analysis_elev_min_cutoff'])+'m')

    if 'meanmonth_temp_avg_livneh2013_wrf2014bc_met_daily' in dictionary.keys():
        # WRF2014
        plt.plot(wy_index, dictionary['meanmonth_maxelev_precip_livneh2013_wrf2014bc_met_daily'][wy_numbers],'g^--',linewidth=1, label='WRFbc Precip- Max Elev='+str(dictionary['analysis_elev_max_cutoff'])+"-"+str(dictionary['analysis_elev_max'])+'m')
        
        plt.plot(wy_index, dictionary['meanmonth_midelev_precip_livneh2013_wrf2014bc_met_daily'][wy_numbers],'g-',linewidth=1, label='WRFbc Precip- Mid Elev='+str(dictionary['analysis_elev_min_cutoff'])+"-"+str(dictionary['analysis_elev_max_cutoff'])+'m')
        
        plt.plot(wy_index, dictionary['meanmonth_minelev_precip_livneh2013_wrf2014bc_met_daily'][wy_numbers],'go--',linewidth=1, label='WRFbc Precip- Min Elev='+str(dictionary['analysis_elev_min'])+"-"+str(dictionary['analysis_elev_min_cutoff'])+'m')

    # add reference line at y=0
    plt.plot([1, 12],[0, 0], 'k-',linewidth=1)

    plt.ylabel('Precip (mm)',fontsize=14)
    plt.xlabel('Month',fontsize=14)
    plt.xlim(1,12);
    plt.xticks(wy_index, month_strings);
        
    plt.tick_params(labelsize=12)
    plt.legend(loc='best')
    plt.grid(which='both')
    plt.title(str(loc_name)+'\nAverage Precipitation\n Years: '+str(start_date.year)+'-'+str(end_date.year)+'; Elevation: '+str(dictionary['analysis_elev_min'])+'-'+str(dictionary['analysis_elev_max'])+'m', fontsize=16)
    plt.savefig('avg_monthly_precip'+str(loc_name)+'.png')
    plt.show()
    
def sand_between_your_toes(gridclim_folder,
                  gridclimname,
                  loc_name,
                  mappingfile,
                  min_elev=None,
                  max_elev=None,
                  file_start_date=None, 
                  file_end_date=None,
                  subset_start_date=None,
                  subset_end_date=None,
                  df_dict=None):
    # pipelined operation for assimilating data, processing it, and standardizing the plotting
    
    # generate the climate locations and n_stations
    climate_locations_df, n_stations = mappingfileToDF(mappingfile)
    
    # generate the climate station info
    if min_elev is None:
        min_elev = climate_locations_df.elevation.min()
    
    if max_elev is None:
        max_elev = climate_locations_df.elevation.max()
    
    # take liv2013 date set date range as default if file reference dates are not given
    if file_start_date is None:
        file_start_date = datetime.datetime(1915,1,1)
        
    if file_end_date is None:
        file_end_date = datetime.datetime(2011,12,31)
        
    # take all dates as default if subset reference dates are not given
    if subset_start_date is None:
        subset_start_date = file_start_date
    
    if subset_end_date is None:
        subset_end_date = file_end_date
        
    if df_dict is None:
        df_dict = dict()
    
    # assemble the stations
    analysis_stations_info = climate_locations_df[(climate_locations_df.elevation >= min_elev) & (climate_locations_df.elevation <= max_elev)].sort_values(by='elevation', ascending=False)

    # the number of stations to include into the top and bottom elevation ranges
    x = np.ceil(len(analysis_stations_info.elevation)*.33)

    # Extract list of station numbers for indexing. Alternative, you can set the list of stations manually!
    analysis_elev_max_station = analysis_stations_info.station.head(int(x)).tolist()
    analysis_elev_min_station = analysis_stations_info.station.tail(int(x)).tolist()
    analysis_elev_mid_station = [k for k in analysis_stations_info.station if k not in analysis_elev_max_station + analysis_elev_min_station]

    analysis_elev_max = analysis_stations_info.elevation.max() # maximum elevaiton of stations in analysis
    analysis_elev_max_cutoff = analysis_stations_info.elevation.head(int(x)).min()
    analysis_elev_min_cutoff = analysis_stations_info.elevation.tail(int(x)).max()    
    analysis_elev_min = analysis_stations_info.elevation.min() # minimum elevaiton of stations in analysis
    
    # create list of dataframe
    all_daily = read_in_all_met_files(file_names=filesWithPath(gridclim_folder),
                                          file_start_date=file_start_date,
                                          file_end_date=file_end_date, 
                                          subset_start_date=subset_start_date,
                                          subset_end_date=subset_end_date, 
                                          n_stations=n_stations)
    
    # create datetime series
    all_daily_dates=list(all_daily[0].index)
    
    # generate the variable dataframes for Livneh 2013
    [temp_min,
     temp_max,
     precip,
     wind,
     temp_avg] = generateVarTables(all_daily_dates, all_daily, n_stations)
    
    # assemble the objects to the dictionary
    df_dict['temp_min_'+gridclimname] = temp_min
    df_dict['temp_max_'+gridclimname] = temp_max
    df_dict['precip_'+gridclimname] = precip
    df_dict['wind_'+gridclimname] = wind
    df_dict['temp_avg_'+gridclimname] = temp_avg
    
    # loop through the dictionary to compute each aggregate_space_time_average object
    for eachvardf in df_dict.keys():
        if eachvardf.endswith(gridclimname): # and not eachvardf.startswith('precip'):
            [month_,
             meanmonth_,
             meanmonth_maxelev_,
             meanmonth_midelev_,
             meanmonth_minelev_,
             year_,
             meanyear_,
             meanallyear_,
             anom_year_] = aggregate_space_time_average(df_dict[eachvardf], n_stations, analysis_elev_min_station, analysis_elev_mid_station, analysis_elev_max_station, subset_start_date, subset_end_date)

            df_dict['month_'+eachvardf] = month_
            df_dict['meanmonth_'+eachvardf] = meanmonth_
            df_dict['meanmonth_maxelev_'+eachvardf] = meanmonth_maxelev_
            df_dict['meanmonth_midelev_'+eachvardf] = meanmonth_midelev_
            df_dict['meanmonth_minelev_'+eachvardf] = meanmonth_minelev_
            df_dict['year_'+eachvardf] = year_
            df_dict['meanyear_'+eachvardf] = meanyear_
            df_dict['meanallyear_'+eachvardf] = meanallyear_
            df_dict['anom_year_'+eachvardf] = anom_year_

    # generate plots
    df_dict['analysis_elev_min'] = analysis_elev_min
    df_dict['analysis_elev_max'] = analysis_elev_max
    df_dict['analysis_elev_max_cutoff'] = analysis_elev_max_cutoff
    df_dict['analysis_elev_min_cutoff'] = analysis_elev_min_cutoff
    #plotTavg(df_dict, loc_name,start_date=subset_start_date, end_date=subset_end_date)
    #plotPavg(df_dict, loc_name,start_date=subset_start_date, end_date=subset_end_date)
    
    return df_dict

def compute_diffs(df_dict, df_str, gridclimname1, gridclimname2, prefix1, prefix2='meanmonth_'):
    #Compute difference between monthly means for some data (e.g,. Temp and precip) for two different gridded datasets (e.g., Liv, WRF)
    
    comp_dict=dict()
    for each1 in prefix1:
        for each2 in prefix2:
            comp_dict[str(each1)+df_str] = df_dict[each2+each1+gridclimname1]-df_dict[each2+each1+gridclimname2]
    return comp_dict

def compute_ratios(df_dict, df_str, gridclimname1, gridclimname2, prefix1, prefix2='meanmonth_'):
    #Compute difference between monthly means for some data (e.g,. Temp and precip) for two different gridded datasets (e.g., Liv, WRF)
    
    comp_dict=dict()
    for each1 in prefix1:
        for each2 in prefix2:
            comp_dict[str(each1)+df_str] = df_dict[each2+each1+gridclimname1]/df_dict[each2+each1+gridclimname2]
    return comp_dict

def compute_elev_diffs(df_dict, df_str, gridclimname1, prefix1, prefix2a='meanmonth_minelev_', prefix2b='meanmonth_maxelev_'):
    comp_dict=dict()
    for each1 in prefix1:
        comp_dict[str(each1)+df_str] = df_dict[prefix2a+each1+gridclimname1]-df_dict[prefix2b+each1+gridclimname1]
    return comp_dict

def monthlyBiasCorrection_deltaTratioP_Livneh_METinput(homedir, mappingfile, BiasCorr,
                                                 lowrange='0to1000m', LowElev=range(0,1000),
                                                 midrange='1000to1500m', MidElev=range(1001,1501),
                                                 highrange='1500to3000m', HighElev=range(1501,3000),
                                                 data_dir=None, file_start_date=None, file_end_date=None):

    np.set_printoptions(precision=3)
    
    # take liv2013 date set date range as default if file reference dates are not given
    if file_start_date is None:
        file_start_date = datetime(1915,1,1)
        
    if file_end_date is None:
        file_end_date = datetime(2011,12,31)
    
    # generate the month vector
    month = pd.date_range(start=file_start_date, end=file_end_date).month
    month = pd.DataFrame({'month':month})
    
    # create NEW directory
    dest_dir = os.path.join(homedir, 'biascorrWRF_liv')
    if not os.path.exists(dest_dir):
        os.mkdir(dest_dir)
        print('destdir created')
    
    # read in the Elevation table
    zdiff = pd.read_table(mappingfile, sep=',', header='infer')
    zdiff = zdiff.rename(columns={'RASTERVALU':'Elev','ELEV':'Elev'})
    zdiff = zdiff[['LAT','LONG_', 'Elev']]
    zdiff['filename'] = zdiff[['LAT','LONG_']].apply(lambda x: '_'.join(['Meteorology_Livneh_CONUSExt_v.1.2_2013',str(x[0]), str(x[1])]), axis=1)
    #print(zdiff[0:10])

    # lapse rate vector by month
    # temperature adjustment vector by month

    # identify the files to read
    print('reading in data_long_lat files')
    data_files = [os.path.join(data_dir,dat) for dat in os.listdir(data_dir) if os.path.basename(dat).startswith('Meteorology_Livneh_CONUSExt_v.1.2_2013')]
    print('done reading data_long_lat files')
    
    # loop through each file
    for eachfile in data_files:
        
        # subset the zdiff table using the eachfile's filename, then assign Elevation to equal the Elev value
        Elevation = zdiff[zdiff['filename']==os.path.basename(eachfile)]['Elev'].reset_index(drop=True)
        print(Elevation)
                
        # decide on the elevation-based Tcorr
        #print('Convert BiasCorr to a df')
        if Elevation.iloc[0] in LowElev:
            BiasCorr_sub = {k: v for k, v in BiasCorr.items() if k.endswith('_'+lowrange)}
            #BiasCorr_sub_df = pd.DataFrame.from_dict(BiasCorr_sub, orient='columns').reset_index()
            #BiasCorr_sub_df.columns = ['month', 'precip', 'Tmax', 'Tmin']
            
        elif Elevation.iloc[0] in MidElev:
            BiasCorr_sub = {k: v for k, v in BiasCorr.items() if k.endswith('_'+midrange)}
            #BiasCorr_sub_df = pd.DataFrame.from_dict(BiasCorr_sub, orient='columns').reset_index()
            #BiasCorr_sub_df.columns = ['month', 'precip', 'Tmax', 'Tmin']
        
        elif Elevation.iloc[0] in HighElev:
            BiasCorr_sub = {k: v for k, v in BiasCorr.items() if k.endswith('_'+highrange)}
            #BiasCorr_sub_df = pd.DataFrame.from_dict(BiasCorr_sub, orient='columns').reset_index()
            #BiasCorr_sub_df.columns = ['month', 'precip', 'Tmax', 'Tmin']
       
        #print('reading in eachfile')
        read_dat = pd.read_table(eachfile, delimiter='\s+', header=None)
        read_dat.columns = ['precip', 'temp_max','temp_min','wind']
        # print('done reading eachfile')

        # extrapolate monthly values for each variable
        for eachvar in ['precip', 'temp_max', 'temp_min']:
            BiasCorr_sub_df = [BiasCorr_sub[eachkey] for eachkey in BiasCorr_sub.keys() if eachkey.startswith(eachvar)]
            
            # subset the column for the eachfile station number
            BiasCorr_sub_df = BiasCorr_sub_df.loc[:,zdiff[zdiff['filename']==eachfile].index]
            BiasCorr_sub_df.columns = ['var']
            
            # regenerate the month
            BiasCorr_sub_df = BiasCorr_sub_df.reset_index().rename(columns={'index':'month'})
            
            # generate s-vectors
            month_obj = month.merge(BiasCorr_sub_df, how='left', on='month')
            
            # generate the s-vector
            s = pd.Series(month_obj.var)
            
            #
            if eachvar=='precip':
                read_dat[eachvar] = np.array(read_dat[eachvar])*np.array(s)
            else:
                read_dat[eachvar] = np.array(read_dat[eachvar])+np.array(s)    
        
        #print('grabbing the S vector of monthlapse after the merge between month and Tcorr_df')
        #print('number of corrections to apply: '+str(len(month)))
        
        # write it out to the new destination location
        read_dat.to_csv(os.path.join(dest_dir, os.path.basename(eachfile)), sep='\t', header=None, index=False)
        print(os.path.join(dest_dir, os.path.basename(eachfile)))
    
    print('mission complete.')
    print('this device will now self-destruct.')
    print('just kidding.')
    
def makebelieve_global(homedir,
                                                 mappingfile,
                                                 BiasCorr,
                                                 lowrange='0to1000m', LowElev=range(0,1000),
                                                 midrange='1000to1500m', MidElev=range(1001,1501),
                                                 highrange='1500to3000m', HighElev=range(1501,3000),
                                                 data_dir=None,
                                                 file_start_date=None,
                                                 file_end_date=None):

    np.set_printoptions(precision=3)
    
    # take liv2013 date set date range as default if file reference dates are not given
    if file_start_date is None:
        file_start_date = datetime(1950,1,1)
        
    if file_end_date is None:
        file_end_date = datetime(2010,12,31)
    
    # generate the month vector
    month = pd.date_range(start=file_start_date, end=file_end_date).month
    month = pd.DataFrame({'month':month})
    
    # create NEW directory
    dest_dir = os.path.join(homedir, 'biascorr_WRF_ltm')
    if not os.path.exists(dest_dir):
        os.mkdir(dest_dir)
        print('destdir created')
    
    # read in the Elevation table
    zdiff = pd.read_table(mappingfile, sep=',', header='infer')
    zdiff = zdiff.rename(columns={'RASTERVALU':'Elev','ELEV':'Elev'})
    zdiff = zdiff[['LAT','LONG_', 'Elev']]
    zdiff['filename'] = zdiff[['LAT','LONG_']].apply(lambda x: '_'.join(['data',str(x[0]), str(x[1])]), axis=1)
    #print(zdiff[0:10])

    # lapse rate vector by month
    # temperature adjustment vector by month

    # identify the files to read
    print('reading in data_long_lat files')
    data_files = [os.path.join(data_dir,dat) for dat in os.listdir(data_dir) if os.path.basename(dat).startswith('data')]
    #print('done reading data_long_lat files')
    
    # loop through each file
    for eachfile in data_files:
        
        # subset the zdiff table using the eachfile's filename, then assign Elevation to equal the Elev value
        Elevation = zdiff[zdiff['filename']==os.path.basename(eachfile)]['Elev'].reset_index(drop=True)
        print(Elevation)
                
        # decide on the elevation-based Tcorr
        #print('Convert BiasCorr to a df')
        if Elevation.iloc[0] in LowElev:
            BiasCorr_sub = {k: v for k, v in BiasCorr.items() if k.endswith('_'+lowrange)}
            BiasCorr_sub_df = pd.DataFrame.from_dict(BiasCorr_sub, orient='columns').reset_index()
            BiasCorr_sub_df.columns = ['month', 'precip', 'Tmax', 'Tmin']
            
        elif Elevation.iloc[0] in MidElev:
            BiasCorr_sub = {k: v for k, v in BiasCorr.items() if k.endswith('_'+midrange)}
            BiasCorr_sub_df = pd.DataFrame.from_dict(BiasCorr_sub, orient='columns').reset_index()
            BiasCorr_sub_df.columns = ['month', 'precip', 'Tmax', 'Tmin']
        
        elif Elevation.iloc[0] in HighElev:
            BiasCorr_sub = {k: v for k, v in BiasCorr.items() if k.endswith('_'+highrange)}
            BiasCorr_sub_df = pd.DataFrame.from_dict(BiasCorr_sub, orient='columns').reset_index()
            BiasCorr_sub_df.columns = ['month', 'precip', 'Tmax', 'Tmin']
       
        print('reading in eachfile')
        read_dat = pd.read_table(eachfile, delimiter='\s+', header=None)
        read_dat.columns = ['precip', 'Tmax','Tmin','wind']
        print('done reading eachfile')

        # extrapolate monthly values
        month_obj = month.merge(BiasCorr_sub_df, how='left', on='month')
        #print('merged month with Tcorr_df')
        #print(month_obj.head(35))
        
        # generate s-vectors
        s1 = pd.Series(month_obj.Tmin)
        s2 = pd.Series(month_obj.Tmax)
        s3 = pd.Series(month_obj.precip)
        
        #print('grabbing the S vector of monthlapse after the merge between month and Tcorr_df')
        #print('number of corrections to apply: '+str(len(month)))

        read_dat['Tmin'] = np.array(read_dat.Tmin)+np.array(s1)
        read_dat['Tmax'] = np.array(read_dat.Tmax)+np.array(s2)
        read_dat['precip'] = np.array(read_dat.precip)*np.array(s3)
        
        # write it out to the new destination location
        read_dat.to_csv(os.path.join(dest_dir, os.path.basename(eachfile)), sep='\t', header=None, index=False)
        print(os.path.join(dest_dir, os.path.basename(eachfile)))
    
    print('mission complete.')
    print('this device will now self-destruct.')
    print('just kidding.')
    
def switchUpVICSoil(input_file=None,
                    output_file='soil',
                    mappingfile=None,
                    homedir=None):
    #Read in table of VIC soil inputs -- assumes all Lat/Long set to zero
    soil_base = pd.read_table(input_file,header=None)

    #Make a list of all lat/long values
    latlong=soil_base.apply(lambda x:tuple([x[2],x[3]]), axis=1)

    #Read in mappingfile from TreatGeoSelf()
    maptable = pd.read_table(mappingfile,sep=",")

    #Make a list Lat/Long files that need to switched up 
    latlong_1=maptable.apply(lambda x:tuple([x[2],x[1]]), axis=1)

    #Switch up from 0 to 1 so VIC will run for this Lat/Long point - print new output file (VIC model input file)
    soil_base[0] = latlong.apply(lambda x: 1 if x in set(latlong_1) else 0)        
    soil_base.to_csv(output_file, header=False,index=False,sep="\t")
    print(str(soil_base[0].sum()) +' VIC grid cells have successfully been switched up.') 
    print('Check your home directory for your new VIC soil model input set to your list of Lat/Long grid centroids.')
    
def makebelieve(homedir, mappingfile, BiasCorr,
                lowrange='0to1000m', LowElev=range(0,1000),
                midrange='1000to1500m', MidElev=range(1001,1501),
                highrange='1500to3000m', HighElev=range(1501,3000),
                data_dir=None, file_start_date=None, file_end_date=None,
                dest_dir_suffix=None):
    np.set_printoptions(precision=6)
   
    # take liv2013 date set date range as default if file reference dates are not given
    if file_start_date is None:
        file_start_date = datetime(1915,1,1)
        
    if file_end_date is None:
        file_end_date = datetime(2011,12,31)
    
    # generate the month vector
    month = pd.date_range(start=file_start_date, end=file_end_date).month
    month = pd.DataFrame({'month':month})
    
    # create NEW directory
    if dest_dir_suffix is None:
        dest_dir_suffix = 'biascorr_output/'
        
    dest_dir = os.path.join(homedir, dest_dir_suffix+'/')
    if not os.path.exists(dest_dir):
        os.mkdir(dest_dir)
        print('destdir created')
    
    # read in the Elevation table
    zdiff = pd.read_table(mappingfile, sep=',', header='infer')
    zdiff = zdiff.rename(columns={'RASTERVALU':'Elev','ELEV':'Elev'})
    zdiff = zdiff[['LAT','LONG_', 'Elev']]
    zdiff['filename'] = zdiff[['LAT','LONG_']].apply(lambda x: '_'.join(['Meteorology_Livneh_CONUSExt_v.1.2_2013',str(x[0]), str(x[1])]), axis=1)
    #print(zdiff[0:10])

    # lapse rate vector by month
    # temperature adjustment vector by month

    # identify the files to read
    print('reading in data_long_lat files')
    data_files = [os.path.join(data_dir,dat) for dat in os.listdir(data_dir) if os.path.basename(dat).startswith('Meteorology_Livneh_CONUSExt_v.1.2_2013')]
    print('done reading data_long_lat files')
    
    # loop through each file
    for eachfile in data_files:
        
        # subset the zdiff table using the eachfile's filename, then assign Elevation to equal the Elev value
        Elevation = zdiff[zdiff['filename']==os.path.basename(eachfile)]['Elev'].reset_index(drop=True)
                
        # decide on the elevation-based Tcorr
        #print('Convert BiasCorr to a df')
        if Elevation.iloc[0] in LowElev:
            BiasCorr_sub = {k: v for k, v in BiasCorr.items() if k.endswith('_'+lowrange)}
            #BiasCorr_sub_df = pd.DataFrame.from_dict(BiasCorr_sub, orient='columns').reset_index()
            #BiasCorr_sub_df.columns = ['month', 'precip', 'Tmax', 'Tmin']
            
        elif Elevation.iloc[0] in MidElev:
            BiasCorr_sub = {k: v for k, v in BiasCorr.items() if k.endswith('_'+midrange)}
            #BiasCorr_sub_df = pd.DataFrame.from_dict(BiasCorr_sub, orient='columns').reset_index()
            #BiasCorr_sub_df.columns = ['month', 'precip', 'Tmax', 'Tmin']
        
        elif Elevation.iloc[0] in HighElev:
            BiasCorr_sub = {k: v for k, v in BiasCorr.items() if k.endswith('_'+highrange)}
            #BiasCorr_sub_df = pd.DataFrame.from_dict(BiasCorr_sub, orient='columns').reset_index()
            #BiasCorr_sub_df.columns = ['month', 'precip', 'Tmax', 'Tmin']
       
        #print('reading in eachfile')
        read_dat = pd.read_table(eachfile, delimiter='\s+', header=None)
        read_dat.columns = ['precip', 'temp_max','temp_min','wind']
        # print('done reading eachfile')

        # extrapolate monthly values for each variable
        for eachvar in ['precip', 'temp_max', 'temp_min']:
            
            for eachkey in BiasCorr_sub.keys():
                if eachkey.startswith(eachvar):
                    BiasCorr_sub_df = BiasCorr_sub[eachkey]
            
            # subset the column for the eachfile station number
            BiasCorr_sub_df = BiasCorr_sub_df.loc[:,zdiff[zdiff['filename']==os.path.basename(eachfile)].index]
            BiasCorr_sub_df.columns = ['var']
                        
            # regenerate the month
            BiasCorr_sub_df = BiasCorr_sub_df.reset_index().rename(columns={'index':'month'})
                        
            # generate s-vectors
            month_obj = month.merge(BiasCorr_sub_df, how='left', on='month')
                        
            # generate the s-vector
            s = month_obj.loc[:,'var']
            
            if eachvar=='precip':
                read_dat[eachvar] = np.multiply(np.array(read_dat[eachvar]), np.array(s))
            else:
                read_dat[eachvar] = np.array(read_dat[eachvar])+np.array(s)    
        
        #print('grabbing the S vector of monthlapse after the merge between month and Tcorr_df')
        #print('number of corrections to apply: '+str(len(month)))
        
        # write it out to the new destination location
        read_dat.to_csv(os.path.join(dest_dir, os.path.basename(eachfile)), sep='\t', header=None, index=False, float_format='%.4f')
        print(os.path.join(dest_dir, os.path.basename(eachfile)))
    
    print('mission complete.')
    print('this device will now self-destruct.')
    print('just kidding.')
    
    return dest_dir

def plot_meanP(dictionary, loc_name, start_date, end_date):
    # Plot 1: Monthly temperature analysis of Livneh data
    if 'meanmonth_precip_liv2013_met_daily' and 'meanmonth_precip_wrf2014_met_daily' not in dictionary.keys():
        pass

    # generate month indices
    wy_index=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    wy_numbers=[10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    month_strings=[ 'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sept']
    
    # initiate the plot object
    fig, ax=plt.subplots(1,1,figsize=(10, 6))

    if 'meanmonth_precip_liv2013_met_daily' in dictionary.keys():
        # Liv2013

        plt.plot(wy_index, dictionary['meanmonth_precip_liv2013_met_daily'][wy_numbers],'r-', linewidth=1, label='Liv Precip')  
    
    if 'meanmonth_precip_wrf2014_met_daily' in dictionary.keys():
        # WRF2014
        plt.plot(wy_index, dictionary['meanmonth_precip_wrf2014_met_daily'][wy_numbers],'b-',linewidth=1, label='WRF Precip')
 
    if 'meanmonth_precip_livneh2013_wrf2014bc_met_daily' in dictionary.keys():
        # WRF2014
        plt.plot(wy_index, dictionary['meanmonth_precip_livneh2013_wrf2014bc_met_daily'][wy_numbers],'g-',linewidth=1, label='WRFbc Precip')
 
    # add reference line at y=0
    plt.plot([1, 12],[0, 0], 'k-',linewidth=1)

    plt.ylabel('Precip (mm)',fontsize=14)
    plt.xlabel('Month',fontsize=14)
    plt.xlim(1,12);
    plt.xticks(wy_index, month_strings);
        
    plt.tick_params(labelsize=12)
    plt.legend(loc='best')
    plt.grid(which='both')
    plt.title(str(loc_name)+'\nAverage Precipitation\n Years: '+str(start_date.year)+'-'+str(end_date.year)+'; Elevation: '+str(dictionary['analysis_elev_min'])+'-'+str(dictionary['analysis_elev_max'])+'m', fontsize=16)
    plt.savefig('monthly_precip'+str(loc_name)+'.png')
    plt.show()
    
def plot_meanTavg(dictionary, loc_name, start_date, end_date):
    # Plot 1: Monthly temperature analysis of Livneh data
    if 'meanmonth_temp_avg_liv2013_met_daily' and 'meanmonth_temp_avg_wrf2014_met_daily' not in dictionary.keys():
        pass

    # generate month indices
    wy_index=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    wy_numbers=[10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    month_strings=[ 'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sept']
    
    # initiate the plot object
    fig, ax=plt.subplots(1,1,figsize=(10, 6))

    if 'meanmonth_temp_avg_liv2013_met_daily' in dictionary.keys():
        # Liv2013

        plt.plot(wy_index, dictionary['meanmonth_temp_avg_liv2013_met_daily'][wy_numbers],'r-', linewidth=1, label='Liv Temp Avg')  
    
    if 'meanmonth_temp_avg_wrf2014_met_daily' in dictionary.keys():
        # WRF2014
        plt.plot(wy_index, dictionary['meanmonth_temp_avg_wrf2014_met_daily'][wy_numbers],'b-',linewidth=1, label='WRF Temp Avg')
 
    if 'meanmonth_temp_avg_livneh2013_wrf2014bc_met_daily' in dictionary.keys():
        # WRF2014
        plt.plot(wy_index, dictionary['meanmonth_temp_avg_livneh2013_wrf2014bc_met_daily'][wy_numbers],'g-',linewidth=1, label='WRFbc Temp Avg')
 
    # add reference line at y=0
    plt.plot([1, 12],[0, 0], 'k-',linewidth=1)

    plt.ylabel('Temp (C)',fontsize=14)
    plt.xlabel('Month',fontsize=14)
    plt.xlim(1,12);
    plt.xticks(wy_index, month_strings);
        
    plt.tick_params(labelsize=12)
    plt.legend(loc='best')
    plt.grid(which='both')
    plt.title(str(loc_name)+'\nAverage Temperature\n Years: '+str(start_date.year)+'-'+str(end_date.year)+'; Elevation: '+str(dictionary['analysis_elev_min'])+'-'+str(dictionary['analysis_elev_max'])+'m', fontsize=16)
    plt.savefig('monthly_Tavg'+str(loc_name)+'.png')
    plt.show()
    
def plot_meanTmin(dictionary, loc_name, start_date, end_date):
    # Plot 1: Monthly temperature analysis of Livneh data
    if 'meanmonth_temp_min_liv2013_met_daily' and 'meanmonth_temp_min_wrf2014_met_daily' not in dictionary.keys():
        pass

    # generate month indices
    wy_index=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    wy_numbers=[10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    month_strings=[ 'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sept']
    
    # initiate the plot object
    fig, ax=plt.subplots(1,1,figsize=(10, 6))

    if 'meanmonth_temp_min_liv2013_met_daily' in dictionary.keys():
        # Liv2013

        plt.plot(wy_index, dictionary['meanmonth_temp_min_liv2013_met_daily'][wy_numbers],'r-', linewidth=1, label='Liv Temp min')  
    
    if 'meanmonth_temp_min_wrf2014_met_daily' in dictionary.keys():
        # WRF2014
        plt.plot(wy_index, dictionary['meanmonth_temp_min_wrf2014_met_daily'][wy_numbers],'b-',linewidth=1, label='WRF Temp min')
 
    if 'meanmonth_temp_min_livneh2013_wrf2014bc_met_daily' in dictionary.keys():
        # WRF2014
        plt.plot(wy_index, dictionary['meanmonth_temp_min_livneh2013_wrf2014bc_met_daily'][wy_numbers],'g-',linewidth=1, label='WRFbc Temp min')
 
    # add reference line at y=0
    plt.plot([1, 12],[0, 0], 'k-',linewidth=1)

    plt.ylabel('Temp (C)',fontsize=14)
    plt.xlabel('Month',fontsize=14)
    plt.xlim(1,12);
    plt.xticks(wy_index, month_strings);
        
    plt.tick_params(labelsize=12)
    plt.legend(loc='best')
    plt.grid(which='both')
    plt.title(str(loc_name)+'\nMinimum Temperature\n Years: '+str(start_date.year)+'-'+str(end_date.year)+'; Elevation: '+str(dictionary['analysis_elev_min'])+'-'+str(dictionary['analysis_elev_max'])+'m', fontsize=16)
    plt.savefig('monthly_Tmin'+str(loc_name)+'.png')
    plt.show()
    
def plot_meanTmax(dictionary, loc_name, start_date, end_date):
    # Plot 1: Monthly temperature analysis of Livneh data
    if 'meanmonth_temp_max_liv2013_met_daily' and 'meanmonth_temp_max_wrf2014_met_daily' not in dictionary.keys():
        pass

    # generate month indices
    wy_index=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    wy_numbers=[10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    month_strings=[ 'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sept']
    
    # initiate the plot object
    fig, ax=plt.subplots(1,1,figsize=(10, 6))

    if 'meanmonth_temp_max_liv2013_met_daily' in dictionary.keys():
        # Liv2013

        plt.plot(wy_index, dictionary['meanmonth_temp_max_liv2013_met_daily'][wy_numbers],'r-', linewidth=1, label='Liv Temp max')  
    
    if 'meanmonth_temp_max_wrf2014_met_daily' in dictionary.keys():
        # WRF2014
        plt.plot(wy_index, dictionary['meanmonth_temp_max_wrf2014_met_daily'][wy_numbers],'b-',linewidth=1, label='WRF Temp max')
 
    if 'meanmonth_temp_max_livneh2013_wrf2014bc_met_daily' in dictionary.keys():
        # WRF2014
        plt.plot(wy_index, dictionary['meanmonth_temp_max_livneh2013_wrf2014bc_met_daily'][wy_numbers],'g-',linewidth=1, label='WRFbc Temp max')
 
    # add reference line at y=0
    plt.plot([1, 12],[0, 0], 'k-',linewidth=1)

    plt.ylabel('Temp (C)',fontsize=14)
    plt.xlabel('Month',fontsize=14)
    plt.xlim(1,12);
    plt.xticks(wy_index, month_strings);
        
    plt.tick_params(labelsize=12)
    plt.legend(loc='best')
    plt.grid(which='both')
    plt.title(str(loc_name)+'\nMaximum Temperature\n Years: '+str(start_date.year)+'-'+str(end_date.year)+'; Elevation: '+str(dictionary['analysis_elev_min'])+'-'+str(dictionary['analysis_elev_max'])+'m', fontsize=16)
    plt.savefig('monthly_Tmax'+str(loc_name)+'.png')
    plt.show()