Instructions for running detectSCV.

detectSCV is an object that combines data storage (organized into a structure) and data functions (object methods). 
It can be used on either argo or meop data after they are formatted correctly (see below).

Follow the initial steps below first, then run detectSCV_master.m, which performs the object methods in order.

Initial steps:
1. Obtain processing routines
- Download the GSW scripts and place them in a folder (http://www.teos-10.org/pubs/gsw/html/gsw_contents.html)
- Update the folder name in detectSCV_settings.m 

2. Obtain climatological product
- If using Andrew's WOA2018.mat, run the following before processing:
	>> fdir = '/directory/where/WOA2018.mat/exists' <<-- update with correct dir
	>> fname = 'WOA2018.mat'
	>> detectSCV.woa18_format(fdir,fname)
- This will create woa18_format.mat which is correctly formatted
- If using a separate climatological product, save the .mat file in the following structure format:
	- (climname).lon   = array of longitudes (1:M)
	- (climname).lat   = array of latitutes  (1:N)
	- (climname).month = array of months     (1:T) -- (1:1:12)
	- (climname).pres  = array of pressures  (1:Z)
	- (climname).temp  = (M x N x Z x T) matrix of temperatures
	- (climname).salt  = (M x N x Z x T) matrix of salinities
- Prefer to have longitudes from 0 - 360, but code should fix it

3. Reformat SealData.mat
- If using Andrew's SealData.mat', run the following before processing:
	>> fdir = '/directory/where/SealData.mat/exists' <<-- update with correct dir
	>> fname = 'SealData.mat'
	>> detectSCV.meop_format(fdir,fname)

4. Update detectSCV.getSettings
- getSettings reproduces settings.mat, which is loaded into objInit
- You can find getSettings at the very end of the detectSCV.m
------------------------------------------------------
-                Data directories                    -
------------------------------------------------------
- settings.().rawdir   = raw data location
- settings.().rawfile  = raw data filename
- settings.().qcdir    = where to put qc'd data
- settings.().qcfile   = what to name qc'd data
- settings.().procdir  = where to put proc'd data
- settings.().procfile = what to name proc'd data
- settings.().anomdir  = where to put anomalies data 
- settings.().anomfile = what to name anomalies data
- settings.().anomclim = what to name float climatologies data
- settings.().iqrdir   = where to put IQR data
- settings.().iqrfile  = what to name IQR data
- settings.().climdir  = where to find climatological data
- settings.().climfile = filename of climatological data
- settings.().climvar  = variable name of climatological data
------------------------------------------------------
-                Script directories                  -
------------------------------------------------------
- settings.().gswdir   = GSW script location
------------------------------------------------------
-                Quality Control                     -
------------------------------------------------------
- settings.().qcmatr = vertical resolution thresholds
	- format: [(lower pressure limit) (upper pressure limit) (minimum resolution)];
	- adding additional rows will add additional thresholds
	- example: [0 700 65; 700 inf 105] rejects profiles with a 65dbar gap in data
		   between 0 - 700 dbar and a 105dbar gap from 700 dbar to max depth 
- settings.().deep    = casts must have data below this pressure
- settings.().shallow = casts must have data above this pressure
- settings.().depths  = casts must have this many depths with good data
------------------------------------------------------
-                Processing                          -
------------------------------------------------------
- settings.().ksrbandwidth   = bandwidth (dbar) to use in ksr.m routine. See detectSCV.ksr.m 
	- ksr.m is used to smooth the buoyancy frequency signal, which can take on 
	- small-scale noise due to instrumentation limitations (accuracy of reported measurements).
- settings.().density_thresh = allowable density inversion (will be sorted).
	- Small density inversions in the surface ocean (i.e. mixed layer) will cause issues
	  during interpolation. This threshold will allow any inversions smaller than this 
	  value to be sorted. In practice, this fixes very small mixed-layer inversions.
------------------------------------------------------
-                Climatology                         -
------------------------------------------------------
- settings.().min_clim_pres = Climatology must reach this pressure
------------------------------------------------------
-                Nearby floats                       -
------------------------------------------------------
- settings.().search_dist = distance (km) to look for nearby casts (in a circle) 
- settings.().search_time = number of days (+-) to look for neartime casts
- settings.().min_nearby  = minimum number of 'good' nearby floats (space, time) for IQR calc
------------------------------------------------------
-                IQR Thresholds                      -
------------------------------------------------------
- settings.().iqr_mult = multiplier used in setting IQR thresholds
	- Typically, IQR is used to identify outliers by calculating 25%, 75% percentiles
	- IQR is defined as the difference between Q3 (75%) and Q1 (25%), or Q3 - Q1
	- An outlier is flagged when data is less than Q1 - 1.5xIQR or greater thatn Q3 + 1.5xIQR
	- iqr_mult is thus defaulted to 1.5
