%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef detectSCV % start define detectSCV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties % start define detectSCV object properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	info;     % provides info on data type and where things are stored
	flags;    % where bad profiles are identified
	profile;  % where loaded and processed data exists
			  % profile.raw    = raw data
			  % profile.qc     = quality controlled data
			  % profile.proc   = processed data (interpolated, extra fields)
			  % profile.anom   = anomalies from climatology data
			  % profile.clim   = climatology for each float
			  % profile.nearby = IDs of nearby floats for IQR calculation
			  % profile.iqr    = interquartile range, Q1 + Q3 values, and
			  %                  IQR thresholds for spiciness/N2 anomalies
	scv;      % where detected SCV data exists
			  % scv.initial    = profiles containing spice + N2 anomalies
			  %	             which exceed IQR thresholds
end % end define object properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods % start define detectSCV object methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	function obj = objInit(obj,choice);
	% ----------------------------------------------------------------------
	% Gather data directories containing raw meop or argo data
	%  
	% usage:
	%	obj = objInit(obj,choice)
	%
	% inputs:
	%       choice = 'argo','meop', or 'all'
	% 	If choosing 'all', note that different settings can be applied to each
	% 	data type (argo, meop), see detectSCV.getSettings 
	%
	% example:
	%	obj = objInit(obj,'argo')
	% ----------------------------------------------------------------------
		
		% Announce routine
		disp(' ');
		disp('------------------------- ');
		disp(['Initiating detectSCV:',choice]);
		disp('------------------------- ');
		disp(' ');

		% Load setings.mat
		try
			load('settings.mat');
		catch
			detectSCV.getSettings;
		end

		% addpath to GSW routines
		addpath(settings.gswdir)
		
		% Get info
		disp('Grabbing settings...');
		if strcmp(choice,'argo');
			obj.info            = settings.argo; 
			obj.info.datatype   = 'argo';
		elseif strcmp(choice,'meop');
			obj.info            = settings.meop; 
			obj.info.datatype   = 'meop';
		elseif strcmp(choice,'all');
			obj.info(1)            = settings.argo; 
			obj.info(1).datatype   = 'argo';
			obj.info(2)            = settings.meop; 
			obj.info(2).datatype   = 'meop';
		end

	end % end object method objInfo

	function obj = objQC(obj)
	% ------------------------------------------------------------
	% Perform quality control routines based on settings
	%
	% usage:
	%	obj = objQC(obj)
	%
	% flags applied:
	% 	.flag.empty      = no data
	%	.flag.shallow    = no data above settings.().shallow
	%       .flag.deep       = no data below settings.().deep
	%       .flag.depths     = less than settings.().depths z-levels
	% ------------------------------------------------------------
		
		% Announce routine
		disp(' ');
		disp('------------------------- ');
		disp('Performing Quality Control');
		disp('------------------------- ');
		disp(' ');

		% Set up datatype loop
		for i = 1:length(obj.info)

			% Load raw data
			disp('Loading raw data...');
			tmpdata = load([obj.info(i).rawdir,obj.info(i).rawfile]);
			if strcmp(obj.info(i).datatype,'argo');
				tmpdata = tmpdata.argo;
			elseif strcmp(obj.info(i).datatype,'meop');
				tmpdata = tmpdata.meop;
			end

			% Set up flags	
			obj.flags(i).empty      = []; % code == 91
			obj.flags(i).deep       = []; % code == 92
			obj.flags(i).shallow    = []; % code == 93
			obj.flags(i).depths     = []; % code == 94
			flag_idx                = zeros(length(tmpdata),1);

			% Separate raw data into different types, perform quality control routines
			disp('Start profile loop...')
			cnt   = 0; % initiate disp_prog
			for j = 1:length(tmpdata)

				% Announce loop progress
				cnt = detectSCV.disp_prog(j,length(tmpdata),cnt);

				% Get meta data from each profile, create ID system
				if strcmp(obj.info(i).datatype,'argo')
					obj.profile(i).meta(j).ID = 1e8 + j;
				elseif strcmp(obj.info(i).datatype,'meop')
					obj.profile(i).meta(j).ID = 2e8 + j;
				end
				obj.profile(i).meta(j).lat  = tmpdata(j).LAT;
				obj.profile(i).meta(j).lon  = tmpdata(j).LON;
				obj.profile(i).meta(j).time = tmpdata(j).TIME;

				% Get raw ctd data (+ QC codes) from each profile
				obj.profile(i).raw(j).ID      = obj.profile(i).meta(j).ID;
				obj.profile(i).raw(j).lon     = tmpdata(j).LON;
				obj.profile(i).raw(j).lat     = tmpdata(j).LAT;
				obj.profile(i).raw(j).time    = tmpdata(j).TIME;
				obj.profile(i).raw(j).temp    = tmpdata(j).TEMP;
				obj.profile(i).raw(j).salt    = tmpdata(j).SALT;
				obj.profile(i).raw(j).pres    = tmpdata(j).PRES;
				obj.profile(i).raw(j).temp_qc = tmpdata(j).TEMP_QC;
				obj.profile(i).raw(j).salt_qc = tmpdata(j).SALT_QC;
				obj.profile(i).raw(j).pres_qc = tmpdata(j).PRES_QC;
				
				% Copy raw to qc
				obj.profile(i).qc(j) = obj.profile(i).raw(j);

				% Check bad flags (3+), make them NaN
				obj.profile(i).qc(j).temp(find(str2num(obj.profile(i).raw(j).temp_qc)>2)) = NaN;
				obj.profile(i).qc(j).salt(find(str2num(obj.profile(i).raw(j).salt_qc)>2)) = NaN;
				obj.profile(i).qc(j).pres(find(str2num(obj.profile(i).raw(j).pres_qc)>2)) = NaN;

				% Remove bad rows
				badline                     = [obj.profile(i).qc(j).pres + ...
							       obj.profile(i).qc(j).temp + ...
						               obj.profile(i).qc(j).salt];
				ind                         = find(isnan(badline)==1);
				obj.profile(i).qc(j).temp(ind) = [];
				obj.profile(i).qc(j).salt(ind) = [];
				obj.profile(i).qc(j).pres(ind) = [];

				% Check for empty profiles
				if isempty(obj.profile(i).qc(j).pres) == 1 | ...
				   isempty(obj.profile(i).qc(j).temp) == 1 | ...
				   isempty(obj.profile(i).qc(j).salt) == 1
					flag_idx(j) = [91];
					continue
				end
				
				% Go through each vertical resolution threshold
				% Ignore any data below a big gap in resolution
				[rws,cls] = size(obj.info(i).qcmatr);
				for k     = 1:rws
					idx  = find(obj.info(i).qcmatr(k,1) <= obj.profile(i).qc(j).pres & ...
					       obj.profile(i).qc(j).pres < obj.info(i).qcmatr(k,2));
					didx = find(diff(obj.profile(i).qc(j).pres(idx))>obj.info(i).qcmatr(k,3));
					if isempty(didx) == 0
						obj.profile(i).qc(j).temp([idx(didx(1))+1]:end) = [];
						obj.profile(i).qc(j).salt([idx(didx(1))+1]:end) = [];
						obj.profile(i).qc(j).pres([idx(didx(1))+1]:end) = [];
					end
				end

				% Check maximum / minimum pressure
				if max(obj.profile(i).qc(j).pres) < obj.info(i).deep;
					% No deep data, flag
					flag_idx(j) = [92];
					continue
				elseif min(obj.profile(i).qc(j).pres) > obj.info(i).shallow;
					% No shallow data, flag
					flag_idx(j) = [93];
					continue
				end

				% Check minimum # of data depths (obj.info(i).depths)
				if length(unique(obj.profile(i).qc(j).pres)) < obj.info(i).depths
					% Not enough data, flag
					flag_idx(j) = [94];
					continue
				end
			end	

			% Log flagged profiles
			ID                   = [obj.profile(i).qc.ID];
			obj.flags(i).empty   = ID(find(flag_idx == 91));
			obj.flags(i).deep    = ID(find(flag_idx == 92));
			obj.flags(i).shallow = ID(find(flag_idx == 93));
			obj.flags(i).depths  = ID(find(flag_idx == 94));
			obj.flags(i).total   = ID(find(flag_idx >   0)); 
			disp(['Removing ',num2str(length(obj.flags(i).empty)),' profiles due to empty data'])
			disp(['Removing ',num2str(length(obj.flags(i).deep)),' profiles due to no deep data'])
			disp(['Removing ',num2str(length(obj.flags(i).shallow)),' profiles due to no shallow data'])
			disp(['Removing ',num2str(length(obj.flags(i).depths)),' profiles due small # of samples'])
			disp(['Removing ',num2str(length(obj.flags(i).total)),' total profiles due to the above QC'])

			%// Remove flagged profiles
			obj.profile(i).qc(flag_idx>0) = [];

			%// Save data for later
			data  = obj.profile(i).qc;
			flags = obj.flags(i);
			save([obj.info(i).qcdir,obj.info(i).qcfile],'data','flags','-v7.3')
		end
	end % end object method objQC

	function obj = objProc(obj)
	% ----------------------------------------------------------------
	% Process data to uniform pressure grid
	%
	% usage:
	%	obj = objProc(obj)
	%
	% flags applied:
	%	.bad_density = potential density is inverted (unstable) and 
	%                      shows an increase in density with depth that
	%		       is greater than settings.().density_thresh
	% -----------------------------------------------------------------

		% Announce routine
		disp(' ');
		disp('------------------------- ');
		disp('Interpolating data') 
		disp('------------------------- ');
		disp(' ');
		
		% Set up grids
		dbar2        = [0:2:2000]';         % 2-dbar grid
		dbar10       = [0:10:2000]';        % 10-dbar grid

		% Setup datatype loop
		for i = 1:length(obj.info)	
		
			% Check that qc data exists
			try
				obj.profile(i).qc(1).pres; % data is loaded
			catch
				tmpdata           = load([obj.info(i).qcdir,obj.info(i).qcfile]); % try loading data
				obj.profile(i).qc = tmpdata.data; 
				if i == 1
					obj.flags = tmpdata.flags;
				else
					obj.flags(i) = tmpdata.flags;
				end
				clear tmpdata
			end
	
			% Process settings
			ksrbandwidth = obj.info(i).ksrbandwidth; 
			
			% Initiate flag
			obj.flags(i).bad_density = []; % code == 91
			flag_idx = zeros(length(obj.profile(i).qc),1);
			
			% Copy data
			obj.profile(i).proc = obj.profile(i).qc;
			obj.profile(i).proc = rmfield(obj.profile(i).proc,{'temp_qc','salt_qc','pres_qc'});

			% Start profile loop
			disp('Start profile loop...')
			cnt   = 0; % iniitate disp_prog
			for j = 1:length(flag_idx)
			
                                % Announce loop progress
                                cnt = detectSCV.disp_prog(j,length(flag_idx),cnt);

				% Initiate sigma0, N2, spice
				obj.profile(i).proc(j).sigma0 = [];
				obj.profile(i).proc(j).N2     = [];
				obj.profile(i).proc(j).spice  = [];

			    	% Remove duplicate pressure entries (rare)
			   	[a,b]                       = unique(obj.profile(i).proc(j).pres);
			    	obj.profile(i).proc(j).pres = obj.profile(i).proc(j).pres(b);
			    	obj.profile(i).proc(j).temp = obj.profile(i).proc(j).temp(b);
			    	obj.profile(i).proc(j).salt = obj.profile(i).proc(j).salt(b);
			    	a = []; b = [];

			    	% Remove NaNs from P,T,S where they appear
				dat          = [obj.profile(i).proc(j).pres + ...
						obj.profile(i).proc(j).temp + ...
						obj.profile(i).proc(j).salt];
				obj.profile(i).proc(j).temp = obj.profile(i).proc(j).temp(~isnan(dat));
				obj.profile(i).proc(j).salt = obj.profile(i).proc(j).salt(~isnan(dat));
				obj.profile(i).proc(j).pres = obj.profile(i).proc(j).pres(~isnan(dat));
				dat = [];

			   	% Interpolate temporary salinity/temperature to 2dbar grid
			   	smooth_temp = obj.profile(i).proc(j).temp;
			    	smooth_salt = obj.profile(i).proc(j).salt;
			    	smooth_temp = interp1(obj.profile(i).proc(j).pres,smooth_temp,dbar2);
			    	smooth_salt = interp1(obj.profile(i).proc(j).pres,smooth_salt,dbar2);

			    	% Remove NaNs, apply smoothing filter
			    	smooth_pres = dbar2(~isnan(smooth_temp));
			    	smooth_temp = smooth_temp(~isnan(smooth_temp));
			    	smooth_salt = smooth_salt(~isnan(smooth_salt));
			    	rt          = detectSCV.ksr(smooth_pres,smooth_temp,ksrbandwidth,length(smooth_temp));
			    	rs          = detectSCV.ksr(smooth_pres,smooth_salt,ksrbandwidth,length(smooth_salt));
			    
			   	% Interpolate all data to 10dbar grid (including smoothed T/S)
			    	obj.profile(i).proc(j).temp = interp1(obj.profile(i).proc(j).pres,...
							      obj.profile(i).proc(j).temp,dbar10);
			    	obj.profile(i).proc(j).salt = interp1(obj.profile(i).proc(j).pres,...
							      obj.profile(i).proc(j).salt,dbar10);
			    	obj.profile(i).proc(j).pres = dbar10;
			    	smooth_temp  = interp1(rt.x,rt.f,dbar10);
			    	smooth_salt  = interp1(rs.x,rs.f,dbar10);
			    	rs 	     = [];
				rt 	     = [];

			    	% Add absolute salinity and conservative temperature (also do this for smoothed T/S)
			    	dat1            = [obj.profile(i).proc(j).pres + ...
						   obj.profile(i).proc(j).temp + ...
						   obj.profile(i).proc(j).salt];
			    	dat2            = [obj.profile(i).proc(j).pres + smooth_temp + smooth_salt];
			    	salt_abs        = gsw_SA_from_SP(obj.profile(i).proc(j).salt(~isnan(dat1)),...
						  obj.profile(i).proc(j).pres(~isnan(dat1)),...
						  obj.profile(i).meta(j).lon,obj.profile(i).meta(j).lat);
			    	smooth_salt_abs = gsw_SA_from_SP(smooth_salt(~isnan(dat2)),...
						  obj.profile(i).proc(j).pres(~isnan(dat2)),...
						  obj.profile(i).meta(j).lon,obj.profile(i).meta(j).lat);
			    	theta           = gsw_CT_from_t(salt_abs,obj.profile(i).proc(j).temp(~isnan(dat1)),...
						  obj.profile(i).proc(j).pres(~isnan(dat1)));
			    	smooth_theta    = gsw_CT_from_t(smooth_salt_abs,smooth_temp(~isnan(dat2)),...
						  obj.profile(i).proc(j).pres(~isnan(dat2)));

			    	% Add density and spiciness
			    	obj.profile(i).proc(j).sigma0 = gsw_sigma0(salt_abs,theta);
			    	obj.profile(i).proc(j).spice  = gsw_spiciness0(salt_abs,theta);

			    	% Get N2 from smoothed data
			    	[N2_mid,pres_mid] = gsw_Nsquared(smooth_salt_abs,smooth_theta,...
						    obj.profile(i).proc(j).pres(~isnan(dat2)),...
						    obj.profile(i).meta(j).lat);

			    	% Interpolate new variables to 10-dbar
			    	ndat           = isnan(N2_mid);
			    	obj.profile(i).proc(j).N2     = interp1(pres_mid(ndat==0),N2_mid(ndat==0),...
								obj.profile(i).proc(j).pres);
			    	obj.profile(i).proc(j).spice  = interp1(obj.profile(i).proc(j).pres(~isnan(dat1)),...
								obj.profile(i).proc(j).spice,obj.profile(i).proc(j).pres);
			    	obj.profile(i).proc(j).sigma0 = interp1(obj.profile(i).proc(j).pres(~isnan(dat1)),...
								obj.profile(i).proc(j).sigma0,obj.profile(i).proc(j).pres);
			    	smooth_temp  = []; smooth_salt     = [];
			 	smooth_theta = []; smooth_salt_abs = [];
				N2_mid       = []; pres_mid        = [];

			    	% Check that density is ascending with depth
			    	% For profiles that aren't, sort the data and see if there are very small differences.
			    	% If profiles show difference > obj.info(i).density_thresh, then reject them
				% as they are probably spikes in data. This small bit of code usually detectssmall
				% mixed-layer inversions from otherwise 'good' density profiles.
			    	if issorted(obj.profile(i).proc(j).sigma0(~isnan(obj.profile(i).proc(j).sigma0)))==0
					sigma0_orig = obj.profile(i).proc(j).sigma0(~isnan(obj.profile(i).proc(j).sigma0));
					sigma0_sort = sort(sigma0_orig);
					if max(abs(sigma0_orig-sigma0_sort))>obj.info(i).density_thresh;
				    		flag_idx(j) = 91;
					        continue
					else
				    		ind = find(isnan(obj.profile(i).proc(j).sigma0)==0);
				    		obj.profile(i).proc(j).sigma0(ind) = sort(obj.profile(i).proc(j).sigma0(ind));
					end
			    	end

			    	%// Only keep levels where all fields exist
			    	dat = [obj.profile(i).proc(j).sigma0 + ...
				       obj.profile(i).proc(j).temp   + ...
				       obj.profile(i).proc(j).salt   + ...
				       obj.profile(i).proc(j).spice  + ...
				       obj.profile(i).proc(j).N2];
			    	obj.profile(i).proc(j).temp(isnan(dat)==1)   = NaN;
			    	obj.profile(i).proc(j).salt(isnan(dat)==1)   = NaN;
			    	obj.profile(i).proc(j).pres(isnan(dat)==1)   = NaN;
			    	obj.profile(i).proc(j).spice(isnan(dat)==1)  = NaN;
			    	obj.profile(i).proc(j).N2(isnan(dat)==1)     = NaN;
			    	obj.profile(i).proc(j).sigma0(isnan(dat)==1) = NaN;
			end
			
			% Log flagged profiles
			ID               	 = [obj.profile(i).proc.ID];
			obj.flags(i).bad_density = ID(find(flag_idx == 91));
			obj.flags(i).total       = sort([obj.flags(i).total  obj.flags(i).bad_density],'ascend');
			disp(['Removing ',num2str(length(find(flag_idx == 91))),' profiles due to bad density'])
	
			% Remove bad data
			obj.profile(i).proc(find(flag_idx == 91)) = [];

			%// Save data for later
			data  = obj.profile(i).proc;
			flags = obj.flags(i);
			save([obj.info(i).procdir,obj.info(i).procfile],'data','flags','-v7.3')
		end
	end % end object method objProc

	function obj = objAnom(obj)
	% ----------------------------------------------------------------
	% Process climatology and get anomalies along isopyncals
	%
	% usage:
	%	obj = objAnom(obj)
	%
	% flags applied:
	%	.climatology = unable to assign climatological profile 
	%                      due to cast position (lon, lat) or due to a
	%		       shallow climatological profile based on
	%		       settings.().min_clim_pres
	% -----------------------------------------------------------------

		% Announce routine
		disp(' ');
		disp('------------------------- ');
		disp('Calculating anomalies')
		disp('------------------------- ');
		disp(' ');

		% Start datatype loop
		for i = 1:length(obj.info)
			
			% Check that proc data exists
			try
				obj.profile(i).proc(1).pres; % data is loaded
			catch
				tmpdata             = load([obj.info(i).procdir,obj.info(i).procfile]); % try loading data
				obj.profile(i).proc = tmpdata.data; 
				if i == 1
					obj.flags = tmpdata.flags;
				else
					obj.flags(i) = tmpdata.flags;
				end
				clear tmpdata
			end

			% Initiate flag
			obj.flags(i).climatology = []; % code == 91
			flag_idx                 = zeros(length(obj.profile(i).proc),1);
			
			% Load climatology
			load([obj.info(i).climdir,obj.info(i).climfile]);
	
			% Rename according to info.climvar
			eval(['climatology = ',obj.info(i).climvar,';'])
			
			% Clear original
			eval(['clear ',obj.info(i).climvar,';']);

			% Convert lon to 0 - 360
			climatology.lon(climatology.lon>360) = climatology.lon(climatology.lon>360)-360;

			% Get data lon and lat (fix lon);
			data.lon             = [obj.profile(i).proc.lon];
			data.lat             = [obj.profile(i).proc.lat];
			data.lon(data.lon<0) = data.lon(data.lon<0)+360;

			% Make matrix data (easier to vectorize)
			disp('Get matrix of CTD data...')
			data.temp   = [obj.profile(i).proc.temp];
			data.salt   = [obj.profile(i).proc.salt];
			data.pres   = [obj.profile(i).proc.pres]; % Same for all floats
			data.spice  = [obj.profile(i).proc.spice];
			data.N2     = [obj.profile(i).proc.N2];
			data.sigma0 = [obj.profile(i).proc.sigma0];
			
			% Start filling in climatological data for each float
			time_ind = nan(1,length(obj.profile(i).proc));
			lon_ind  = nan(1,length(obj.profile(i).proc));
			lat_ind  = nan(1,length(obj.profile(i).proc));

			% Find index of climatological pressure grid that exceeds info.min_clim_pres
			min_clim_pres_ind = find(climatology.pres >= obj.info(i).min_clim_pres);
			min_clim_pres_ind = min_clim_pres_ind(1);

			% Start profile loop
			disp('Start profile loop...')
			cnt   = 0; % initialize disp_prog
			for j = 1:length(obj.profile(i).proc)
                                
				% Announce loop progress
                                cnt = detectSCV.disp_prog(j,length(obj.profile(i).proc),cnt);

				% Grab time and nearest lon/lat on climatological grid
				time_ind(j) = obj.profile(i).proc(j).time(2);
			        lon_ind(j)  = min(find(abs(climatology.lon-data.lon(j)) == min(abs(climatology.lon-data.lon(j)))));
			        lat_ind(j)  = min(find(abs(climatology.lat-data.lat(j)) == min(abs(climatology.lat-data.lat(j)))));

				% Start climatological matrix
				obj.profile(i).clim(j).lon  = climatology.lon(lon_ind(j));
				obj.profile(i).clim(j).lat  = climatology.lat(lat_ind(j));
				obj.profile(i).clim(j).pres = climatology.pres';

    				% Check if temp/salinity data exists and at least info.min_clim_pres of climatology is available
				% ...if not get average from nearby grid cells
    				check = squeeze(climatology.temp(lon_ind(j),lat_ind(j),:,time_ind(j)));
    				if isnan(nanmean(check)) == 1 | length(check(~isnan(check))) < min_clim_pres_ind

					% Get climatological profiles of nearest 9 cells
					lon_rng = [lon_ind(j)-1:1:lon_ind(j)+1];
					lat_rng = [lat_ind(j)-1:1:lat_ind(j)+1];
					lon_rng(lon_rng > length(climatology.lon)) = [];
					lon_rng(lon_rng < 0) = [];
					lat_rng(lat_rng > length(climatology.lat)) = [];
					lat_rng(lat_rng < min(lat_ind)) = [];
					tmp.temp = climatology.temp(lon_rng,lat_rng,:,time_ind(j));
					tmp.salt = climatology.salt(lon_rng,lat_rng,:,time_ind(j));
					
					% Reshape 3D matrix into to 2D
					[a,b,c]  = size(tmp.temp);
					tmp.temp = reshape(tmp.temp,[a*b,c]);
					tmp.salt = reshape(tmp.salt,[a*b,c]);

					% Get average
					tmp.temp = nanmean(tmp.temp,1)';
					tmp.salt = nanmean(tmp.salt,1)';

					% Ignore attempts that still don't have at least 100m of climatological data 
					if length(tmp.temp(~isnan(tmp.temp))) < min_clim_pres_ind
	    					tmp.temp = nan(size(climatology.pres));
	    					tmp.salt = nan(size(climatology.pres));
						flag_idx(j) = 91;
					end
					
					% Save averaged profile
					obj.profile(i).clim(j).temp = tmp.temp;
					obj.profile(i).clim(j).salt = tmp.salt;
    				else
					% Data exists in region, grab T/S profile for that time of year
					obj.profile(i).clim(j).temp = squeeze(climatology.temp(lon_ind(j),lat_ind(j),:,time_ind(j)));
					obj.profile(i).clim(j).salt = squeeze(climatology.salt(lon_ind(j),lat_ind(j),:,time_ind(j)));
    				end

  				% Fix orientation (if necessary)
    				[a,b] = size(obj.profile(i).clim(j).temp);
    				if b > a
					obj.profile(i).clim(j).temp = obj.profile(i).clim(j).temp';
					obj.profile(i).clim(j).salt = obj.profile(i).clim(j).salt';
    				end
			end
	
			% Get matrix of climatolgoical CTD data for each float
			disp('Get matrix of climatological CTD data...')
			clim.lon  = [obj.profile(i).clim.lon];
			clim.lat  = [obj.profile(i).clim.lat];
			clim.temp = [obj.profile(i).clim.temp];
			clim.salt = [obj.profile(i).clim.salt];
			clim.pres = [obj.profile(i).clim.pres];

			% Switch back to -180:180 lon for calculations
 			clim.lon(clim.lon>360) = clim.lon(clim.lon>360)-360;

			% Add other fields
			disp('Grab salt_abs (patience!)...')
			clim.salt_abs = gsw_SA_from_SP(clim.salt,clim.pres,clim.lon,clim.lat);
			disp('Grab theta...')
			clim.theta = gsw_CT_from_t(clim.salt_abs,clim.temp,clim.pres);
			disp('Grab sigma0...')
			clim.sigma0 = gsw_sigma0(clim.salt_abs,clim.theta);
			disp('Grab spice...')
			clim.spice = gsw_spiciness0(clim.salt_abs,clim.theta);
			disp('Grab N2...')
			[clim.N2_mid,clim.pres_mid] = gsw_Nsquared(clim.salt_abs,clim.theta,clim.pres,clim.lat);
			dat = isnan(clim.N2_mid);

			% Get N2 on same pressure grid
			disp('Interpolate N2 back to regular pressure grid..')
			cnt   = 0; % initialize disp_prog
			for j = 1:length(obj.profile(i).proc)

                                % Announce loop progress
                                cnt = detectSCV.disp_prog(j,length(obj.profile(i).proc),cnt);
				
				% Try and interpolate
				try
					clim.N2(:,j) = interp1(clim.pres_mid(dat(:,j)==0,j),clim.N2_mid(dat(:,j)==0,j),clim.pres(:,j));
			    	catch
					clim.N2(:,j) = nan(length(climatology.pres),1);
			    	end
			end
			clear clim.N2_mid clim.pres_mid

			% Calculate anomalies
			disp('Calculate anomalies along isopycnals...')
			cnt   = 0; % initialize disp_prog
			for j = 1:length(obj.profile(i).proc)

                                % Announce loop progress
                                cnt = detectSCV.disp_prog(j,length(obj.profile(i).proc),cnt);

				% Copy meta fields
				obj.profile(i).anom(j).ID   = obj.profile(i).proc(j).ID;
				obj.profile(i).anom(j).lon  = obj.profile(i).proc(j).lon;
				obj.profile(i).anom(j).lat  = obj.profile(i).proc(j).lat;
				obj.profile(i).anom(j).time = obj.profile(i).proc(j).time;

				% Check for flags
				if flag_idx(j) == 91
					continue;
				end
			    
			        % Check for inverted climatological density (at surface)
			    	tmp_sigma0 = clim.sigma0(:,j);
			    	if issorted(tmp_sigma0(~isnan(tmp_sigma0)))==0
					ind                = find(isnan(clim.sigma0(:,j))==0);
					clim.sigma0(ind,j) = sort(clim.sigma0(ind,j));
			    	end

			        % Find where float and climatology potential density exists
			        data_ind = find(isnan(obj.profile(i).proc(j).sigma0)==0);
			   	clim_ind = [clim.sigma0(:,j) + clim.temp(:,j) + clim.salt(:,j) + clim.spice(:,j) + clim.N2(:,j)];

			    	% Interpolate climatology to float sigma0 grid
			   	clim_sigma0      = clim.sigma0(:,j);
			   	clim_temp        = clim.temp(:,j);
			    	clim_temp        = interp1(clim_sigma0(~isnan(clim_ind)),clim_temp(~isnan(clim_ind)),...
					           obj.profile(i).proc(j).sigma0(data_ind));
			    	filler           = nan(length(obj.profile(i).proc(j).sigma0),1);
			    	filler(data_ind) = clim_temp; clim_temp = filler;
			    	clim_salt        = clim.salt(:,j);
			    	clim_salt        = interp1(clim_sigma0(~isnan(clim_ind)),clim_salt(~isnan(clim_ind)),...
					           obj.profile(i).proc(j).sigma0(data_ind));
			    	filler           = nan(length(obj.profile(i).proc(j).sigma0),1);
			    	filler(data_ind) = clim_salt; clim_salt = filler;
			    	clim_spice       = clim.spice(:,j);
			    	clim_spice       = interp1(clim_sigma0(~isnan(clim_ind)),clim_spice(~isnan(clim_ind)),...
					           obj.profile(i).proc(j).sigma0(data_ind));
			    	filler           = nan(length(obj.profile(i).proc(j).sigma0),1);
			    	filler(data_ind) = clim_spice; clim_spice = filler;
			    	clim_N2          = clim.N2(:,j);
			    	clim_N2          = interp1(clim_sigma0(~isnan(clim_ind)),clim_N2(~isnan(clim_ind)),...
					           obj.profile(i).proc(j).sigma0(data_ind));
			    	filler           = nan(length(obj.profile(i).proc(j).sigma0),1);
			    	filler(data_ind) = clim_N2; clim_N2 = filler;
			    	clim_pres        = clim.pres(:,j);
			    	clim_pres        = interp1(clim_sigma0(~isnan(clim_ind)),clim_pres(~isnan(clim_ind)),...
					           obj.profile(i).proc(j).sigma0(data_ind));
			    	filler           = nan(length(obj.profile(i).proc(j).sigma0),1);
			    	filler(data_ind) = clim_pres; clim_pres = filler;

			        % Build anomaly structure
			        obj.profile(i).anom(j).temp_anom  = obj.profile(i).proc(j).temp - clim_temp;
			        obj.profile(i).anom(j).salt_anom  = obj.profile(i).proc(j).salt - clim_salt;
			        obj.profile(i).anom(j).spice_anom = obj.profile(i).proc(j).spice - clim_spice;
			        obj.profile(i).anom(j).N2_anom    = obj.profile(i).proc(j).N2 - clim_N2;
			        obj.profile(i).anom(j).pres_anom  = obj.profile(i).proc(j).pres - clim_pres;
			        obj.profile(i).anom(j).sigma0     = obj.profile(i).proc(j).sigma0;

				% Build final clim structure
				obj.profile(i).clim(j).ID     = obj.profile(i).proc(j).ID;
				obj.profile(i).clim(j).lon    = obj.profile(i).proc(j).lon;
				obj.profile(i).clim(j).lat    = obj.profile(i).proc(j).lat;
				obj.profile(i).clim(j).time   = obj.profile(i).proc(j).time;
				obj.profile(i).clim(j).temp   = clim_temp;
				obj.profile(i).clim(j).salt   = clim_salt;
				obj.profile(i).clim(j).pres   = clim_pres;
				obj.profile(i).clim(j).N2     = clim_N2;
				obj.profile(i).clim(j).spice  = clim_spice;
				obj.profile(i).clim(j).sigma0 = obj.profile(i).proc(j).sigma0;

				% Data check
				check = obj.profile(i).anom(j).temp_anom + ...
					obj.profile(i).anom(j).salt_anom + ...
					obj.profile(i).anom(j).pres_anom + ...
					obj.profile(i).anom(j).spice_anom + ...
					obj.profile(i).anom(j).N2_anom;
				if length(check(~isnan(check))) == 0
					flag_idx(j) = 91;
				end
			end

			% Log flagged profiles
			ID               	 = [obj.profile(i).anom.ID];
			obj.flags(i).climatology = ID(find(flag_idx == 91));
			obj.flags(i).total       = sort([obj.flags(i).total  obj.flags(i).climatology],'ascend');
			disp(['Removing ',num2str(length(find(flag_idx == 91))),' profiles due to no comparable climatology'])
			
			% Remove bad data
			obj.profile(i).anom(find(flag_idx == 91)) = [];
			obj.profile(i).clim(find(flag_idx == 91)) = [];
			
			% Save anomalies for later
			data  = obj.profile(i).anom;
			flags = obj.flags(i);
			save([obj.info(i).anomdir,obj.info(i).anomfile],'data','flags','-v7.3')

			% Save climatology for later
			data = obj.profile(i).clim;
			save([obj.info(i).anomdir,obj.info(i).climfile],'data','-v7.3');
		end
	end % end object method objAnom

	function obj = objNearby(obj)
	% ----------------------------------------------------------------
	% Find nearby profiles for calculating IQR
	%
	% usage:
	%	obj = objNearby(obj)
	% ----------------------------------------------------------------

		% Announce routine
		disp(' ');
		disp('------------------------- ');
		disp('Finding nearby floats')
		disp('------------------------- ');
		disp(' ');
		
		% Specify earth radius in meters
		earthRad = 6371000;

		% Initialize data structure
		tmp.lon  = [];
		tmp.lat  = [];
		tmp.ID   = [];
		tmp.mnth = [];
		tmp.day  = [];
		tmp.date = [];
		
		% Assemble lon, lat, time, and ID from every float
		for i = 1:length(obj.info)	
			% Check that anom data exists
			try
				obj.profile(i).anom(1).lon; % data is loaded
			catch
				tmpdata             = load([obj.info(i).anomdir,obj.info(i).anomfile]); % try loading data
				obj.profile(i).anom = tmpdata.data; 
				if i == 1
					obj.flags = tmpdata.flags;
				else
					obj.flags(i) = tmpdata.flags;
				end
				clear tmpdata
			end

			% Gather data, add to existing array if i > 1
			lon      = [obj.profile.anom.lon];
			tmp.lon  = [tmp.lon lon];
			lat      = [obj.profile.anom.lat];
			tmp.lat  = [tmp.lat lat];
			ID       = [obj.profile.anom.ID];
			tmp.ID   = [tmp.ID ID];
			mnth     = nan(size(ID));
			day      = nan(size(ID));
			date     = nan(size(ID));
			for i = 1:length(tmp.lon)
	    			mnth(1,i) = obj.profile.anom(i).time(2);
    				day(1,i)  = obj.profile.anom(i).time(3);
    				date(1,i) = datenum(obj.profile.anom(i).time);
			end
			tmp.mnth = [tmp.mnth mnth];
			tmp.day  = [tmp.day day];
			tmp.date = [tmp.date date];
		end 	

		% Now compare every float with its neighbors in space / time
		disp('Start profile loop...')
		cnt   = 0; % initiate disp_prog
		for i = 1:length(tmp.lon)

                        % Announce loop progress
                        cnt = detectSCV.disp_prog(i,length(tmp.lon),cnt);
    			
			% Find distance
			lon   = tmp.lon(i)*ones(size(tmp.lon));
			lat   = tmp.lat(i)*ones(size(tmp.lat));
			dist  = distance(lat,lon,tmp.lat,tmp.lon,earthRad);
			dist  = dist / 1000; % m to km

    			% Re-do if longitude is near zero
			if lon <= 2 | lon >= 358
				lon           = lon + 360;
				alon          = tmp.lon;
				alon(alon<=2) = alon(alon<=2)+360;
				dist = distance(lat,lon,tmp.lat,alon,earthRad);
				dist = dist / 1000; % m to km
			end

		    	% Get distance index based on settings
		    	didx = find(dist <= obj.info(1).search_dist);

		    	% Setup date_span to find when dates are in range (based on settings)
			month     = tmp.mnth(i);
		    	day       = tmp.day(i);
		    	yr_span   = [1995:2025]'; %'
		    	mn_span   = month*ones(size(yr_span));
		    	dy_span   = day*ones(size(yr_span));
		    	date_span = datenum([yr_span mn_span dy_span]); 
		    	date_span = [date_span-obj.info(1).search_time date_span+obj.info(1).search_time];

		    	% Get time index
			tidx = [];
		    	for ii = 1:length(date_span)
				di   = find(date_span(ii,1) < tmp.date & tmp.date < date_span(ii,2));
				tidx = [tidx di]; 
		    	end

			% Find where the intersect, save as nearby floats
		    	if isempty(didx) == 0 & isempty(tidx) == 0
				aidx = intersect(didx,tidx);
				if length(obj.info) == 1
					obj.profile.nearby(i).IDs = tmp.ID(aidx);
				else
					if tmp.ID(i) < 2e8 % argo ID
						obj.profile(1).nearby(i).IDs = tmp.ID(aidx);
					elseif tmp.ID(i) >= 2e8 % meop ID
						j = 1-length(obj.profile(1).nearby);
						obj.profile(2).nearby(j).IDs = tmp.ID(aidx);
					end
				end
		    	else
				if tmp.ID(i) < 2e8 % argo ID
					obj.profile(1).nearby(i).IDs = NaN;
				elseif tmp.ID(i) >= 2e8 % meop ID
					j = 1-length(obj.profile(1).nearby);
					obj.profile(2).nearby(j).IDs = NaN;
				end
		    	end
		end

		% Save results for later
		for i = 1:length(obj.info)
			% Save nearby IDs for later
			data  = obj.profile(i).nearby;
			save([obj.info(i).nearbydir,obj.info(i).nearbyfile],'data','-v7.3')
		end
	end % end object method objNearby

	function obj = objIQR(obj)
	% ----------------------------------------------------------------
	% Calculate IQR of anomalies along isopyncals
	%
	% usage:
	%	obj = objIQR(obj)
	%
	% flags applied:
	%	.min_nearby = Not enough nearby/neartime floats to 
	%		      calculate IQR, based on settings.().min_nearby 
	% -----------------------------------------------------------------
		
		% Announce routine
		disp(' ');
		disp('------------------------- ');
		disp('Finding IQR + thresholds')
		disp('------------------------- ');
		disp(' ');

		% Check that data is loaded
		for i = 1:length(obj.info)

			% Check that anom data, nearby data exists
			try
				obj.profile(i).anom(1).pres;  % anom data is loaded
				obj.profile(i).nearby(1).IDs; % nearby data is loaded
			catch
				tmpdata             = load([obj.info(i).anomdir,obj.info(i).anomfile]); % try loading data
				obj.profile(i).anom = tmpdata.data; 
				tmpdata             = load([obj.info(i).nearbydir,obj.info(i).nearbyfile]); % try loading data
				obj.profile(i).anom = tmpdata.data; 
				if i == 1
					obj.flags = tmpdata.flags;
				else
					obj.flags(i) = tmpdata.flags;
				end
				clear tmpdata
			end
		end

		% Define variables for IQR calculation
		vars      = {'spice_anom','N2_anom'};
		vars_fill = {'spice','N2'};
 
		% Get matrix of potential densities from every float
		if length(obj.info) > 1 % both argo/meop
			sigma0_1 = [];
			sigma0_2 = [];
			sigma0   = [];
			sigma0_1 = [obj.profile(1).anom.sigma0];
			sigma0_2 = [obj.profile(2).anom.sigma0];
			sigma0   = [sigma0_1 sigma1_2];
			clear sigma0_1 sigma0_2
		else
			sigma0 = [];
			sigma0 = [obj.profile.anom.sigma0];
		end

		% Get cell array of nearby float IDs
		if length(obj.info) > 1 % both argo/meop
			ID_1 = [];
			ID_2 = [];
			ID   = [];	
			ID_1 = [obj.profile(1).anom.ID];
			ID_2 = [obj.profile(2).anom.ID];
			ID   = [ID_1 ID_2];
			% Nearby floats
			nearby_1 = [];
			nearby_2 = [];
			nearby   = [];
			for i = 1:length(obj.profile(1).anom)
				nearby_1{i} = obj.profile(1).nearby(i).IDs;
			end
			for i = 1:length(obj.profile(2).anom)
				nearby_2{i} = obj.profile(2).nearby(i).IDs;
			end
			nearby = [nearby_1 nearby_2];
			clear nearby_1 nearby_2
		else
			ID     = [];
			nearby = [];
			ID = [obj.profile.anom.ID];
			for i = 1:length(obj.profile.anom)
				nearby{i} = obj.profile.nearby(i).IDs;
			end
		end
		
		% Initiate flag index
		flag_idx = zeros(1,length(nearby));

		% Loop through vars and calculate IQR + thresholds for each float
		for i = 1:length(vars)
			
			% Grab matrix of data
			if length(obj.info) > 1 % both argo/meop
				tmpdata_1 = [obj.profile(1).anom.(vars{i})];
				tmpdata_2 = [obj.profile(2).anom.(vars{i})];
				tmpdata   = [tmpdata_1 tmpdata_2];
				clear tmpdata_1 tmpdata_2
			else
				tmpdata = [obj.profile.anom.(vars{i})];
			end
	
			% Initiate matrices
			tmpdata_IQR    = nan(size(tmpdata));
			tmpdata_p25    = nan(size(tmpdata));
			tmpdata_p75    = nan(size(tmpdata));
			tmpdata_lim_hi = nan(size(tmpdata));
			tmpdata_lim_lo = nan(size(tmpdata));
			
			% Start profile loop
			disp('Start profile loop...')
			cnt   = 0; % initiate disp_prog
			for j = 1:length(tmpdata)
                        	
				% Announce loop progress
                        	cnt = detectSCV.disp_prog(j,tmpdata,cnt);

				% Check for minimum number of nearby profiles
				if length(nearby{j}) < obj.info(1).min_nearby
					flag_idx(j) = 91;
					continue
				end

				% Grab indices
				ID_idx = find(ismember(ID,nearby{j})==1);

				% Grab nearby data and densities
				nearbydata = [];
				nearbydens = [];
				nearbydata = tmpdata(:,ID_idx);
				nearbydens = sigma0(:,ID_idx);
				
				% Find where data exists for sample cast
				data_idx = find(isnan(sigma0(:,j))==0);

				% Initiate data_grid matrix
				data_grid = nan(size(tmpdata,1),length(ID_idx));

				% Interpolate each nearby cast to castdens levels
				for k = 1:length(ID_idx)
				
					% Grab cast data, find where it exists
					castdata = []; 
					castdens = [];
					castdata = nearbydata(:,k);
					castdens = nearbydens(:,k);
					dat      = [castdata + castdens];
					castdata = castdata(~isnan(dat),:);
					castdens = castdens(~isnan(dat),:);
				
					% Only allow unique densities
					[~,ia,~] = unique(castdens);
					castdata = castdata(ia,:);
					castdens = castdens(ia,:);

					% Interpolate to sample cast densities
					filler           = nan(length(sigma0(:,j)),1);
					if isempty(castdata) | length(castdata)==1
						continue
					else
						data_int         = interp1(castdens,castdata,sigma0(data_idx,j));
						filler(data_idx) = data_int;
					end
					data_grid(:,k)   = filler;
				end
				
				% Grab interquartile range and percentiles
				[a,b] = size(data_grid);
				for k = 1:b
					ind = find(isnan(data_grid(:,k))==1);
					if length(ind) < obj.info(1).min_nearby							
						data_grid(:,k) == NaN;
					end
				end
				variqr = [vars{i},'_iqr'];
				varp25 = [vars{i},'_p25'];
				varp75 = [vars{i},'_p75'];
				varlim = [vars{i},'_limits'];

				% Set IQR, percentiles, and IQR threshold (p25 - iqr_mult*IQR, p75 + iqr_mult*IQR)
				tmpdata_iqr(:,j)    = iqr(data_grid');
				tmpdata_p25(:,j)    = prctile(data_grid',25);
				tmpdata_p75(:,j)    = prctile(data_grid',75);
				if length(obj.info) > 1 % both argo/meop
					if ID(j) < 2e8 % argo ID
						tmpdata_lim_hi(:,j) = tmpdata_p75(:,j) + obj.info(1).iqr_mult.*tmpdata_iqr(:,j);
						tmpdata_lim_lo(:,j) = tmpdata_p25(:,j) - obj.info(1).iqr_mult.*tmpdata_iqr(:,j);
					elseif ID(j) >= 2e8 % meop ID
						tmpdata_lim_hi(:,j) = tmpdata_p75(:,j) + obj.info(2).iqr_mult.*tmpdata_iqr(:,j);
						tmpdata_lim_lo(:,j) = tmpdata_p25(:,j) - obj.info(2).iqr_mult.*tmpdata_iqr(:,j);
					end
				else
					tmpdata_lim_hi(:,j) = tmpdata_p75(:,j) + obj.info.iqr_mult.*tmpdata_iqr(:,j);
					tmpdata_lim_lo(:,j) = tmpdata_p25(:,j) - obj.info.iqr_mult.*tmpdata_iqr(:,j);
				end	
			end
	
			% Assign tmpdata to profile
			for j = 1:length(tmpdata)
				if length(obj.info) > 1 % both argo/meop
					if ID(j) < 2e8 % argo ID
						obj.profile(1).iqr(j).ID       = obj.profile(1).anom(j).ID;
						obj.profile(1).iqr(j).lon      = obj.profile(1).anom(j).lon;
						obj.profile(1).iqr(j).lat      = obj.profile(1).anom(j).lat;
						obj.profile(1).iqr(j).time     = obj.profile(1).anom(j).time;
						obj.profile(1).iqr(j).(variqr) = tmpdata_iqr(:,j);
						obj.profile(1).iqr(j).(varp25) = tmpdata_p25(:,j);
						obj.profile(1).iqr(j).(varp75) = tmpdata_p75(:,j);
						obj.profile(1).iqr(j).(varlim) = [tmpdata_lim_lo(:,j) tmpdata_lim_hi(:,j)];
					elseif ID(j) >= 2e8 % meop ID
						k = j - length(obj.profile(1).iqr);
						obj.profile(2).iqr(k).ID       = obj.profile(2).anom(k).ID;
						obj.profile(2).iqr(k).lon      = obj.profile(2).anom(k).lon;
						obj.profile(2).iqr(k).lat      = obj.profile(2).anom(k).lat;
						obj.profile(1).iqr(k).time     = obj.profile(2).anom(k).time;
						obj.profile(2).iqr(k).(variqr) = tmpdata_iqr(:,j);
						obj.profile(2).iqr(k).(varp25) = tmpdata_p25(:,j);
						obj.profile(2).iqr(k).(varp75) = tmpdata_p75(:,j);
						obj.profile(2).iqr(k).(varlim) = [tmpdata_lim_lo(:,j) tmpdata_lim_hi(:,j)];
					end
				else
					obj.profile.iqr(j).ID       = obj.profile.anom(j).ID;
					obj.profile.iqr(j).lon      = obj.profile.anom(j).lon;
					obj.profile.iqr(j).lat      = obj.profile.anom(j).lat;
					obj.profile.iqr(j).time     = obj.profile.anom(j).time;
					obj.profile.iqr(j).(variqr) = tmpdata_iqr(:,j);
					obj.profile.iqr(j).(varp25) = tmpdata_p25(:,j);
					obj.profile.iqr(j).(varp75) = tmpdata_p75(:,j);
					obj.profile.iqr(j).(varlim) = [tmpdata_lim_lo(:,j) tmpdata_lim_hi(:,j)];
				end
			end
		end

		% Save data
		if length(obj.info) > 1 % both argo/meop
			flag_1 = flag_idx(1:length(obj.profile(1).anom));
			flag_2 = flag_idx(length(obj.profile(1).anom)+1:end);

			% Log flagged profiles
			obj.flags(1).min_nearby = ID_1(find(flag_1==91));
			obj.flags(1).total      = sort([obj.flags(1).total obj.flags(1).min_nearby],'ascend') 
			obj.flags(2).min_nearby = ID_2(find(flag_2==91));
			obj.flags(2).total      = sort([obj.flags(2).total obj.flags(2).min_nearby],'ascend') 
			
			% Remove bad data
			obj.profile(1).iqr(find(flag_1==91)) = [];
			obj.profile(2).iqr(find(flag_2==91)) = [];

			% Save data for later
			data  = obj.profile(1).iqr;
			flags = obj.flags(1);
			save([obj.info(1).iqrdir,obj.info(2).iqrfile],'data','flags','-v7.3')
			data  = obj.profile(2).iqr;
			flags = obj.flags(2);
			save([obj.info(2).iqrdir,obj.info(2).iqrfile],'data','flags','-v7.3')
		else
			% Add flags
			obj.flags.min_nearby = ID(find(flag_idx==91));
			obj.flags.total      = sort([obj.flags.total obj.flags.min_nearby],'ascend');
			
			% Remove bad data
			obj.profile.iqr(find(flag_idx==91)) = [];

			% Save data for later
			data  = obj.profile.iqr;
			flags = obj.flags;
			save([obj.info.iqrdir,obj.info.iqrfile],'data','flags','-v7.3')
		end
	end % end object method objIQR

	function obj = objDetect(obj)
	% ----------------------------------------------------------------
	% Detect isopycnals with anomalous spice/N2 anomalies
	%
	% usage:
	%	obj = objDetect(obj)
	% -----------------------------------------------------------------

		% Announce routine
		disp(' ');
		disp('------------------------- ');
		disp('Detecting spice + N2 outliers')
		disp('------------------------- ');
		disp(' ');

		% Load iqr and anomaly data
		for i = 1:length(obj.info)	
			
			% Check that anom data exists
			try
				obj.profile(i).anom(1).lon;           % data is loaded
				obj.profile(i).iqr(1).spice_anom_iqr; % data is loaded 
			catch
				tmpdata             = load([obj.info(i).anomdir,obj.info(i).anomfile]); % load anom data
				obj.profile(i).anom = tmpdata.data; 
				tmpdata             = load([obj.info(i).iqrdir,obj.info(i).iqrfile]);   % load iqr data
				obj.profile(i).iqr  = tmpdata.data;
				if i == 1
					obj.flags = tmpdata.flags;
				else
					obj.flags(i) = tmpdata.flags;
				end
				clear tmpdata
			end
			
			% Remove flagged data
			IDa = [obj.profile(i).anom.ID];
			IDi = [obj.profile(i).iqr.ID]; 
			IDf = [obj.flags(i).total];
			inda = find(ismember(IDa,IDf));
			indi = find(ismember(IDi,IDf));
			obj.profile(i).anom(inda) = [];
			obj.profile(i).iqr(inda)  = []; 
		end	
	end % end object method objDetect
end % end define object methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static) % start define Static methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	function cnt = disp_prog(i,li,cnt)
	% ------------------------------------------------------------
	% disp_prog announces progress in a for-loop
	% i is the current iteration
	% li is the number of total iterations
	% cnt is the current progress
	% ------------------------------------------------------------
		progress = 0:floor((li/10)):li;
		if ismember(i,progress)
			cnt = cnt + 10;
			disp([num2str(cnt),'% complete']);
		end
	end % end static method disp_prog

	function params=parse_pv_pairs(params,pv_pairs)
	% ------------------------------------------------------------
	% parse_pv_pairs: parses sets of property value pairs, allows defaults
	% usage: params=parse_pv_pairs(default_params,pv_pairs)
	%
	% arguments: (input)
	%  default_params - structure, with one field for every potential
	%             property/value pair. Each field will contain the default
	%             value for that property. If no default is supplied for a
	%             given property, then that field must be empty.
	%
	%  pv_array - cell array of property/value pairs.
	%             Case is ignored when comparing properties to the list
	%             of field names. Also, any unambiguous shortening of a
	%             field/property name is allowed.
	%
	% arguments: (output)
	%  params   - parameter struct that reflects any updated property/value
	%             pairs in the pv_array.
	%
	% Example usage:
	% First, set default values for the parameters. Assume we
	% have four parameters that we wish to use optionally in
	% the function examplefun.
	%
	%  - 'viscosity', which will have a default value of 1
	%  - 'volume', which will default to 1
	%  - 'pie' - which will have default value 3.141592653589793
	%  - 'description' - a text field, left empty by default
	%
	% The first argument to examplefun is one which will always be
	% supplied.
	%
	%   function examplefun(dummyarg1,varargin)
	%   params.Viscosity = 1;
	%   params.Volume = 1;
	%   params.Pie = 3.141592653589793
	%
	%   params.Description = '';
	%   params=parse_pv_pairs(params,varargin);
	%   params
	%
	% Use examplefun, overriding the defaults for 'pie', 'viscosity'
	% and 'description'. The 'volume' parameter is left at its default.
	%
	%   examplefun(rand(10),'vis',10,'pie',3,'Description','Hello world')
	%
	% params = 
	%     Viscosity: 10
	%        Volume: 1
	%           Pie: 3
	%   Description: 'Hello world'
	%
	% Note that capitalization was ignored, and the property 'viscosity'
	% was truncated as supplied. Also note that the order the pairs were
	% supplied was arbitrary.
	% ------------------------------------------------------------

		npv = length(pv_pairs);
		n = npv/2;

		if n~=floor(n)
		  error 'Property/value pairs must come in PAIRS.'
		end
		if n<=0
		  % just return the defaults
		  return
		end

		if ~isstruct(params)
		  error 'No structure for defaults was supplied'
		end

		% there was at least one pv pair. process any supplied
		propnames = fieldnames(params);
		lpropnames = lower(propnames);
		for i=1:n
		  p_i = lower(pv_pairs{2*i-1});
		  v_i = pv_pairs{2*i};
		  
		  ind = strmatch(p_i,lpropnames,'exact');
		  if isempty(ind)
		    ind = find(strncmp(p_i,lpropnames,length(p_i)));
		    if isempty(ind)
		      error(['No matching property found for: ',pv_pairs{2*i-1}])
		    elseif length(ind)>1
		      error(['Ambiguous property name: ',pv_pairs{2*i-1}])
		    end
		  end
		  p_i = propnames{ind};
		  
		  % override the corresponding default in params
		  params = setfield(params,p_i,v_i);
		  
		end
	end % end static method parse_pv_pairs.m

	function r=ksr(x,y,h,N)
	% ------------------------------------------------------------
	% KSR   Kernel smoothing regression
	%
	% r=ksr(x,y) returns the Gaussian kernel regression in structure r such that
	%   r.f(r.x) = y(x) + e
	% The bandwidth and number of samples are also stored in r.h and r.n
	% respectively.
	%
	% r=ksr(x,y,h) performs the regression using the specified bandwidth, h.
	%
	% r=ksr(x,y,h,n) calculates the regression in n points (default n=100).
	%
	% Without output, ksr(x,y) or ksr(x,y,h) will display the regression plot.
	%
	% Algorithm
	% The kernel regression is a non-parametric approach to estimate the
	% conditional expectation of a random variable:
	%
	% E(Y|X) = f(X)
	%
	% where f is a non-parametric function. Based on the kernel density
	% estimation, this code implements the Nadaraya-Watson kernel regression
	% using the Gaussian kernel as follows:
	%
	% f(x) = sum(kerf((x-X)/h).*Y)/sum(kerf((x-X)/h))
	%
	% See also gkde, ksdensity

	% Example 1: smooth curve with noise
	%{
	x = 1:100;
	y = sin(x/10)+(x/50).^2;
	yn = y + 0.2*randn(1,100);
	r=ksr(x,yn);
	plot(x,y,'b-',x,yn,'co',r.x,r.f,'r--','linewidth',2)
	legend('true','data','regression','location','northwest');
	title('Gaussian kernel regression')
	%}
	% Example 2: with missing data
	%{
	x = sort(rand(1,100)*99)+1;
	y = sin(x/10)+(x/50).^2;
	y(round(rand(1,20)*100)) = NaN;
	yn = y + 0.2*randn(1,100);
	r=ksr(x,yn);
	plot(x,y,'b-',x,yn,'co',r.x,r.f,'r--','linewidth',2)
	legend('true','data','regression','location','northwest');
	title('Gaussian kernel regression with 20% missing data')
	%}
	% By Yi Cao at Cranfield University on 12 March 2008.
	% ------------------------------------------------------------
		
		% Check input and output
		error(nargchk(2,4,nargin));
		error(nargoutchk(0,1,nargout));
		if numel(x)~=numel(y)
		    error('x and y are in different sizes.');
		end

		x=x(:);
		y=y(:);
		% clean missing or invalid data points
		inv=(x~=x)|(y~=y);
		x(inv)=[];
		y(inv)=[];

		% Default parameters
		if nargin<4
		    N=100;
		elseif ~isscalar(N)
		    error('N must be a scalar.')
		end
		r.n=length(x);
		if nargin<3
		    % optimal bandwidth suggested by Bowman and Azzalini (1997) p.31
		    hx=median(abs(x-median(x)))/0.6745*(4/3/r.n)^0.2;
		    hy=median(abs(y-median(y)))/0.6745*(4/3/r.n)^0.2;
		    h=sqrt(hy*hx);
		    if h<sqrt(eps)*N
			error('There is no enough variation in the data. Regression is meaningless.')
		    end
		elseif ~isscalar(h)
		    error('h must be a scalar.')
		end
		r.h=h;

		% Gaussian kernel function
		kerf=@(z)exp(-z.*z/2)/sqrt(2*pi);

		r.x=linspace(min(x),max(x),N);
		r.f=zeros(1,N);
		for k=1:N
		    z=kerf((r.x(k)-x)/h);
		    r.f(k)=sum(z.*y)/sum(z);
		end

		% Plot
		if ~nargout
		    plot(r.x,r.f,'r',x,y,'bo')
		    ylabel('f(x)')
		    xlabel('x')
		    title('Kernel Smoothing Regression');
		end


	end % end static method ksr.m

	function woa18_format(fdir,fname)
	% ------------------------------------------------------------
	% Format WOA18 for SCV detection
	% ------------------------------------------------------------
		
		% Load WOA18 data
		load([fdir,fname])
		
		% Save in easier-to-use structure
		woa18.lon    = woa_lons;
		woa18.lat    = woa_lats;
		woa18.month  = woa_months;
		woa18.temp   = woa_temp;
		woa18.salt   = woa_salt;
		woa18.pres   = woa_depths; 

		% Save new mat file
		save([fdir,'woa18_format.mat'],'woa18','-v7.3');
	end % end static method woa18_format

	function meop_format(fdir,fname)
	% ------------------------------------------------------------
	% Format MEOP CTD data for SCV detection
	% ------------------------------------------------------------
		
		% Load SealData.mat
		load([fdir,fname])
		
		% Save as meop
		meop = sealdata_all;

		% Save new mat file
		save([fdir,'meop_raw_data.mat'],'meop','-v7.3')
	end

	function getSettings
	% ----------------------------------------------------------------------
	% Script to update settings.mat
	% ----------------------------------------------------------------------
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Raw data directories and filenames
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Used in objInit, objQC
		% Argo
		settings.argo.rawdir  = '/data/project1/demccoy/data/argo/update/';
		settings.argo.rawfile = 'RAW_global_argo.mat';
		% Meop
		settings.meop.rawdir  = '/data/project1/demccoy/data/meop/';
		settings.meop.rawfile = 'meop_raw_data.mat';

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Quality controlled data directories and filenames
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Used in objInit, objQC
		% Argo
		settings.argo.qcdir  = '/data/project1/demccoy/data/argo/update/';
		settings.argo.qcfile = 'argo_proc_data.mat';
		% Meop
		settings.meop.qcdir  = '/data/project1/demccoy/data/meop/';
		settings.meop.qcfile = 'argo_meop_data.mat';

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Processed data directories
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Used in objInit, objProc
		% Argo
		settings.argo.procdir  = '/data/project1/demccoy/data/argo/update/';
		settings.argo.procfile = 'argo_proc_data.mat';
		% Meop
		settings.meop.procdir  = '/data/project1/demccoy/data/meop/';
		settings.meop.procfile = 'argo_meop_data.mat';

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Anomalies data directories
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Used in objInit, objAnom
		% Argo
		settings.argo.anomdir  = '/data/project1/demccoy/data/argo/update/';
		settings.argo.anomfile = 'argo_anom_data.mat';
		settings.argo.anomclim = 'argo_clim_data.mat';
		% Meop
		settings.meop.anomdir  = '/data/project1/demccoy/data/meop/';
		settings.meop.anomfile = 'meop_anom_data.mat';
		settings.meop.anomclim = 'meop_clim_data.mat';

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Nearby ID data directories
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Used in objInit, objNearby, objIQR
		% Argo
		settings.argo.neardir  = '/data/project1/demccoy/data/argo/update/';
		settings.argo.nearfile = 'argo_nearby_data.mat';
		% Meop
		settings.meop.neardir  = '/data/project1/demccoy/data/meop';
		settings.meop.nearfile = 'meop_nearby_data.mat';

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% IQR data directories
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Used in objInit, objIQR
		% Argo
		settings.argo.iqrdir  = '/data/project1/demccoy/data/argo/update/';
		settings.argo.iqrfile = 'argo_iqr_data.mat';
		% Meop
		settings.meop.iqrdir  = '/data/project1/demccoy/data/meop/';
		settings.meop.iqrfile = 'meop_iqr_data.mat';

		%%%%%%%%%%%%%%%
		% GSW directory
		%%%%%%%%%%%%%%%
		% Used in objInit, objProc, objAnom
		settings.gswdir = '/data/project1/matlabpathfiles/gsw_matlab_v3_06_11/'; 

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Climatology directory and filename
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Used in objInit, objAnom
		% Argo
		settings.argo.climdir  = '/data/project1/demccoy/Argo_SCVs/update/scripps_clim/';
		settings.argo.climfile = 'RG_ArgoClim_TS_Climatologies.mat';
		settings.argo.climvar  = 'RG_ArgoClim_TS_Climatologies';
		% Meop
		settings.meop.climdir  = '/data/project1/demccoy/data/meop/';
		settings.meop.climfile = 'woa18_format.mat';
		settings.meop.climvar  = 'woa18';

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Quality control thresholds
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Used in objInit, objQC
		% Argo
		settings.argo.qcmatr  = [100  1000  65; 
					 1000  inf 105];
		settings.argo.deep    = 200; 
		settings.argo.shallow = 50; 
		settings.argo.depths  = 10; 
		% Meop
		settings.meop.qcmatr  = settings.argo.qcmatr;  % Copy argo settings
		settings.meop.deep    = settings.argo.deep;    % Copy argo settings
		settings.meop.shallow = settings.argo.shallow; % Copy argo settings
		settings.meop.depths  = settings.argo.depths;  % Copy argo settings

		%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Interpolation thresholds
		%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Used in objInit, objProc
		% Argo
		settings.argo.ksrbandwidth   = 10; 
		settings.argo.density_thresh = 0.1;
		% Meop
		settings.meop.ksrbandwidth   = settings.argo.ksrbandwidth;   % Copy argo settings
		settings.meop.density_thresh = settings.argo.density_thresh; % Copy argo settings

		%%%%%%%%%%%%%%%%%%%%%%
		% Climatology settings
		%%%%%%%%%%%%%%%%%%%%%%
		% Used in objAnom
		% Argo
		settings.argo.min_clim_pres = 100;
		% Meop
		settings.meop.min_clim_pres = settings.argo.min_clim_pres; % Copy argo settings

		%%%%%%%%%%%%%%%%%%%%%%%
		% Nearby float settings
		%%%%%%%%%%%%%%%%%%%%%%%
		% Use in objInit, objNearby, objIQR
		% Argo
		settings.argo.search_dist   = 200;
		settings.argo.search_time   = 15;
		settings.argo.min_nearby    = 30;
		% Meop
		settings.meop.search_dist   = settings.argo.search_dist;   % Copy argo settings
		settings.meop.search_time   = settings.argo.search_time;   % Copy argo settings
		settings.meop.min_nearby    = settings.argo.min_nearby;    % Copy argo settings

		%%%%%%%%%%%%%%%%%%%%%%%
		% IQR settings
		%%%%%%%%%%%%%%%%%%%%%%%
		% Argo
		settings.argo.iqr_mult = 1.5;		
		% Meop
		settings.meop.iqr_mult = settings.argo.iqr_mult;           % Copy argo settings

		% Save
		save('settings.mat','settings');
	end % end method getSettings
end % end define Static methods
end % end define detectSCV object

