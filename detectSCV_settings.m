% Set data paths, thresholds for quality control,
% and other processing settings here.
%
% See detectSCV_readme.txt

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
% Meop
settings.meop.anomdir  = '/data/project1/demccoy/data/meop/';
settings.meop.anomfile = 'meop_anom_data.mat';

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
settings.meop.qcmatr  = settings.argo.qcmatr;  % Copy settings
settings.meop.deep    = settings.argo.deep;    % Copy settings
settings.meop.shallow = settings.argo.shallow; % Copy settings
settings.meop.depths  = settings.argo.depths;  % Copy settings

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolation thresholds
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Used in objInit, objProc
% Argo
settings.argo.ksrbandwidth   = 10; 
settings.argo.density_thresh = 0.1;
% Meop
settings.meop.ksrbandwidth   = settings.argo.ksrbandwidth;   % Copy settings
settings.meop.density_thresh = settings.argo.density_thresh; % Copy settings

%%%%%%%%%%%%%%%%%%%%%%
% Climatology settings
%%%%%%%%%%%%%%%%%%%%%%
% Used in objAnom
% Argo
settings.argo.min_clim_pres = 100;
% Meop
settings.meop.min_clim_pres = settings.argo.min_clim_pres; % Copy settings

%%%%%%%%%%%%%%%%%%%%%%%
% Nearby float settings
%%%%%%%%%%%%%%%%%%%%%%%
% Use in objInit, objNearby, objIQR
% Argo
settings.argo.search_dist   = 200;
settings.argo.search_time   = 15;
settings.argo.min_nearby    = 30;

% Meop
settings.meop.search_dist   = settings.argo.search_dist;   % Copy settings
settings.meop.search_time   = settings.argo.search_time;   % Copy settings
settings.meop.min_nearby    = settings.argo.min_nearby;    % Copy settings

% Save
save('settings.mat','settings');
