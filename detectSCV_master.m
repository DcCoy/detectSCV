% Master script to run detectSCV in the correct order

% Reload settings.mat
detectSCV.getSettings;

% Choose which dataset to process
choice = 'meop';

% Initialize detectSCV object
obj = detectSCV();

% Load settings
obj = objInit(obj,choice);

% Quality control 
obj = objQC(obj);

% Process data
obj = objProc(obj);

% Find anomalies
obj = objAnom(obj);

% Find nearby floats
obj = objNearby(obj);

% Calculate IQR
obj = objIQR(obj);

