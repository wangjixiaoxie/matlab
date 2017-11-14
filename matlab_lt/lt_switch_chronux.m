function lt_switch_chronux(version)
%% lt 11/9/17 - add or remove chronux package from path

% inputs
% version = 0; no chronux
% = 1; adds chronux
%% Add + all subfolders

if version==1;
    % FOR EvTAFv4
    addpath(genpath('/bluejay4/lucas/Dropbox/SCIENCE/code/MATLAB/chronux_2_12'));
    
elseif version ==0;
    
    rmpath(genpath('/bluejay4/lucas/Dropbox/SCIENCE/code/MATLAB/chronux_2_12'));
    
end

