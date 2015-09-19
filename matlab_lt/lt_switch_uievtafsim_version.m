function lt_switch_uievtafsim_version(version);
%% LT 4/9/15 - run this before running uievtafsim. this adds or removes dir for evtafv4 version to paths.(0=amp, 1=v4)
% This is necessary because both versions have same name, and also same
% name for called functions (but actually diff functions).

% inputs
% version = 0; evtaf amp
% = 1; evtafv4
%% Add + all subfolders

if version==1;
    % FOR EvTAFv4
    addpath(genpath('/bluejay4/lucas/Dropbox/SCIENCE/code/MATLAB/EvTAFv4_good_path_badremoved/'));
    
elseif version ==0;
    
    rmpath(genpath('/bluejay4/lucas/Dropbox/SCIENCE/code/MATLAB/EvTAFv4_good_path_badremoved/'));
    
end

