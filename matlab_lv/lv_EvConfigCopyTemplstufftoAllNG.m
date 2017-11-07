function lv_EvConfigCopyTemplstufftoAllNG(EvConfigFile,copyfrom,copyto)

%  copies CntRng (i.e. threshold, min max birdtaf etc and counter logic
%from first notegroup to other notegroups
% template and column names should already be set by 
%AddTemplatesToEvConfig(EvConfigFile,varargin

[ND,OP]=ReadEvTAFv4ConfigFile(EvConfigFile,0);

nngs = size(ND,2);

log = ND(copyfrom).CntLog;
cntrng = ND(copyfrom).CntRng;

for i = copyto
   
  ND(i).CntLog = log;
  ND(i).CntRng = cntrng;
  
end

%lena
%modify name
% EvConfigFile = [EvConfigFile '_mod'];

dotdot = strfind(EvConfigFile, '.ev');
config_fname = [EvConfigFile(1:dotdot-1) '_mod.evconfig2'];

disp(['Rewriting config file named : ',config_fname]);
WriteEvTAFv4ConfigFile(config_fname,ND,OP);
return;