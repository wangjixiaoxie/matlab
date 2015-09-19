%to make uigetfile work
%setappdata(0,'UseNativeSystemDialogs',false);
f=filesep;
if strcmp('\',f);
elseif strcmp('/',f);
   disp('Platform = Unix');
   path('/home/twarren/matlab',path);
   
   %AApath('/net/swallow7/blab/matlab/msb/altfeed',path);
   %AApath('/net/swallow7/blab/matlab/msb/altfeed/bird',path);
   %AApath('/net/swallow7/blab/matlab/msb/altfeed/data',path);
   %AApath('/net/swallow7/blab/matlab/msb/altfeed/data_z05',path);
   %AApath('/net/swallow7/blab/matlab/msb/cardanal',path);
   %AApath('/net/swallow7/blab/matlab/msb/noteanal',path);
   %AApath('/net/swallow7/blab/matlab/msb/songanal/utils',path);
   %AApath('/net/swallow7/blab/matlab/msb/songanal/xcorr',path);
   %AApath('/net/swallow7/blab/matlab/msb/songanal/utils/sylanal',path);
   %AApath('/net/swallow7/blab/matlab/msb/sound_utils',path);
   %AApath('/net/swallow7/blab/matlab/msb/utils',path);
   %AApath('/net/swallow7/blab/matlab/msb',path);
   %CHANGE THIS path('/net/swallow7/matlab/clint',path); TO ->
   %AApath('/net/swallow7/blab/matlab/clint',path);
   %AApath('/net/swallow7/blab/matlab/uisonganal/spect',path);
   %AApath('/net/swallow7/blab/matlab/uisonganal',path);
   % REMOVED path('/net/swallow7/blab/matlab/user');
  
end
clear f;
