f=filesep;
if strcmp('\',f);
elseif strcmp('/',f);
   disp('Platform = Unix');
   path(path,'/home/twarren/matlab/');
   path(path,'/home/twarren/matlab/evren');
   path(path,'/home/twarren/matlab/evsonganaly');
   path(path,'/home/twarren/matlab/mfiles');
   path(path,'/home/twarren/matlab/evren/uievtafsim');
   
   %path('/home/evren/evsonganaly/',path);
   %path('/home/evren/matlab/uievtafsim/',path);
   %path('/home/evren/matlab/uievautolabel/',path);
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
   %path('/home/evren/matlab',path);
end
clear f;
