function load_n_plot_spec_batch(filt)

% function load_n_plot_spec_batch(filt)
% 
% This function uses a batch file to tell it the names of files to
% load.  These files can be in any directory.  It then plots all
% the spectrograms in their own figure.
%
% Dependancies: plotspec filecb.
%



if nargin==0
 [bfile, bpath]=uigetfile('*','Select batch file of song names');
else
  batch_file_filt = filt;
  [bfile, bpath]=uigetfile(batch_file_filt,'Select batch file of songnames');
end

% Get initial directory.  
initdir = pwd;

% Get the results directory
btitle = sprintf('Dir to load songs from batch file %s',bfile);
song_dir = filecb('browsedir',initdir,btitle);



bfid = fopen([bpath,bfile]);
ifile = 1;
disp('What is type of sound file? [w]')
disp(' b = binary from mac')
disp(' w = wavefile (i.e. cbs)')
disp(' d = dcpfile')
disp(' f = foosong/gogo (sets Fs = 32000)')
disp(' o = observer file (last/song channel)')
disp(' o1r = observer file (second to last)')

filetype = 'null';
while strcmp(filetype,'null')
  temp=input(' ','s');
  if strncmpi(temp,'b',1);
    filetype = 'b';
  elseif strncmpi(temp,'w',1) | isempty(temp)
    filetype = 'w';
  elseif strncmpi(temp,'d',1)
    filetype = 'd';  
  elseif strncmpi(temp,'f',1)
    filetype = 'f';  
  elseif strcmpi(temp,'o')
    filetype = 'obs0r';
  elseif strcmpi(temp,'o1r')
    filetype = 'obs1r';
  else
    disp('Unacceptable! Pick again')
    filetype = 'null';
  end
end

filetype



while 1
    % To ** written by BD Wright.
    fname = fgetl(bfid);
    % Check for whitespace here!
    spaceflag = 0;
    if isspace(fname)
      spaceflag = 1;
    end
    if (~ischar(fname)) | isempty(fname) | spaceflag
      disp('End of batch file reached.');
      break
    end
    % **

    % get .not.mat filename 
    song_file = fullfile(song_dir,fname)
    
    % Load this file (if it exists) and extract it's contents...
    if exist(song_file,'file')
      
     filetype
    [song,fs]=soundin(song_dir,fname,filetype);
    
     fignum = extractfignum(song_file)
     eval(['figure(',fignum,')']);
     plotspec(song,fs,fname)
     
   end
 end