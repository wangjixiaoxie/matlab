function [song,Fs,name] = loadraw(filt)

% function [song,Fs,name] = loadraw(filt)
%
% This function loads raw data from a file and returns the data,
% the sampling rate and the name of the file.  FILT is a filter
% string for use in uigetfile.

if isempty(filt)
 [songfile, path_songfile]=uigetfile('*','select file')
else
  [songfile, path_songfile]=uigetfile(filt,'select file')
end
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


[song,Fs]=soundin(path_songfile, songfile, filetype);

name = songfile;