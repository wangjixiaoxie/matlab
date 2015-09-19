function handseg2binary(results)

% function handseg2binary(results)
%
% This function writes binary files of zeros and one representing
% the onsets and offsets for the
% data set in the structs in RESULTS.  It writes ones between each
% onset and offset and zeros everywhere else.  
%



% Get the raw song directory
initdir = pwd;
btitle = sprintf('Which dir are the raw songs in? ');
raw_dir = filecb('browsedir',initdir,btitle);

% Get the dir where we will store the binary files
btitle = sprintf('Which dir do you want to store the binary files in? ');
bin_dir = filecb('browsedir',initdir,btitle);


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




for i = 1:length(results)

  [song,fs]=soundin(raw_dir,results{i}.filename,filetype);
  len =length(song);
  
  datfile = results{i}.filename

  onsets = results{i}.onsets'*fs/1000;
  offsets = results{i}.offsets'*fs/1000;

  boundaries = sort(round([onsets,offsets]));

  segs = bound2seg_sym(boundaries,len,1);

  clf
  plot((song-mean(song))/std(song))
  hold on
  plot(segs,'g')
  hold off
  pause;
  
  segfn = fullfile(bin_dir,[results{i}.filename,'handseg']);
  fid = fopen(segfn,'w')
  fwrite(fid,segs,'int')
  
end
