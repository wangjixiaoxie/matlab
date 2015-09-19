function [ff_mat]=ffscreen(note, Fwin_low, Fwin_high,start_time,sample_size);

% function [ff_mat]=ffscreen(note, Fwin_low, Fwin_high,start_time,sample_size);
% when running melcorr_ff_pinterp or related functions, i want to screen
% the data to make sure the ff calculations match what your eye tells you.
% this should be run after pick_notes; therefore, there should be a
% batchfile of filt files and spect files. 

% messily take code from spect_sheet, obsscreen and melcorr_ff_pinterp
 
% arguments similar to those in melcorr_ff_pinterp (in ms); only use ms from start
% of syllable (not percent of syllable)


ff_mat=[];

%convert to seconds
sample_size=sample_size/1000;
start_time=start_time/1000;

mark='';   

%for bandpass
F_low=500;
F_high=8000;

% for resampling
d_Fs=1000;

%window length for smoothing prior to segmentation (in ms)
sm_win=2;

%set up card_sheet window
h_main=figure;
 

%get name of metafile containing file names

meta_fid = -1;
metafile = 0;

   disp('select batchfile');
   [metafile, pathname]=uigetfile('*','select batchfile')
   meta_fid=fopen([pathname, metafile]);
   if meta_fid == -1 | metafile == 0
      disp('cannot open file' )
      disp (metafile)
   end
Ten_notes=[];
for i=1:10

  %get soundfile name
   soundfile = fscanf(meta_fid,'%s',1);
   if isempty(soundfile)
      disp('End of soundfiles')
      break
   end

  data=[];  

    %if song exists, get it,and load song
   if exist(soundfile) 
      disp('loading files...');
      %filetype=eval(['channel',num2str(i)]);
      [datain, Fs] = soundin('',soundfile,'filt');
      datain=makerow(datain);
      %[data] = [data;datain];
            
   else
      disp(['file ', soundfile,' doesn''t exist on this path?']);
      
   end

   filt_file=[soundfile];
   disp('loading filtered song...');
     [filtsong, Fs] = read_filt(filt_file);
     default_Fs = Fs;   %if Fs was read, it is the new default value (takes precedence over notefile if discrepant
   Ten_notes=[Ten_notes filtsong];
   ff_mat=Ten_notes;
   
end

