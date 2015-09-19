function []=quickscreen(markfile, filetype)
%quick and dirty look at observer data:
%  open %  read each file,and display


if nargin==1
  filetype='obs0r';   
end


savefile=[markfile,'.short'];
delfile=[markfile,'.del'];
callfile=[markfile,'.long'];


%open files for saving
if exist(savefile)
  disp([savefile,' already exists']);
  disp(['will append new marks']);
end
save_fid=fopen(savefile,'a');

if exist(delfile)
  disp([delfile,' already exists']);
  disp(['will append new marks']);
end
del_fid=fopen(delfile,'a');

if exist(callfile)
  disp([callfile,' already exists']);
  disp(['will append new marks']);
end
call_fid=fopen(callfile,'a');



mark='';

%smooth_data?  1=yes anything else = no
sm_flag=1;
%sm_win specifies amount of smoothing in ms; square shaped filter
sm_win=25;
%plot stuff? 1=yes
plot_flag=0;

%resample data (speeds things up, reduces memory problems yes=1
resample_flag=1;
d_Fs=1000;  %samples per second after resampling  




%matrix for output: each row = one rectified waveform
data=[];  


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


while ~strcmp(mark,'q')

   if exist('h_main')
      delete(h_main);
   end
    %get soundfile name
   soundfile = fscanf(meta_fid,'%s',1);
   if isempty(soundfile)
      disp('End of soundfiles')
      break
   end


  

  %if song exists, get it,
   if exist(soundfile) 
      disp('loading file...');
      [sounddata, Fs] = soundin('',soundfile,filetype);
   else
      disp(['file ', soundfile,' doesn''t exist on this path?']);
      
   end
  
 
 
   
  % rectify the filtfile
  %filtdata=abs(sounddata);


  %smooth the filtfile?
  %if sm_flag==1 | resample_flag==1;
   %len=round(Fs*sm_win/1000);                      
   %h=ones(1,len)/len;
   %smooth=conv(h, filtdata);
   %offset=round((length(smooth)-length(filtdata))/2); %get rid of convolution induced offset
   %smooth=smooth(offset:length(filtdata)+offset);
   %filtdata=smooth;
  %end 
   
  %if resample_flag==1
   %resample the filtfile?
   %resample song for more rapid display
   %disp('resampling for display...');
   %step=round(Fs/d_Fs);
   %d_Fs=Fs/step;
   %filtdata=filtdata(1:step:length(filtdata));
   %Fs=d_Fs;
  %end

   
  
  %calculate time indeces for sample points

  time_d=[0:length(sounddata)-1]*1/Fs;  %vector for time axis

  %plot file
  h_main=figure;
  plot(time_d,sounddata);
  title(soundfile);
  input_flag=1;
  while input_flag==1
   figure(h_main);
   [x,y,mark]=ginput(1);
   mark=char(mark);
        if strcmp(mark,'g')
            pos_main=get(h_main,'position');
            %make a spectrogram in a separate window
            h_spect=figure;
            v_offset=pos_main(4)+75;
            v_size=pos_main(4);
            h_size=pos_main(3);
            h_location=pos_main(1);
            pos_spect=[h_location pos_main(2)-v_offset h_size v_size];
            set(h_spect,'position',pos_spect);
            specgram(sounddata,512,Fs,[],384);
        elseif strcmp (mark,'p')
            play(sounddata,Fs);
        elseif strcmp(mark,'s')
            fprintf(save_fid,'%s\n',soundfile);input_flag=0;        
        elseif strcmp(mark,'d')
            fprintf(del_fid,'%s\n',soundfile);input_flag=0;    
        elseif strcmp(mark,'l')
            fprintf(call_fid,'%s\n',soundfile);input_flag=0;    
        elseif strcmp(mark,'q')
            disp('quit');input_flag=0;
        end     
  end
  if exist('h_spect'); delete(h_spect); clear('h_spect'); end
end
fclose(save_fid);
fclose(del_fid);
fclose(call_fid);
