function []=quickscreenjw(batchf, filetype, markfile)
% screen files:
% quickscreenjw(markfile, filetype, batchf)
% markfile is prefix for output batchfiles/imagefiles, filetype=
% 'w' for wav, 'obs0r' for observer, batchf is input batchfile
%  


fprintf('%');
if nargin==1
  filetype='obs0r';   
end


savefile=[markfile,'.save'];
delfile=[markfile,'.del'];
callfile=[markfile,'.call'];
image_array=dir([markfile,'*-image*jpg']);
if size(image_array)>0
icount=str2num(image_array(length(image_array)).name((findstr('image',image_array(length(image_array)).name)+5):(findstr('jpg',image_array(length(image_array)).name)-2)));
else
icount=0;
end
clear image_array;



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

%   disp('select batchfile');
%   [metafile, pathname]=uigetfile('*','select batchfile')
%   meta_fid=fopen([pathname, metafile]);
%   if meta_fid == -1 | metafile == 0
%      disp('cannot open file' )
%      disp (metafile)
%   end
meta_fid=fopen(batchf);

while ~strcmp(mark,'q')
   if exist('h_spect')
	   delete(h_spect);
   end
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
      disp(['loading ',soundfile]);
      [sounddata, Fs] = evsoundin('',soundfile,filetype);
      sounddata = highpass(sounddata,200,Fs);
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

  %plot file ^M
  h_main=figure;
  plot(time_d,sounddata);
  title(soundfile);
   pos_main=get(h_main,'position');
  set(h_main,'position',[pos_main(1) pos_main(2)+10 pos_main(3) pos_main(4)]);
  input_flag=1;
  while input_flag==1
    figure(h_main);
    stupid=waitforbuttonpress; %this command must run before keypresses, so that characters are not passed to command line
    mark=get(h_main,'CurrentCharacter');
    if strcmp(mark,'.')    %indicates that a spectrogram has been generated and that keypress has been done with the spectrogram as the current fig. In which case mark needs to get its value from h_spect instead of h_main.
	if exist('h_spect')
	    mark=get(h_spect,'CurrentCharacter');
    end

        end

%   if exist('h_spect')
%   if strcmp(get(h_spect,'CurrentCharacter'),'i')
%	   disp('specgram-i')
%   end
%   end

   %mark=char(mark);
        if strcmp(mark,'g')
	set(h_main,'CurrentCharacter','.')
            %make a spectrogram in a separate window
	    %set(h_main,'position',[pos_main(1) pos_main(2)+10 pos_main(3) pos_main(4)]);
	    h_spect=figure;
	    set (h_spect,'CurrentCharacter','.');
            v_offset=pos_main(4)+25;
            v_size=pos_main(4)-70;
            h_size=pos_main(3)+500;
            h_location=pos_main(1)-250;
            pos_spect=[h_location pos_main(2)-v_offset h_size v_size];
            set(h_spect,'position',pos_spect);
            specgram(sounddata,1024,Fs,[],768); ylim([0,1.6e4]);
	    %specgram(sounddata,882,Fs,882,780); ylim([0,1.6e4]);
	    caxis('auto');
	    [cax]=caxis;
	    caxis([cax(1)+(cax(2)-cax(1))*.2,cax(2)]);
	    %the below caxis works great if you are using wavread16 in evsoundin instead of wavread
	    %caxis(((2^8)-1)*[-.1, .46]);
	   %from RK's newer quickscreenr2: set(gca,'YD','n','Ylim',([0 10000]));
        elseif strcmp (mark,'i')
		set(h_main,'CurrentCharacter','.')
		if exist('h_spect')==1
		figure(h_spect)
		end
		icount=icount+1;
		imgpath=eval('pwd');
		
		if imgpath(findstr(imgpath,'/song')+6)=='/'
			loctag=imgpath(1+findstr(imgpath,'/song'):findstr(imgpath,'/song')+5);
		elseif imgpath(findstr(imgpath,'/song')+6)=='0'|imgpath(findstr(imgpath,'/song')+6)=='1'|imgpath(findstr(imgpath,'/song')+6)=='3'
			loctag=imgpath(1+findstr(imgpath,'/song'):findstr(imgpath,'/song')+6);
		else
		disp('song computer name not found in pwd path');
		loctag=['rigmachine'];
		end
		if soundfile((length(soundfile)-3):length(soundfile))=='.wav'
		image_name=[markfile,'-',soundfile((length(soundfile)-11):length(soundfile)-9),'-',soundfile((length(soundfile)-7):length(soundfile)-4),'-on',loctag,'-Image',num2str(icount)];
		elseif soundfile((length(soundfile)-3):length(soundfile))=='cbin'
		%image_name=[markfile,'-',soundfile((length(soundfile)-11):length(soundfile)-9),'-',soundfile((length(soundfile)-8):length(soundfile)-5),'-on',loctag,'-Image',num2str(icount)];
		dateind=findstr(soundfile,'_');
		dotind=findstr(soundfile,'.');
		image_name=[markfile,'-',soundfile(dateind(1)+1:dateind(2)-1),'-',soundfile(dotind(1)+1:dotind(2)-1),'-on',loctag,'-Image',num2str(icount)];
		else
		image_name=[markfile,'-',soundfile((length(soundfile)-7):length(soundfile)-5),'-',soundfile((length(soundfile)-3):length(soundfile)),'-on',loctag,'-Image',num2str(icount)];
		end
		eval(['print -djpeg ',image_name]);
        %elseif strcmp (mark,'p')
        %   disp([' ',num2str(Fs)]);
        %    play(sounddata,Fs);
        elseif strcmp(mark,'s')
            fprintf(save_fid,'%s\n',soundfile);input_flag=0;        
        elseif strcmp(mark,'d')
            fprintf(del_fid,'%s\n',soundfile);input_flag=0;    
        elseif strcmp(mark,'c')
            fprintf(call_fid,'%s\n',soundfile);input_flag=0;    
        elseif strcmp(mark,'q')
            disp('quit');input_flag=0;delete(h_main);
        end     
  end %for inputflag
  if exist('h_spect'); delete(h_spect); clear('h_spect'); end
end %for ~strcmp q
fclose(save_fid);
fclose(del_fid);
fclose(call_fid);
