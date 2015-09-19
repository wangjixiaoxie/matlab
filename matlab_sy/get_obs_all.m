function [data, time_d]=get_obs_all(filetype1);

%  7/9/01 modified to work for filt files with corresponding .not.mat files; old version for foosong now called:
% get_filts_labelled_foo
% like get_filts, except that it pulls out sound surrounding 1st match to specified note label (+/- pre/post)
%  requires corresponding .mat file
%  quick and dirty look at neural (or song) data:
%  open batch_file w/ filt_file names
%  read each file, rectify it, (smooth it? resample?) store in data
%  calculate average, std, etc
%  plot something... if plot_flag==1;
%data_sfx is file out of which data should be pulled (.1  .3, .4, etc);
% optional_batch is optional name of batch file of .filt files; if present, then user is not queried for input batch


%should we ask user for batchfile name?
%if nargin==2  % not if batchfile is specified as arg
 % menu_in=0;
 %else
 menu_in=1;
%end


%smooth_data?  1=yes anything else = no
sm_flag=1;
%smooth_win specifies amount of smoothing in ms; square shaped filter
smooth_win=5;
%plot stuff? 1=yes
plot_flag=0;

%which occurence of match to 'notes' should be pulled out?
number=1;

%SND_SFX='1';  %suffix for .not.mat files (corresponds to soundfiles saved w/ .1 ending


%rectify data?
rectify_flag=1;


%resample data (speeds things up, reduces memory problems yes=1
resample_flag=1;
d_Fs=1000;  %samples per second after resampling  




%matrix for output: each row = one rectified waveform
data=[];  
trials=[];



%get name of metafile containing filtfile names
if menu_in==1
  meta_fid = -1;
  metafile = 0;
  while meta_fid == -1 | metafile == 0 | metafile == ''
   disp('select batchfile');
   [metafile, pathname]=uigetfile('*','select batchfile')
   meta_fid=fopen([pathname, metafile]);
   if meta_fid == -1 | metafile == 0
      disp('cannot open file' )
      disp (metafile)
   end
  end
else
  meta_fid=fopen(optional_batch);
  if meta_fid == -1
      disp('cannot open file' )
      disp(optional_batch)
  end    
end

while 1

  %get soundfile name
   filename = fscanf(meta_fid,'%s',1)
   if isempty(filename)
      disp('End of filtfiles')
      break
   end

  %check if filename ends in '.filt' ; if so, strip .filt
   
   filt_idx=findstr('.filt',filename);
   if ~isempty(filt_idx)
     filename=filename(1:filt_idx(1)-1);         
   end
   
   datafile=filename;
  
  
  
   %get sound
   [filtdata, Fs]=soundin('',datafile, filetype1);
   %[filtdata2, Fs2]=soundin('',datafile, filetype2);
   
   % rectify the filtfile
   if rectify_flag==1
          filtdata=abs(filtdata);
      else
             
   end
    


   %smooth the filtfile?
   if sm_flag==1 | resample_flag==1;
    len=round(Fs*smooth_win/1000);                      
    h=ones(1,len)/len;
    smooth=conv(h, filtdata);
    offset=round((length(smooth)-length(filtdata))/2); %get rid of convolution induced offset
    smooth=smooth(offset:length(filtdata)+offset);
    filtdata=smooth;
   end 
   
   if resample_flag==1
    %resample the filtfile?
    %resample song for more rapid display
    %disp('resampling for display...');
    step=round(Fs/d_Fs);
    d_Fs=Fs/step;
    filtdata=filtdata(1:step:length(filtdata));
    Fs=d_Fs;
   end

   %add filtfile to the data structure
   %make sure we are working w/ a row vector:
   filtdata=makerow(filtdata);
   
  
  
   data=[data;filtdata];
   
  end

end

%calculate time indeces for sample points

if resample_flag==1
    Fs=d_Fs;
end


time_d=[0:length(filtdata)-1]*1000/Fs;  %vector for time axis


%calculate mean, std, stderr, cov, etc...

%[nrows,ncols]=size(data);

%mean_data=mean(data);
%std_data=std(data);
%se_data=std_data/sqrt(nrows);
%cov_data=std_data./mean(data);


%plot something useful...
if plot_flag==1;
	h_raw=figure;
	 hold on;
 	mm_data=mean(mean_data);
 	p_offset=5*median(abs(mean_data-mm_data));
 	for i=1:nrows
 	  plot(time_d,data(i,:)+i*p_offset); 
 	end


	h_summary=figure;
	hold on;
	plot(time_d,mean_data,'b');
	plot(time_d,mean_data+se_data,'r:');
	plot(time_d,mean_data-se_data,'r:');
	xlabel('time(ms)')
	title([metafile,'; n=',num2str(nrows)]);
end
