function   make_current(h_fig, status)

%makes the specified figure into the current figure
% and updates data variables depending on status setting
% in all cases, resets data_status of all open figures appropriately
%
% for status = 1, 
% -> brings figure to front
% if figure is not already current
% -> saves notedata from old figure
% -> resets global variables:
%    h_main, h_main_amp, h_main_spect
% -> reads in notefile for current figure which resets global variables:
%    onsets, offsets, labels, Fs, threshold,  min_int, min_dur, sm_win
% -> retrieves d_song and time_d_song from amplitude plot
%
% for status = 2, does above plus 
% -> retrieves filtsong and smooth
%
% following option not yet implemented
% for status = 3, does above plus
% -> retrieves spectrogram data:
%   idx_spect etc 

%global variables
global h_main h_main_amp h_main_spect
global onsets offsets labels Fs threshold  min_int min_dur sm_win
global main_win_code
global d_song d_time_song
global soundfile path_notefile
global filtsong smooth path_filtfile path_songfile path_spectfile
global do_zc d_zc
global zc zcwin

%indexes into main figure userdata
%these just facilitate mneumonics of storage and retrieval of values in userdata
%values are set in songanal
global imain_win_code idata_status ih_main_amp ih_main_spect ih_fname 
global ih_path_songfile ih_path_notefile ih_path_filtfile ih_path_spectfile ih_amp_plot ih_zc_plot

watchon;
%---------------------------------------------------
%  Ensure that old data is saved and note data is current
%---------------------------------------------------
if status == 1  &  h_fig == h_main
  figure(h_main);
  watchoff;
  return
  %if fig is not current or need more than notes, first update notes
  % save old figure's data, get new figure's data
else
  %save notedata from old figure
  if (ischar(soundfile) & ~isempty(soundfile) & ~isspace(soundfile) )
    note_file=[soundfile,'.not.mat'];
    save_data(note_file, path_notefile, Fs, onsets, offsets, labels, threshold, min_int, min_dur, sm_win);
  end
  %get userdata from h_fig (contains figure id and handles)
  userdata = get(h_fig,'userdata');
  %get current data status
  current_status = userdata(idata_status);
  %set h_main to current figure and bring to front
  h_main = h_fig;
  figure(h_main); 
  %set other handles
  h_main_amp = userdata(ih_main_amp);
  h_main_spect = userdata(ih_main_spect);
  %retrieve soundfile
  h_fname = userdata(ih_fname);
  soundfile = get(h_fname,'string');
  %get path to notefile  
  h_path_notefile = userdata(ih_path_notefile);
  path_notefile = get(h_path_notefile,'string');
  %load notefile
  note_file=[soundfile,'.not.mat'];
  if exist(fullfile(path_notefile, note_file),'file')
    disp('loading notefile...');     
    load(fullfile(path_notefile, note_file));
  else 
    %this should only happen if files are deleted/moved while running
    disp('make_current: Warning! unable to load notefile')
  end     
  %retrieve d_song and time_d_song from amplitude axis
  %     h_amp_plot = userdata(ih_amp_plot);
  h_amp_plot = findobj(h_main_amp,'Tag','AmpPlot');
  d_song = get(h_amp_plot,'ydata');
  time_d_song = get(h_amp_plot,'xdata');
  if do_zc
    h_zc_plot = findobj(h_main_amp,'Tag','ZCPlot');
    d_zc = get(h_zc_plot,'ydata');
  end
  %update current status
  current_status = max(1, current_status);   
end

%-------------------------------------------------------
%if necessary, get filtsong and smooth
%-------------------------------------------------------
if current_status < status & status >= 2     
  %retrieve path to filtfile
  h_path_filtfile = userdata(ih_path_filtfile);
  path_filtfile = get(h_path_filtfile,'string');
  %if filtfile exists, read it and smooth it
  filt_file=fullfile(path_filtfile, [soundfile, '.filt']);
  if exist(filt_file,'file')     
    %read filtsong
    disp('loading filtered song...');
    [filtsong, Fs] = read_filt(filt_file);
    %calculate square of signal: proportional to power
    squared_song = filtsong.^2;
    %smooth the rectified song
    len=round(Fs*sm_win/1000);                      
    h=ones(1,len)/len;
    smooth=conv(h, squared_song);
    offset=round((length(smooth)-length(filtsong))/2); %get rid of convolution induced offset
    smooth=smooth(offset:length(filtsong)+offset);

    % Zero crossings
    if do_zc
      zcwinlen=round(Fs*zcwin/1000);                      
      % This gets average # of zero crossings per sample
      zcwindow=ones(1,zcwinlen)/zcwinlen;
      zc = zero_crossings(filtsong,zcwindow);
      offsetzc=round((length(zc)-length(filtsong))/2); %get rid of convolution induced offset
      zc = (Fs/2)*zc(1+offsetzc:length(filtsong)+offsetzc); 
    end

  
  else
    disp('make_current: Warning! filtsong does not exist')
    disp('   filtsong and smooth have not been updated')
    disp('    need to read in raw song and filter in make_current')
  end            
  %update current status
  current_status = max(2, current_status);
end

%------------------------------------------------------
%update spect data if neccessary
%------------------------------------------------------
if current_status < status & status == 3
  disp('make_current: Warning! spectrogram updating not curreently implemented')
  disp('   spectrogram has not been updated')
  disp('    need to add code to make_current')
  disp('    have a nice day!')
end  

%------------------------------------------------------
% update data status of all open windows
%------------------------------------------------------
windows = get(0,'children');
for i = 1:length(windows)
  userdata = get(windows(i),'userdata');
  if length(userdata > 0)
    if userdata(imain_win_code) == main_win_code;
      userdata(idata_status) = 0;
      set(windows(i),'userdata',userdata);
    end
  end
end

%set data status correctly for current window
userdata = get(h_main,'userdata');
userdata(idata_status) = current_status;
set(h_main,'userdata', userdata); 

watchoff;





