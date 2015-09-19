function uipsdanal

%function to analyze psd of selecterd song segment
%

global h_main h_main_amp
global smooth Fs threshold min_int min_dur sm_win
global onsets offsets labels
global filtsong soundfile
global percent_from_start ms_from_start psd_label sampled_dur


%size of song segment ot get (centered on left line=yellow when 's' key is hit insted of return

sample_size=16;

%get x_range for analysis
 %display spect selection message on amplitude
 subplot(h_main_amp);
 title('spectral analysis, select L/R boundaries with mouse/ret, q to cancel','color',[0, 0, 1])
 [left,right,button,flag, percent_from_start, ms_from_start, psd_label]=psd_get_xrange;
 title('');

%if there are changes, proceed
if flag == 1;
  
  %get relevant piece of song

  %case 1 's' pushed -> sample_size msec centered on right cursor
  if strcmp(char(button), 's')    
    right=left+.5*sample_size;
    left=left-.5*sample_size;
    disp(['sampled ',num2str(sample_size),' msec centered on left cursor']);
  end
    
  %case 2: otherwise take what is between cursors
  sampled_dur=right-left;  
 
  %get indexes
  left_indx=ceil(left*Fs/1000);
  right_indx=floor(right*Fs/1000);
  right_indx=max(left_indx+1,right_indx);
 
    
  song_seg=filtsong(left_indx:right_indx);
  

  %now run psdanal
%  [best_freq, best_aratio,best_pratio, modulation, ampspect, f_centers]=psdanal(song_seg,Fs);
  [best_freq, best_aratio,best_pratio, ampspect, f_centers]=psdanal(song_seg,Fs);


  %within psdanal or upon return need to reset active window to uisonganal window.
  
end
 



