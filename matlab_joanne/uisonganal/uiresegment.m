function uiresegment

%allows resegmentation of song
%should check data status here

global h_main_amp
global smooth Fs threshold min_int min_dur sm_win
global onsets offsets labels 
global filtsong

%get new values for segmentation
[threshold, min_int, min_dur, sm_win, change_flag]=get_segvals( Fs, threshold, min_int, min_dur, sm_win)

%if changed, resegment
if nnz(change_flag) > 0    %if something is changed           
  if change_flag(4) == 1  %if sm_win is changed, recalculate smooth 
    %calculate square of signal: proportional to power
    squared_song = filtsong.^2;
    %smooth the rectified song
    len=round(Fs*sm_win/1000);                      
    h=ones(1,len)/len;
    smooth=conv(h, squared_song);
    offset=round((length(smooth)-length(filtsong))/2); %get rid of convolution induced offset
    smooth=smooth(offset:length(filtsong)+offset);
  end
  %resegment
  [onsets, offsets] = segment(smooth, Fs, min_int, min_dur, threshold);               
  %check labels & erase if there has been a change in note number
  if length(labels) ~= length(onsets)
    labels = num2str(zeros(size(onsets)))';   
  end         
end
%redisplay amp data
subplot(h_main_amp);
disp_song;
        

