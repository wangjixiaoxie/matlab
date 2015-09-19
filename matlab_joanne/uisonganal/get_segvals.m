function [threshold, min_int, min_dur, sm_win, change_flag] = get_segvals(Fs, threshold, min_int, min_dur, sm_win)

%queries for changes to values used in segmentation of songfile

%flags to keep track of changes
thresh_i=1;
int_i=2;
dur_i=3;
win_i=4;
change_flag=[0 0 0 0];
t_min_dur=min_dur;
t_min_int=min_int;
t_thresh=threshold;
t_sm_win=sm_win;
  

while 1

  disp('Modify which parameters?')
  disp(' t = get new threshold')
  disp([' d = min_dur [' num2str(t_min_dur) ']'])
  disp([' i = min_int [' num2str(t_min_int) ']'])
  disp([' s = sm_win [' num2str(t_sm_win) ']'])
  disp(' r = recalculate')
  disp(' q = quit (no changes)')
  
  [x,y,action] = ginput(1);
     
  if strcmp(char(action),'d')
    %change min_dur
    disp(['enter min_dur (in ms): [' num2str(t_min_dur) '] '])
    [x,y,in]=ginput;
    if (isempty(in))   %no number entered 
      disp('no change to min_dur')
    else            %convert to number
      in=char(in);
      in=eval(in);
      if (in == t_min_dur)  %check whether number has been changed
         disp('no change to min_dur')
      else                  %if changed, assign value and set flag
         t_min_dur = in;
         change_flag(dur_i)=1;
         disp(['min_dur = ' num2str(t_min_dur)])
      end               
    end   
    

  elseif strcmp(char(action),'i')
    %change int_dur
    disp(['enter int_dur (in ms): [' num2str(t_min_int) '] '])
    [x,y,in]=ginput;
    if (isempty(in))   %no number entered 
      disp('no change to min_int')
    else            %convert to number
      in=char(in);
      in=eval(in);
      if (in == t_min_int)  %check whether number has been changed
         disp('no change to min_int')
      else 
         t_min_int = in;
         change_flag(int_i)=1;
         disp(['int_dur = ' num2str(t_min_int)])
      end 
    end   
    
  elseif strcmp (char(action),'s')
    %change smoothing window
    disp(['enter sm_win (in ms): [' num2str(t_sm_win) '] '])
    [x,y,in]=ginput;
    if (isempty(in))   %no number entered 
      disp('no change to sm_win')
    else            %convert to number
      in=char(in);
      in=eval(in); 
      if (in == sm_win)  %check whether number has been changed
         disp('no change to sm_win')
      else 
        t_sm_win = in;
        change_flag(win_i) = 1;
        disp(['sm_win = ' num2str(t_sm_win)])
      end
    end   
    
  elseif strcmp(char(action),'t')
    %get new value of threshold from interactive program
    [threshold, flag] = get_thresh_val;
    if flag == 1
       change_flag(thresh_i)=1;
       return   
    end
  elseif strcmp(char(action), 'q')
    change_flag = [0 0 0 0];
    return    
  
  elseif strcmp(char(action), 'r')
    min_dur=t_min_dur;
    min_int=t_min_int;
    threshold=t_thresh;
    sm_win=t_sm_win;
    return      
  
  else
    disp('Invalid option')
  
  end

end
