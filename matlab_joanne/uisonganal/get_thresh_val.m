function   [threshold, change_flag]=get_thresh_val

%interactive program for selecting threshold for note segmentation

%global variables
global smooth Fs threshold d_song time_d_song
global win_percent

%other variables
num_bins = 200;    %number of bins for histogram display
amp_max = .01;  %how much of hist amp range to use in calculation (smaller = faster)
n_points = 32000;  %maximum number of points to use for histogram calculation
num_sdev = 5; %number of standard deviation lines to plot above mean
change_flag=0;   %indicates whether threshold has been updated
base = [];   %vector of user selected baseline values
h_amp_base = [];  %handle vector for baseline markers
cur_thresh = threshold;  % current threshold
h_amp_thresh = [];  %handle for display of current threshold in amp window
h_hist_thresh = [];  %handle for display of current threshold in hist window
h_amp_sdev = [];  %handle for plotting mean and sdev lines in amp window
h_hist_sdev = [];  %handle for plotting mean and sdev lines in hist window
mean_base = [];
sdev_base = [];

%display amplitude info
   h_thresh=figure;
   h_thresh_amp=subplot(2, 1, 2);
   %get global values to store with axes
   g_xmin=time_d_song(1);
   g_xmax=time_d_song(length(time_d_song));
   g_ymin= 0;
   g_ymax=max(d_song)*1.2;
   
   %store the global values of data limits with axes
   set(h_thresh_amp,'userdata',[g_xmin g_xmax g_ymin g_ymax]);
   %set the current display to show zoom to .5*amp_max % of data
   set(h_thresh_amp,'xlim',[g_xmin g_xmax]);
   set(h_thresh_amp,'ylim',[g_ymin g_ymax*amp_max/2]);
   disp_amp;

%draw input threshold
   h_amp_thresh=line([g_xmin, g_xmax], [threshold, threshold], 'color', [1,0,0]);

%display histogram data
   h_thresh_hist = subplot(2, 1, 1);
   %caclulate histogram over range of interest
   cut_smooth=smooth(find(smooth < max(smooth)*amp_max));
   %get about n_points from cut_smooth approximately evenly distributed
   if 2*n_points > length(cut_smooth)  %only shorten for long songs
      i_step = round(length(cut_smooth)/n_points);
      cut_index = [1:i_step:length(cut_smooth)];
      cut_smooth = cut_smooth(cut_index);
   end
   [n, x] = hist(cut_smooth, num_bins);
   bar(x,n)
   %get and store current and absolute axis limits for hist plot
   hlim = get(h_thresh_hist,'ylim');
   alim = get(h_thresh_hist,'xlim');
   g_amin=g_xmin;
   g_amax=g_xmax;
   g_hmin=hlim(1);
   g_hmax=hlim(2);
   set(h_thresh_hist,'userdata',[g_amin g_amax g_hmin g_hmax])
   %display current threshold
   h_hist_thresh=line([threshold threshold],[g_hmin g_hmax], 'color', [1,0,0]);

%allow interactive setting in either window
action='null';
subplot(h_thresh_amp)

while ~strcmp(char(action),'q') & ~isempty(action) & (action ~= 13)

  disp('options are')

  disp(' b = select baseline interval')
  disp(' c = clear baseline intervals')
  disp(' n = change number of histogram bins')
  disp(' s = set threshold as number of sd above mean')
  disp(' zx (ux)= zoom in (out) on x axis')
  disp(' zy (uy) = zoom in (out) on y axis')
  disp(' zh (uh) = zoom in (out) on histogram frequency axis') 
  disp(' button 1 to select threshold with mouse')
  disp( ' return to exit with current threshold (red)')
  disp( ' q = quit with no changes')
  disp(' ') 


  [x ,y, action] = ginput(1);

  if strcmp(char(action), 'b')
     subplot(h_thresh_amp)
     title('Select left and right boundaries of baseline interval/return')
     [left, right, button, flag] = get_xrange;
     if flag ~= 0;    %if change   
        %convert to index from ms
        left_index = ceil(left*Fs/1000);
        right_index= floor(right*Fs/1000);          
        base_in = smooth(left_index:right_index)';           
        base = [base, base_in];
        mean_base = mean(base);
        sdev_base = std(base);
        %add baseline markers at left and right with handles
        h_amp_base(length(h_amp_base)+1)=line([left left], [0 mean_base], 'color', [0, 1, 0]);
        h_amp_base(length(h_amp_base)+1)=line([right right], [0 mean_base], 'color', [0, 1, 0]);
        h_amp_base(length(h_amp_base)+1)=line([left right], [mean_base mean_base], 'color', [0, 1, 0]);   
        h_amp_base(length(h_amp_base)+1)=line([left right], [0 0], 'color', [0, 1, 0],'linewidth',3);   
        %plot line at some number of standard deviations above mean in amp window
        if ~isempty(h_amp_sdev)   %first remove old sdev lines
           delete(h_amp_sdev)        
           h_amp_sdev = []; %reset handle value
        end
        %plot mean and standard deviations
        h_amp_sdev(1) = line([g_xmin g_xmax], [mean_base mean_base], 'color', [0, 1, 0], 'linestyle', '--');
        for i = 1:num_sdev %number of standard deviations to plot
           h_amp_sdev(length(h_amp_sdev)+1)=line([g_xmin g_xmax], [(mean_base+i*sdev_base) (mean_base+i*sdev_base)], 'color', [0, 1, 0],'linestyle','--');
        end
        %plot same in hist window
        subplot(h_thresh_hist)
        if ~isempty(h_hist_sdev)   %first remove old sdev lines
           delete(h_hist_sdev)        
           h_hist_sdev = []; %reset handle value
        end
        %plot mean and standard deviations
        h_hist_sdev(1) = line([mean_base mean_base], [g_hmin g_hmax], 'color', [0, 1, 0]);
        for i = 1:num_sdev %nuber of standard deviations to plot
           h_hist_sdev(length(h_hist_sdev)+1)=line([(mean_base+i*sdev_base) (mean_base+i*sdev_base)],...
                                                    [g_hmin g_hmax], 'color', [0, 1, 0],'linestyle','--');
        end
     end
     subplot(h_thresh_amp)
     title('')
   
  elseif strcmp(char(action), 'c')   %reset baseline values, clear standard deviation lines
     subplot(h_thresh_amp);
     base=[];
     if ~isempty(h_amp_sdev)   %first remove old sdev lines from amp window
        delete(h_amp_sdev)
        h_amp_sdev = [];
     end   
     if ~isempty(h_amp_base)
        delete(h_amp_base)
        h_amp_base = [];         
     end
     if ~isempty(h_hist_sdev)  % remove old sdev lines from hist window
        delete(h_hist_sdev)
        h_hist_sdev = [];
     end   
  
  %for now, all control is in amp window 
  %elseif action == 1      %set current threshold based on location of click in histogram window
    % subplot(h_thresh_hist);
    % xlim=get(h_thresh_hist,'xlim');
    % ylim=get(h_thresh_hist,'ylim');    
     %if y < ylim(1)
    % elseif x < xlim(1)          %move left by win_percent
     %   move_left(win_percent)
     %elseif x > xlim(2)      %move right win_percent
     %   move_right(win_percent)
     %elseif y > ylim(1)      %clicked on amplitude window
      %  cur_thresh = max(0, x);   %set current threshold
        %plot line in histogram window
      %  if h_hist_thresh ~= [];   %delete old line
       %    delete(h_hist_thresh);
       % end
        %plot new line for current threshold
       % h_hist_thresh = line([cur_thresh cur_thresh],[g_hmin, g_hmax],'color', [1,0,0])  
        %plot line on amplitude axis
       % subplot(h_thresh_amp);
      %  if h_amp_thresh ~= [];   %delete old line
      %     delete(h_amp_thresh);
      %  end
        %plot new line for current threshold
      %  h_amp_thresh = line([g_xmin, g_xmax],[cur_thresh cur_thresh],'color', [1,0,0]);
    % end
    % subplot(h_thresh_amp);

  elseif ~isempty(action) & action == 1      %set current threshold based on location of click in amplitude window
     subplot(h_thresh_amp);
     xlim=get(h_thresh_amp,'xlim');
     ylim=get(h_thresh_amp,'ylim');    
     if y > ylim(2)              %don't do anything if above window
     elseif x < xlim(1)          %move left by win_percent
        move_left(win_percent)
     elseif x > xlim(2)      %move right win_percent
        move_right(win_percent)
     elseif y < ylim(2)      %clicked on amplitude window
        cur_thresh = max(0, y);   %set current threshold
        if ~isempty(h_amp_thresh)   %delete old line
           delete(h_amp_thresh);
        end
        %plot new line for current threshold
        h_amp_thresh = line([g_xmin, g_xmax],[cur_thresh cur_thresh],'color', [1,0,0])
        %plot line in histogram window
        subplot(h_thresh_hist);
        if ~isempty(h_hist_thresh)   %delete old line
           delete(h_hist_thresh);
        end
        %plot new line for current threshold
        h_hist_thresh = line([cur_thresh cur_thresh],[g_hmin, g_hmax],'color', [1,0,0])  
     end
     subplot(h_thresh_amp);
     
     disp(['current threshold is ', num2str(cur_thresh)])
     if ~isempty(mean_base) & ~isempty(sdev_base)
       disp([num2str((cur_thresh-mean_base)/sdev_base),' sd above mean'])
     end
     
  elseif strcmp (char(action),'z')     
      %get input from keyboard or mouse
      [t,y,key]=ginput(1);        
      if strcmp(char(key), 'x')  %zoom on x axis of amplitude plot 
         %display x zoom message and zoom on amplitude x axis
          subplot(h_thresh_amp);
          title('X-Zoom mode,select with mouse/ret, q to cancel','color',[0, 0, 1])
          [xmin, xmax, button, flag] = zoom_x;
          title('')
      elseif strcmp(char(key),'y') %zoom on y axis of amplitude plot
         %display y zoom message and zoom on amplitude yaxis
          subplot(h_thresh_amp);
          title('Y-Zoom mode,select with mouse/ret, q to cancel','color',[0, 0, 1])
          [ymin, ymax, button, flag] = zoom_y;
          title('')
          %also zoom amplitude portion of hist plot
          subplot(h_thresh_hist);
          set(h_thresh_hist,'xlim',[ymin*2 ymax*2]);                    
      elseif strcmp(char(key), 'h')
         %display x zoom message and zoom on amplitude x axis
          subplot(h_thresh_hist);
          title('H-Zoom mode,select with mouse/ret, q to cancel','color',[0, 0, 1])
          [xmin, xmax, button, flag] = zoom_x;
          title('')
     % elseif strcmp(char(key),'a') 
         %display y zoom message and zoom on amplitude yaxis
         % subplot(h_thresh_hist)
         % title('A-Zoom mode,select with mouse/ret, q to cancel','color',[0, 0, 1])
         % [ymin, ymax, button, flag] = zoom_y;
         % title('')        
      end    
         
  elseif strcmp (char(action),'u')
    %get input from keyboard or mouse
    [t,y,key]=ginput(1);         
    if strcmp(char(key), 'x')             
       set(h_thresh_amp,'xlim',[g_xmin g_xmax])
     elseif strcmp(char(key),'y')  | strcmp(char(key),'a') 
       set(h_thresh_amp,'ylim',[g_ymin g_ymax])
       set(h_thresh_hist,'xlim',[g_amin g_amax])
     elseif strcmp(char(key),'h') 
       set(h_thresh_hist,'ylim',[g_hmin g_hmax])                
     elseif strcmp(char(key),'a') 
       
    end
  
  end

end



%return changed values only if threshold was altered
if (isempty(action) | (action == 13)) & (cur_thresh ~= threshold)     
    threshold = round(cur_thresh*1000)/1000;  %rounding added so that files can be opened on mac
    change_flag = 1;
else 
    change_flag = 0;
end   

%delete window, redisplay after exiting
delete(h_thresh)




