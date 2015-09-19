function [left,right,button,flag, percent_from_start, ms_from_start, curr_label]=psd_get_xrange


%modified from get_xrange to return % of the way through note when
% used while uisonganal is running (to allow consistent picking of note parts
% for spectral analysis

global onsets offsets labels

%allows user to select range of x axis values with mouse from current axis
%select left value with left button
%select right value with middle or right button
%hit any key (except q) when done
%q to quit without changing values
%the key value is returned as button
%flag is set to 1 if legal changes were made

%initial values from current axes
xlim=get(gca,'xlim');
ylim=get(gca,'ylim');
xmin=xlim(1);
xmax=xlim(2);
ymin=ylim(1);
ymax=ylim(2);

%default values
left=xmin;
right=xmax;
button=1;
flag=0;  % no change state

%plot default lines

h_left=line([left left],[ymin ymax],'color',[1,1,0],'LineStyle','--');
h_right=line([right right],[ymin ymax],'color',[1,0,0],'LineStyle','--');

while ~isempty(button) & ((button == 1) | (button == 2) | (button == 3))

   [xval,y,button]=ginput(1)
   
   if ~isempty(button)
     if (button == 1)
       left=max(xval,xmin);
       left=min(left,xmax);
       set(h_left,'xdata',[left left])
       %what percentage through the note are we?
       ons_before_ind=find(onsets<=left);
       curr_note_ind=ons_before_ind(length(ons_before_ind));
       curr_on=onsets(curr_note_ind);
       curr_off=offsets(curr_note_ind);
       curr_label=labels(curr_note_ind);
       ms_from_start=left-curr_on;
       percent_from_start=100*ms_from_start/(curr_off-curr_on);
       disp(['percent from start = ',num2str(percent_from_start),' (',num2str(ms_from_start),' msec)']);
     elseif (button == 2) | (button == 3)
       right=min(xval,xmax);
       right=max(right,xmin);
       set(h_right,'xdata',[right right])     %move line
     end
   end
end


%if quiting, or illegal values, or no change leave values unchanged
%and set flag to 0

if (left == xmin) & (right == xmax)
   flag  = 0;
elseif (left >= right) | strcmp(char(button),'q')
   left=xmin
   right=xmax
   flag = 0;
else
   flag = 1
end

%delete lines

delete(h_left)
delete(h_right)




     
