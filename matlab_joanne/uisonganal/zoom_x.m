function [left, right, button, flag] = zoom_x

%zoom in to x axes of current axes
%select left and right boundaries with left and right mouse buttons
%hit any key when done (q to quit without updating)
%returns flag of 0 if there has been no change, else 1
%returns the key that was used on exit in "button" 
%returns new left and right axis values


%get new left and right values
 [left, right, button, flag] = get_xrange;

%if there are changes, redisplay
 if flag == 1;
   xlim=[left, right];
   set(gca,'xlim', xlim)
 end
    
