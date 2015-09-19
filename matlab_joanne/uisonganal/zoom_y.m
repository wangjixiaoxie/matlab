function [bottom, top, button, flag] = zoom_y

%zoom in to y axes of current axes
%select top and bottom boundaries with left and right mouse buttons
%hit any key when done (q to quit without updating)
%returns flag of 0 if there has been no change, else 1
%returns the key that was used on exit in "button" 
%returns new top and bottom axis values


%get new top and bottom values
 [bottom, top, button, flag] = get_yrange;

%if there are changes, redisplay
 if flag == 1;
   ylim=[bottom, top];
   set(gca,'ylim', ylim)
 end
 

            
   
