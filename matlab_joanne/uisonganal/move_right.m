function move_flag = move_right(percent_win)

%shifts the viewing frame to the right by percent_win percentage of window size
%move_flag=1 if the axis were moved, 0 if absolute axis limit was reached

%get current axis boundaries
xlim=get(gca,'xlim');
xmin=xlim(1);
xmax=xlim(2);

%get absolute axis boundaries (i.e. data limits)
g_lim=get(gca,'userdata');
g_xmin=g_lim(1);
g_xmax=g_lim(2);

% Displaying up to end of data already, don't move any more.
if xmax >= g_xmax
  move_flag = 0
  return
end

win_size=xmax-xmin;
new_xmax=min(xmax+percent_win*win_size,g_xmax);
new_xmin=new_xmax-win_size;
if new_xmax == xmax
  move_flag = 0
else
  move_flag = 1
end 
xlim=[new_xmin new_xmax];
set(gca,'xlim',xlim);
