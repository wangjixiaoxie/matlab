function move_left(percent_win)

%shifts the viewing frame to the left by percent_win percentage of window size

%get current axis boundaries
xlim=get(gca,'xlim');
xmin=xlim(1);
xmax=xlim(2);

%get absolute axis boundaries (i.e. data limits)
g_lim=get(gca,'userdata');
g_xmin=g_lim(1);
g_xmax=g_lim(2);

win_size=xmax-xmin;
xmin=max(g_xmin,xmin-percent_win*win_size);
xmax=xmin+win_size;
xlim=[xmin xmax];
set(gca,'xlim',xlim);
