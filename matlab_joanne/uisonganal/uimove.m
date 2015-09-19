function uimove(option)

%shift figure left or right

global h_main_amp h_main_spect win_percent

ha_curr = gca;
ha = findobj(gcf,'Tag','SAAxis');
na = length(ha);

if strcmp(option,'left')
  for i=1:na
    subplot(ha(i));
    move_left(win_percent);
  end
elseif strcmp(option,'right')         
  for i=1:na
    subplot(ha(i));
    move_right(win_percent);
  end
end

set(gcf,'CurrentAxes',ha_curr);
