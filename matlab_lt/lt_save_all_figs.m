function lt_save_all_figs
%% LT 4/7/15 - saves all open figures as figure1, figure2,...
% note, order might not necessarily be correct, and will start from 1, ...


h = get(0,'children');
for i=1:length(h)
  saveas(h(i), ['figure' num2str(length(h)+1-i)], 'fig');
end


