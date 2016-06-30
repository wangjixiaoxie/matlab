lt_figure; hold on;


%% ==== first, enter data here [MUSC/PBS, separation(2-dir) or sum(1-dir)]
X1=[253/296 281/380, 215/227, 390/401]; % same direction
X1=[256/282 281/380, 215/227, 390/401]; % same direction (using same rd23 syls as in 2 dir)

X2=[-4/190, 179/325, 201/293, 124/188, 74/172, 39/190, 226/297, 168/246]; % diff dir, including all(including pu11 2, which is actually not same type
X2=[-4/190, 201/293, 124/188, 74/172, 39/190, 226/297, 168/246]; % diff dir, excluding pu11 2, which is actually not same type


%% test significance for diff in bin 3, not paired


[~, p]=ttest2(X1, X2);

% ================== all, including pu11-2 which is diff type
lt_subplot(4,2,1); hold on
ylabel('consolidation');
title('unpaired, including pi11-2, diff type');
X1=[253/296 281/380, 215/227, 390/401]; % same direction
X2=[-4/190, 179/325, 201/293, 124/188, 74/172, 39/190, 226/297, 168/246]; % including all(including pu11 2, which is actually not same type
[~, p]=ttest2(X1, X2);

plot(1.1, X1, 'ok');

plot(2.1, X2, 'ok');

lt_plot_bar([1 2], [mean(X1) mean(X2)], {'Errors', [lt_sem(X1) lt_sem(X2)]});
lt_plot_text(1.5, 1.1*max([X1 X2]), num2str(p, '%3.2g'), 'r')
xlim([0 3]);
ylim([0 1.5])
line(xlim, [1 1], 'LineStyle', '--', 'Color','k')

% ======================== all [same type], (i.e excluding pu11-2 which is diff type
lt_subplot(4,2,2); hold on
ylabel('consolidation');
title('unpaired, only same type(excluding pu11-2)');
X1=[253/296 281/380, 215/227, 390/401]; % same direction
X2=[-4/190, 201/293, 124/188, 74/172, 39/190, 226/297, 168/246]; % excluding pu11 2, which is actually not same type
[~, p]=ttest2(X1, X2);

plot(1.1, X1, 'ok');
plot(2.1, X2, 'ok');

lt_plot_bar([1 2], [mean(X1) mean(X2)], {'Errors', [lt_sem(X1) lt_sem(X2)]});
lt_plot_text(1.5, 1.1*max([X1 X2]), num2str(p, '%3.2g'), 'r')
xlim([0 3]);
ylim([0 1.5])
line(xlim, [1 1], 'LineStyle', '--', 'Color','k')


%% test significance for diff in bin 3, PAIRED



% =========================== all, using exact same syls
lt_subplot(4,2,3); hold on
ylabel('consolidation');
title('paired (using exact same syls in both phases');
X1=[256/282 281/380, 215/227, 390/401]; % same direction (using same rd23 syls as in 2 dir)
X2=[124/188, 39/190, 226/297, 168/246]; % 
[~, p]=ttest(X1, X2);

plot([1.1 2.1], [X1; X2], '-ok');

lt_plot_bar([1 2], [mean(X1) mean(X2)], {'Errors', [lt_sem(X1) lt_sem(X2)]});
lt_plot_text(1.5, 1.1*max([X1 X2]), num2str(p, '%3.2g'), 'r')
xlim([0 3]);
ylim([0 1.5])
line(xlim, [1 1], 'LineStyle', '--', 'Color','k')





% ============================= all, rd23 one syl same, one syl diff, between phases
lt_subplot(4,2,4); hold on
ylabel('consolidation');
title('paired (one rd23 syl not overlapped)');
X1=[253/296 281/380, 215/227, 390/401]; % same direction
X2=[124/188, 39/190, 226/297, 168/246]; % 
[~, p]=ttest(X1, X2);

plot([1.1 2.1], [X1; X2], '-ok');

lt_plot_bar([1 2], [mean(X1) mean(X2)], {'Errors', [lt_sem(X1) lt_sem(X2)]});
lt_plot_text(1.5, 1.1*max([X1 X2]), num2str(p, '%3.2g'), 'r')
xlim([0 3]);
ylim([0 1.5])
line(xlim, [1 1], 'LineStyle', '--', 'Color','k')



% % =========================== all, using exact same syls [TEMP, EYEBALLING
% % DATA FOR LAST BIRD, RD28]
% lt_figure; hold on;
% ylabel('consolidation');
% title('paired (using exact same syls in both phases');
% X1=[256/282 281/380, 215/227, 390/401, 373.1/450]; % same direction (using same rd23 syls as in 2 dir)
% X2=[124/188, 39/190, 226/297, 168/246, 70/170]; % 
% [~, p]=ttest(X1, X2);
% 
% plot([1.1 2.1], [X1; X2], '-ok');
% 
% lt_plot_bar([1 2], [mean(X1) mean(X2)], {'Errors', [lt_sem(X1) lt_sem(X2)]});
% lt_plot_text(1.5, 1.1*max([X1 X2]), num2str(p, '%3.2g'), 'r')
% xlim([0 3]);
% ylim([0 1.5])
% line(xlim, [1 1], 'LineStyle', '--', 'Color','k')
% 
% 

%% 
lt_subtitle('1: same dir --- 2: diff dir [day 5-6]');

%%
% --- all paired, same dir, using rd23 same syls
x=[253 281, 215, 390]
y=[296 380, 227, 401]; % same direction

[~, p] =ttest(x, y)

% % --- [TEMP] all paired, same dir, using rd23 same syls
% x=[253 281, 215, 390, 373]
% y=[296 380, 227, 401, 450]; % same direction
% 
% [~, p] =ttest(x, y)
% p =signrank(x, y)

% --- all paired, diff dir, using rd23 same syls
x=[124 39 226 168];
y=[188, 190, 297, 246]; % same direction

[~, p] =ttest(x, y)



