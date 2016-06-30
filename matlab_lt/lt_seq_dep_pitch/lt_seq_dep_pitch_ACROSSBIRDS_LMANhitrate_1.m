% samedir
samedir_PBS=[0.45 0.2 0.37 0.22 0.2];
samedir_MUSC=[0.45 0.4 0.37 0.22 0.3];

% bidir
bidir_PBS=[0.5 0.5 0.22 0.35 0.4];
bidir_MUSC=[0.78 0.77 0.6  0.35 0.8];


%% === plot
lt_figure; hold on;
xlim([0 6])
title('NOTE: these are hand entered values!')
ylabel('hit rate');
xlabel('congr(PBS) -- congr(MUSC) -- incon(PBS) -- incon (MUSC)');

% -- samedir
x=[1 2];
plot(x, [samedir_PBS' samedir_MUSC'], 'ob-');

Ymean=mean([samedir_PBS' samedir_MUSC']);
Ysem=lt_sem([samedir_PBS' samedir_MUSC']);
lt_plot_bar(x, Ymean, {'Errors', Ysem, 'Color','b'});

% -- bidir
x=[4 5];
plot(x, [bidir_PBS' bidir_MUSC'], 'or-');

Ymean=mean([bidir_PBS' bidir_MUSC']);
Ysem=lt_sem([bidir_PBS' bidir_MUSC']);
lt_plot_bar(x, Ymean, {'Errors', Ysem, 'Color','r'});

%% === plot (fractions)
lt_figure; hold on;
xlim([0 3])

% -- samedir
x=[1];
Y=samedir_MUSC./samedir_PBS;
plot(x, Y, 'ob');

Ymean=mean(Y);
Ysem=lt_sem(Y);
lt_plot_bar(x, Ymean, {'Errors', Ysem, 'Color','b'});

% -- samedir
x=[2];
Y=bidir_MUSC./bidir_PBS;
plot(x, Y, 'or');

Ymean=mean(Y);
Ysem=lt_sem(Y);
lt_plot_bar(x, Ymean, {'Errors', Ysem, 'Color','r'});


