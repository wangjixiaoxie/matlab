
%% ============= SOME PLOTS
% ================ plot by neuron
lt_figure; hold on;

lt_subplot(2,2,1); hold on;
title('RA');
xlabel('neurnum');
ylabel('Nspks');

inds = strcmp(AllNeurLocation, 'RA');

frmat = FRmatMotifByNeur(:, inds);
plot(1:size(frmat,2), frmat, 'o-');


lt_subplot(2,2,2); hold on;
title('LMAN');
xlabel('neurnum');
ylabel('Nspks');

inds = strcmp(AllNeurLocation, 'LMAN');

frmat = FRmatMotifByNeur(:, inds);
plot(1:size(frmat,2), frmat, 'o-');



% ================ plot by MOTIF
lt_figure; hold on;

lt_subplot(2,2,1); hold on;
title('RA');
xlabel('motifnum');
ylabel('Nspks');

inds = strcmp(AllNeurLocation, 'RA');

frmat = FRmatMotifByNeur(:, inds);
plot(1:size(frmat,1), frmat, 'o-');


lt_subplot(2,2,2); hold on;
title('LMAN');
xlabel('motifnum');
ylabel('Nspks');

inds = strcmp(AllNeurLocation, 'LMAN');

frmat = FRmatMotifByNeur(:, inds);
plot(1:size(frmat,1), frmat, 'o-');



%% ################################# SEPARATE PLOT FOR EACH SYL [BY NEURON]
figcount=1;
subplotrows=6;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
nummotifs = size(FRmatMotifByNeur,1);
hsplots = [];

locthis = 'RA';
inds = strcmp(AllNeurLocation, locthis);
for j=1:nummotifs
   [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([locthis '-' AllMotifRegexp{j}]);
    hsplots = [hsplots hsplot];
    xlabel('neurnum');
    ylabel('Nspks');
    
    fr = FRmatMotifByNeur(j, inds);
    
    plot(fr, '-ok');    
    lt_plot_zeroline;
end
linkaxes(hsplots, 'xy');


%% ################################# SEPARATE PLOT FOR EACH NEURON [function of motif]
figcount=1;
subplotrows=6;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
nummotifs = size(FRmatMotifByNeur,1);
hsplots = [];

locthis = 'LMAN';
inds = find(strcmp(AllNeurLocation, locthis));
for j=inds'
   [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([locthis '-neur' num2str(j)]);
    hsplots = [hsplots hsplot];
    xlabel('motifnum');
    ylabel('Nspks');
    
    fr = FRmatMotifByNeur(:, j);
    
    plot(fr, '-ok');    
    lt_plot_zeroline;
end

% ---------------- OVERLAY ALL AND PLOT MEAN
frmat = FRmatMotifByNeur(:, inds);
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title([locthis '-allneur']);
plot(1:nummotifs, frmat, '-k');
y = mean(frmat,2);
plot(1:nummotifs, y, '-r');
lt_plot_zeroline;

linkaxes(hsplots, 'xy');

set(gca, 'XTick', 1:nummotifs, 'XTickLabel', AllMotifRegexp);
rotateXLabels(gca, 90);
