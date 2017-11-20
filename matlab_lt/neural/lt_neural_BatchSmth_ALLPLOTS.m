
% ========================= PLOT RAW DAT
dirfield = 'DIR';
mm = 1;
chan = 14;

figcount=1;
subplotrows=5;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

datraw = DATSTRUCT.(dirfield).motifnum(mm).DatAllRaw{chan};
datsm = DATSTRUCT.(dirfield).motifnum(mm).DatAll{chan};
t = DATSTRUCT.(dirfield).motifnum(mm).t;
for k =1 :size(datraw,1)
    motifname = DATSTRUCT.(dirfield).motifnum(mm).motifname;
    
    % --- raw dat    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['[' dirfield ']' motifname '-ch' num2str(chan)]);
    plot(t, datraw(k,:), 'b');
    xlim([pretime-0.05 pretime+0.1]);
    ylim([-200 200])
    line([pretime pretime], ylim);
    lt_plot_zeroline;
    
    % -- smoothed
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['[' dirfield ']' motifname '-ch' num2str(chan)]);
    plot(t, datsm(k, :), 'k', 'LineWidth', 2);
    xlim([pretime-0.05 pretime+0.1]);
    ylim([0 50]);
    line([pretime pretime], ylim);
    lt_plot_zeroline;
end


% ========================= PLOT RAW DAT [compare 2 chans]
dirfield = 'DIR';
mm = 1;
chan1 = 14;
chan2 = 21;

figcount=1;
subplotrows=5;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

% -- chan1
datraw1 = DATSTRUCT.(dirfield).motifnum(mm).DatAllRaw{chan1};
% datsm1 = DATSTRUCT.(dirfield).motifnum(mm).DatAll{chan1};
datraw2 = DATSTRUCT.(dirfield).motifnum(mm).DatAllRaw{chan2};
% datsm2 = DATSTRUCT.(dirfield).motifnum(mm).DatAll{chan2};
t = DATSTRUCT.(dirfield).motifnum(mm).t;
assert(all(size(datraw1) == size(datraw2)), 'diff sample siezes?');
for k =1 :size(datraw1,1)
    motifname = DATSTRUCT.(dirfield).motifnum(mm).motifname;
    
    % --- raw dat (chan1)   
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['[' dirfield ']' motifname '-ch' num2str(chan1)]);
    plot(t, datraw1(k,:), 'b');
    xlim([pretime-0.05 pretime+0.1]);
    ylim([-200 200])
    line([pretime pretime], ylim);
    lt_plot_zeroline;
    
    % -- raw (chan2)
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['[' dirfield ']' motifname '-ch' num2str(chan2)]);
    plot(t, datraw2(k,:), 'b');
    xlim([pretime-0.05 pretime+0.1]);
    ylim([-200 200])
    line([pretime pretime], ylim);
    lt_plot_zeroline;
    
%     % -- smoothed
%     [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%     title(['[' dirfield ']' motifname '-ch' num2str(chan)]);
%     plot(t, datsm(k, :), 'k', 'LineWidth', 2);
%     xlim([pretime-0.05 pretime+0.1]);
%     ylim([0 50]);
%     line([pretime pretime], ylim);
%     lt_plot_zeroline;
end


% ========================= PLOT RAW DAT [compare 3 chans]
dirfield = 'UNDIR';
mm = 1;
chan1 = 9;
chan2 = 14;
chan3 = 21;

figcount=1;
subplotrows=8;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];

% -- chan1
datraw1 = DATSTRUCT.(dirfield).motifnum(mm).DatAllRaw{chan1};
% datsm1 = DATSTRUCT.(dirfield).motifnum(mm).DatAll{chan1};
datraw2 = DATSTRUCT.(dirfield).motifnum(mm).DatAllRaw{chan2};
datraw3 = DATSTRUCT.(dirfield).motifnum(mm).DatAllRaw{chan3};
% datsm2 = DATSTRUCT.(dirfield).motifnum(mm).DatAll{chan2};
t = DATSTRUCT.(dirfield).motifnum(mm).t;
assert(all(size(datraw1) == size(datraw2)), 'diff sample siezes?');
assert(all(size(datraw1) == size(datraw3)), 'diff sample siezes?');
for k =1 :size(datraw1,1)
    motifname = DATSTRUCT.(dirfield).motifnum(mm).motifname;
    
    % --- raw dat (chan1)   
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['[' dirfield ']' motifname '-ch' num2str(chan1)]);
    plot(t, datraw1(k,:), 'b');
    xlim([pretime-0.05 pretime+0.1]);
    ylim([-200 200])
    line([pretime pretime], ylim);
    lt_plot_zeroline;
    
    % -- raw (chan2)
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['[' dirfield ']' motifname '-ch' num2str(chan2)]);
    plot(t, datraw2(k,:), 'b');
    xlim([pretime-0.05 pretime+0.1]);
    ylim([-200 200])
    line([pretime pretime], ylim);
    lt_plot_zeroline;
    
        % --- raw dat (chan3)   
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['[' dirfield ']' motifname '-ch' num2str(chan3)]);
    plot(t, datraw3(k,:), 'b');
    xlim([pretime-0.05 pretime+0.1]);
    ylim([-200 200])
    line([pretime pretime], ylim);
    lt_plot_zeroline;
    

%     % -- smoothed
%     [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%     title(['[' dirfield ']' motifname '-ch' num2str(chan)]);
%     plot(t, datsm(k, :), 'k', 'LineWidth', 2);
%     xlim([pretime-0.05 pretime+0.1]);
%     ylim([0 50]);
%     line([pretime pretime], ylim);
%     lt_plot_zeroline;
end
