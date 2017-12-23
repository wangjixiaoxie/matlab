function lt_seq_dep_pitch_ACROSSBIRDS_LMANbasePkLrn(DATSTRUCT, ...
    SeqDepPitch_AcrossBirds_LMAN, PARAMS, plotgeneralization, ...
    targlearnThresh);
%% lt 12/18/17 - predicting magnitude of learning based on baseline AFP effect
% 
% targlearnThresh=40.34; % hz
% plotgeneralization =0; % if 0, then plots hz ()targ dir);

%% Collect learning related things

numexpts = max(DATSTRUCT.singlesyls.ExptCounter);
IndOfPairWithTarg = [];

for i=1:numexpts
   
    sylinds = find(DATSTRUCT.singlesyls.ExptCounter==i);
    indoftarg = sylinds(logical([DATSTRUCT.singlesyls.IsTargAll(sylinds)]));
    assert(length(indoftarg)==1, 'asdfdas');
    
    % ============== go thru all syls.
    for ind = sylinds
       
        % ======== is this targ?
        istarg = DATSTRUCT.singlesyls.IsTargAll(ind);
        sylstr = DATSTRUCT.singlesyls.Syllables{ind};
        
        % ------- skip if is targ
        if istarg==1
            IndOfPairWithTarg = [IndOfPairWithTarg nan];
            continue
        end
        
        % =============== find ind that is this syl with targ syl (paired)
        indofpair = find(DATSTRUCT.pairedsyls.Pairs_OriginalInds(:,1)==min([ind indoftarg]) ...
            & DATSTRUCT.pairedsyls.Pairs_OriginalInds(:,2)==max([ind indoftarg]));
        assert(length(indofpair)==1, 'dafads')
        
    IndOfPairWithTarg = [IndOfPairWithTarg indofpair];
    
    % ========= CORRELATIONS
    
    
    end
end

        DATSTRUCT.singlesyls.IndOfPairWithTarg = IndOfPairWithTarg;

%%
figcount=1;
subplotrows=2;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];


%% Correlation vs. learning (same motif)

% ========= SAME TYPES
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('SAME SYL, SAME MOTIF, PBS');
xlabel('corr (motif)');
ylabel('learning (rel target)');

[indPair, indSing] = intersect(DATSTRUCT.singlesyls.IndOfPairWithTarg, find(DATSTRUCT.pairedsyls.IsSameSyl==1));

% ------------------ RUN
x = DATSTRUCT.pairedsyls.CorrMotif_PBS_pairs(indPair);
if plotgeneralization==1
    ylearn = DATSTRUCT.singlesyls.LearnHzAll(indSing);
ytarg = DATSTRUCT.singlesyls.LearnHzTargAll(indSing);
y = ylearn./ytarg;

x = x(abs(ytarg)>=targlearnThresh);
y = y(abs(ytarg)>=targlearnThresh);

else
y = DATSTRUCT.singlesyls.LearnHzAll(indSing).*DATSTRUCT.singlesyls.TargLearnDir(indSing);
end

lt_regress(y, x, 1);
xlim([-1 1]);
lt_plot_zeroline;
lt_plot_zeroline_vert;


% ========= SAME TYPES
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('SAME SYL, SAME MOTIF, MUSC');
xlabel('corr (motif)');
ylabel('learning (rel target)');

[indPair, indSing] = intersect(DATSTRUCT.singlesyls.IndOfPairWithTarg, find(DATSTRUCT.pairedsyls.IsSameSyl==1));

% ------------------ RUN
x = DATSTRUCT.pairedsyls.CorrMotif_MUSC_pairs(indPair);
if plotgeneralization==1
    ylearn = DATSTRUCT.singlesyls.LearnHzAll(indSing);
ytarg = DATSTRUCT.singlesyls.LearnHzTargAll(indSing);
y = ylearn./ytarg;

x = x(abs(ytarg)>=targlearnThresh);
y = y(abs(ytarg)>=targlearnThresh);

else
y = DATSTRUCT.singlesyls.LearnHzAll(indSing).*DATSTRUCT.singlesyls.TargLearnDir(indSing);
end

lt_regress(y, x, 1);
xlim([-1 1]);
lt_plot_zeroline;
lt_plot_zeroline_vert;

%% Correlation vs. learning (DIFF motif)

% ========= SAME TYPES
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('SAME SYL, DIFF MOTIF, PBS');
xlabel('corr (SONG)');
ylabel('learning (rel target)');

[indPair, indSing] = intersect(DATSTRUCT.singlesyls.IndOfPairWithTarg, ...
    find(DATSTRUCT.pairedsyls.IsSameSyl==1 & DATSTRUCT.pairedsyls.IsSameMotif==0));

% ------------------ RUN
x = DATSTRUCT.pairedsyls.CorrSong_PBS_pairs(indPair);
if plotgeneralization==1
    ylearn = DATSTRUCT.singlesyls.LearnHzAll(indSing);
ytarg = DATSTRUCT.singlesyls.LearnHzTargAll(indSing);
y = ylearn./ytarg;

x = x(abs(ytarg)>=targlearnThresh);
y = y(abs(ytarg)>=targlearnThresh);

else
y = DATSTRUCT.singlesyls.LearnHzAll(indSing).*DATSTRUCT.singlesyls.TargLearnDir(indSing);
end


lt_regress(y, x, 1);
xlim([-1 1]);
lt_plot_zeroline;
lt_plot_zeroline_vert;

% ========= 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('SAME SYL, DIFF MOTIF, MUSC');
xlabel('corr (SONG)');
ylabel('learning (rel target)');

[indPair, indSing] = intersect(DATSTRUCT.singlesyls.IndOfPairWithTarg, ...
    find(DATSTRUCT.pairedsyls.IsSameSyl==1 & DATSTRUCT.pairedsyls.IsSameMotif==0));

% ------------------ RUN
x = DATSTRUCT.pairedsyls.CorrSong_MUSC_pairs(indPair);
if plotgeneralization==1
    ylearn = DATSTRUCT.singlesyls.LearnHzAll(indSing);
ytarg = DATSTRUCT.singlesyls.LearnHzTargAll(indSing);
y = ylearn./ytarg;

x = x(abs(ytarg)>=targlearnThresh);
y = y(abs(ytarg)>=targlearnThresh);

else
y = DATSTRUCT.singlesyls.LearnHzAll(indSing).*DATSTRUCT.singlesyls.TargLearnDir(indSing);
end


lt_regress(y, x, 1);
xlim([-1 1]);
lt_plot_zeroline;
lt_plot_zeroline_vert;



%% Correlation vs. learning (adjacent) 

% =========
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('SAME SYL, ADJACENT, PBS');
xlabel('corr (motif)');
ylabel('learning (rel target)');

[indPair, indSing] = intersect(DATSTRUCT.singlesyls.IndOfPairWithTarg, ...
    find(DATSTRUCT.pairedsyls.IsSameSyl==1 & DATSTRUCT.pairedsyls.NumSylsInBetween==0));

% ------------------ RUN
x = DATSTRUCT.pairedsyls.CorrMotif_PBS_pairs(indPair);
if plotgeneralization==1
    ylearn = DATSTRUCT.singlesyls.LearnHzAll(indSing);
ytarg = DATSTRUCT.singlesyls.LearnHzTargAll(indSing);
y = ylearn./ytarg;

x = x(abs(ytarg)>=targlearnThresh);
y = y(abs(ytarg)>=targlearnThresh);

else
y = DATSTRUCT.singlesyls.LearnHzAll(indSing).*DATSTRUCT.singlesyls.TargLearnDir(indSing);
end


lt_regress(y, x, 1);
xlim([-1 1]);
lt_plot_zeroline;
lt_plot_zeroline_vert;


% ========= 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('SAME SYL, ADJACENT, MUSC');
xlabel('corr (motif)');
ylabel('learning (rel target)');

[indPair, indSing] = intersect(DATSTRUCT.singlesyls.IndOfPairWithTarg, ...
    find(DATSTRUCT.pairedsyls.IsSameSyl==1 & DATSTRUCT.pairedsyls.NumSylsInBetween==0));

% ------------------ RUN
x = DATSTRUCT.pairedsyls.CorrMotif_MUSC_pairs(indPair);
if plotgeneralization==1
    ylearn = DATSTRUCT.singlesyls.LearnHzAll(indSing);
ytarg = DATSTRUCT.singlesyls.LearnHzTargAll(indSing);
y = ylearn./ytarg;

x = x(abs(ytarg)>=targlearnThresh);
y = y(abs(ytarg)>=targlearnThresh);

else
y = DATSTRUCT.singlesyls.LearnHzAll(indSing).*DATSTRUCT.singlesyls.TargLearnDir(indSing);
end


lt_regress(y, x, 1);
xlim([-1 1]);
lt_plot_zeroline;
lt_plot_zeroline_vert;


%% Correlation vs. learning (adjacent) [Diff type]

% =========
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('DIFF SYL, ADJACENT, PBS');
xlabel('corr (motif)');
ylabel('learning (rel target)');

[indPair, indSing] = intersect(DATSTRUCT.singlesyls.IndOfPairWithTarg, ...
    find(DATSTRUCT.pairedsyls.IsSameSyl==0 & DATSTRUCT.pairedsyls.NumSylsInBetween==0));

% ------------------ RUN
x = DATSTRUCT.pairedsyls.CorrMotif_PBS_pairs(indPair);
if plotgeneralization==1
    ylearn = DATSTRUCT.singlesyls.LearnHzAll(indSing);
ytarg = DATSTRUCT.singlesyls.LearnHzTargAll(indSing);
y = ylearn./ytarg;

x = x(abs(ytarg)>=targlearnThresh);
y = y(abs(ytarg)>=targlearnThresh);

else
y = DATSTRUCT.singlesyls.LearnHzAll(indSing).*DATSTRUCT.singlesyls.TargLearnDir(indSing);
end


lt_regress(y, x, 1);
xlim([-1 1]);
lt_plot_zeroline;
lt_plot_zeroline_vert;



% =========
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('DIFF SYL, ADJACENT, MUSC');
xlabel('corr (motif)');
ylabel('learning (rel target)');

[indPair, indSing] = intersect(DATSTRUCT.singlesyls.IndOfPairWithTarg, ...
    find(DATSTRUCT.pairedsyls.IsSameSyl==0 & DATSTRUCT.pairedsyls.NumSylsInBetween==0));

% ------------------ RUN
x = DATSTRUCT.pairedsyls.CorrMotif_MUSC_pairs(indPair);
if plotgeneralization==1
    ylearn = DATSTRUCT.singlesyls.LearnHzAll(indSing);
ytarg = DATSTRUCT.singlesyls.LearnHzTargAll(indSing);
y = ylearn./ytarg;

x = x(abs(ytarg)>=targlearnThresh);
y = y(abs(ytarg)>=targlearnThresh);

else
y = DATSTRUCT.singlesyls.LearnHzAll(indSing).*DATSTRUCT.singlesyls.TargLearnDir(indSing);
end


lt_regress(y, x, 1);
xlim([-1 1]);
lt_plot_zeroline;
lt_plot_zeroline_vert;


%% Correlation vs. learning (adjacent) [alltype]

% =========
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('ALL SYL, ADJACENT, PBS');
xlabel('corr (motif)');
ylabel('learning (rel target)');

[indPair, indSing] = intersect(DATSTRUCT.singlesyls.IndOfPairWithTarg, ...
    find(DATSTRUCT.pairedsyls.NumSylsInBetween==0));

% ------------------ RUN
x = DATSTRUCT.pairedsyls.CorrMotif_PBS_pairs(indPair);
if plotgeneralization==1
    ylearn = DATSTRUCT.singlesyls.LearnHzAll(indSing);
ytarg = DATSTRUCT.singlesyls.LearnHzTargAll(indSing);
y = ylearn./ytarg;

x = x(abs(ytarg)>=targlearnThresh);
y = y(abs(ytarg)>=targlearnThresh);

else
y = DATSTRUCT.singlesyls.LearnHzAll(indSing).*DATSTRUCT.singlesyls.TargLearnDir(indSing);
end


lt_regress(y, x, 1);
xlim([-1 1]);
lt_plot_zeroline;
lt_plot_zeroline_vert;

% =========
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('ALL SYL, ADJACENT, MUSC');
xlabel('corr (motif)');
ylabel('learning (rel target)');

[indPair, indSing] = intersect(DATSTRUCT.singlesyls.IndOfPairWithTarg, ...
    find(DATSTRUCT.pairedsyls.NumSylsInBetween==0));

% ------------------ RUN
x = DATSTRUCT.pairedsyls.CorrMotif_MUSC_pairs(indPair);
if plotgeneralization==1
    ylearn = DATSTRUCT.singlesyls.LearnHzAll(indSing);
ytarg = DATSTRUCT.singlesyls.LearnHzTargAll(indSing);
y = ylearn./ytarg;

x = x(abs(ytarg)>=targlearnThresh);
y = y(abs(ytarg)>=targlearnThresh);

else
y = DATSTRUCT.singlesyls.LearnHzAll(indSing).*DATSTRUCT.singlesyls.TargLearnDir(indSing);
end


lt_regress(y, x, 1);
xlim([-1 1]);
lt_plot_zeroline;
lt_plot_zeroline_vert;


%% Correlation vs. learning [ALL TYPES, SAME MOTIF]

% =========
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('ALL SYL, SAME MOTIF, PBS');
xlabel('corr (motif)');
ylabel('learning (rel target)');

[indPair, indSing] = intersect(DATSTRUCT.singlesyls.IndOfPairWithTarg, ...
    find(DATSTRUCT.pairedsyls.IsSameMotif==1));

% ------------------ RUN
x = DATSTRUCT.pairedsyls.CorrMotif_PBS_pairs(indPair);
if plotgeneralization==1
    ylearn = DATSTRUCT.singlesyls.LearnHzAll(indSing);
ytarg = DATSTRUCT.singlesyls.LearnHzTargAll(indSing);
y = ylearn./ytarg;

x = x(abs(ytarg)>=targlearnThresh);
y = y(abs(ytarg)>=targlearnThresh);

else
y = DATSTRUCT.singlesyls.LearnHzAll(indSing).*DATSTRUCT.singlesyls.TargLearnDir(indSing);
end


lt_regress(y, x, 1);
xlim([-1 1]);
lt_plot_zeroline;
lt_plot_zeroline_vert;


% =========
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('ALL SYL, SAME MOTIF, MUSC');
xlabel('corr (motif)');
ylabel('learning (rel target)');

[indPair, indSing] = intersect(DATSTRUCT.singlesyls.IndOfPairWithTarg, ...
    find(DATSTRUCT.pairedsyls.IsSameMotif==1));

% ------------------ RUN
x = DATSTRUCT.pairedsyls.CorrMotif_MUSC_pairs(indPair);
if plotgeneralization==1
    ylearn = DATSTRUCT.singlesyls.LearnHzAll(indSing);
ytarg = DATSTRUCT.singlesyls.LearnHzTargAll(indSing);
y = ylearn./ytarg;

x = x(abs(ytarg)>=targlearnThresh);
y = y(abs(ytarg)>=targlearnThresh);

else
y = DATSTRUCT.singlesyls.LearnHzAll(indSing).*DATSTRUCT.singlesyls.TargLearnDir(indSing);
end


lt_regress(y, x, 1);
xlim([-1 1]);
lt_plot_zeroline;
lt_plot_zeroline_vert;


%% Correlation vs. learning (NOT adjacent) 

% =========
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('SAME SYL, NOT ADJACENT, PBS');
xlabel('corr (motif)');
ylabel('learning (rel target)');

[indPair, indSing] = intersect(DATSTRUCT.singlesyls.IndOfPairWithTarg, ...
    find(DATSTRUCT.pairedsyls.IsSameSyl==1 & DATSTRUCT.pairedsyls.NumSylsInBetween>0));

% ------------------ RUN
x = DATSTRUCT.pairedsyls.CorrMotif_PBS_pairs(indPair);
if plotgeneralization==1
    ylearn = DATSTRUCT.singlesyls.LearnHzAll(indSing);
ytarg = DATSTRUCT.singlesyls.LearnHzTargAll(indSing);
y = ylearn./ytarg;

x = x(abs(ytarg)>=targlearnThresh);
y = y(abs(ytarg)>=targlearnThresh);

else
y = DATSTRUCT.singlesyls.LearnHzAll(indSing).*DATSTRUCT.singlesyls.TargLearnDir(indSing);
end


lt_regress(y, x, 1);
xlim([-1 1]);
lt_plot_zeroline;
lt_plot_zeroline_vert;


% ========= 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('SAME SYL, NOT ADJACENT, MUSC');
xlabel('corr (motif)');
ylabel('learning (rel target)');

[indPair, indSing] = intersect(DATSTRUCT.singlesyls.IndOfPairWithTarg, ...
    find(DATSTRUCT.pairedsyls.IsSameSyl==1 & DATSTRUCT.pairedsyls.NumSylsInBetween>0));

% ------------------ RUN
x = DATSTRUCT.pairedsyls.CorrMotif_MUSC_pairs(indPair);
if plotgeneralization==1
    ylearn = DATSTRUCT.singlesyls.LearnHzAll(indSing);
ytarg = DATSTRUCT.singlesyls.LearnHzTargAll(indSing);
y = ylearn./ytarg;

x = x(abs(ytarg)>=targlearnThresh);
y = y(abs(ytarg)>=targlearnThresh);

else
y = DATSTRUCT.singlesyls.LearnHzAll(indSing).*DATSTRUCT.singlesyls.TargLearnDir(indSing);
end


lt_regress(y, x, 1);
xlim([-1 1]);
lt_plot_zeroline;
lt_plot_zeroline_vert;

%% ========= IS THERE GREATER CORR FOR SYLS ADJACENT TO TARG?

% ============= SAME TYPES
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('SAME SYL, SAME MOTIF, PBS');
xlabel('dist from targ (syl pos)');
ylabel('corr with targ (motif, PBS)');
plotcol = 'k';
[indPair, indSing] = intersect(DATSTRUCT.singlesyls.IndOfPairWithTarg, find(DATSTRUCT.pairedsyls.IsSameSyl==1));

x = DATSTRUCT.singlesyls.DistFromTarg(indSing);
% x2 = DATSTRUCT.pairedsyls.NumSylsInBetween(indPair);
% assert(all(isnan(x2)==isnan(x)), 'asfasd');
% assert(all(x2(~isnan(x2)) == abs(x(~isnan(x)))-1), 'asdfdas');
y = DATSTRUCT.pairedsyls.CorrMotif_PBS_pairs(indPair);
ymusc = DATSTRUCT.pairedsyls.CorrMotif_MUSC_pairs(indPair);

X = [x'-0.2 x'+0.2];
Y = [y' ymusc'];
plot(X', Y', 'o-', 'Color', plotcol);
plot(X(:,2), Y(:,2), 'o', 'Color', 'r');
lt_plot_zeroline;
lt_plot_zeroline_vert;

[ymean, ysem] = grpstats(Y, x, {'mean', 'sem'});
lt_plot(unique(x(~isnan(x)))-0.4, ymean(:,1), {'Errors', ysem(:,1), 'LineStyle', '-', 'Color', 'k'})
lt_plot(unique(x(~isnan(x)))-0.4, ymean(:,2), {'Errors', ysem(:,2), 'LineStyle', '-', 'Color', 'r'})


% ============= SAME TYPES
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('DIFF SYL, SAME MOTIF, PBS');
xlabel('dist from targ (syl pos)');
ylabel('corr with targ (motif, PBS)');
plotcol = 'k';
[indPair, indSing] = intersect(DATSTRUCT.singlesyls.IndOfPairWithTarg, find(DATSTRUCT.pairedsyls.IsSameSyl==0));

x = DATSTRUCT.singlesyls.DistFromTarg(indSing);
% x2 = DATSTRUCT.pairedsyls.NumSylsInBetween(indPair);
% assert(all(isnan(x2)==isnan(x)), 'asfasd');
% assert(all(x2(~isnan(x2)) == abs(x(~isnan(x)))-1), 'asdfdas');
y = DATSTRUCT.pairedsyls.CorrMotif_PBS_pairs(indPair);
ymusc = DATSTRUCT.pairedsyls.CorrMotif_MUSC_pairs(indPair);

X = [x'-0.2 x'+0.2];
Y = [y' ymusc'];
plot(X', Y', 'o-', 'Color', plotcol);
plot(X(:,2), Y(:,2), 'o', 'Color', 'r');
lt_plot_zeroline;
lt_plot_zeroline_vert;

[ymean, ysem] = grpstats(Y, x, {'mean', 'sem'});
lt_plot(unique(x(~isnan(x)))-0.4, ymean(:,1), {'Errors', ysem(:,1), 'LineStyle', '-', 'Color', 'k'})
lt_plot(unique(x(~isnan(x)))-0.4, ymean(:,2), {'Errors', ysem(:,2), 'LineStyle', '-', 'Color', 'r'})

