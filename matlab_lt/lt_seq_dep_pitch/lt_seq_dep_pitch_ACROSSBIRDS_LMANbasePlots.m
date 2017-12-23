function lt_seq_dep_pitch_ACROSSBIRDS_LMANbasePlots(DATSTRUCT, SeqDepPitch_AcrossBirds_LMAN, ...
    PARAMS, same_type_thr, NullControl_splitday, useHandLabForSametype, ...
    onlyUseFirstExpt, recalcValues)

%% CONVERT FROM STRUCTURE TO VARIABLES (just use name of field)

fnames = fieldnames(DATSTRUCT.singlesyls);
for i=1:length(fnames)
    
   eval([fnames{i} ' = DATSTRUCT.singlesyls.(fnames{i});']);
end

fnames = fieldnames(DATSTRUCT.pairedsyls);
for i=1:length(fnames)
    
   eval([fnames{i} ' = DATSTRUCT.pairedsyls.(fnames{i});']);
end


%% ################################### SINGLE SYL PLOTS ##########
%% ###############################################################

%% +++++++++++++ PLOTS

figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];



%% ====== DISTRIBUTION OF PBS AND MUSC PITCH
%  ########### 1) distrubitions
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('distr of baseline PBS pitch');

% --- PBS
[~, Xcenters] = lt_plot_histogram(Pitch_PBSraw, '', 1, 1, 0, 1, 'k'); % PBS

% --- MUSC
lt_plot_histogram(Pitch_PBSraw+Pitch_MUSC_minus_PBS, Xcenters, 1, 1, 0, 1, 'r'); % MUSC

% ############## 2) paired
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('all syls');

x = [1 2];
y = [Pitch_PBSraw', Pitch_PBSraw'+Pitch_MUSC_minus_PBS'];
plot(x, y, '-ok');
plot(2, y(:,2), 'or');
xlim([0 3]);
xlabel('PBS -- MUSC');


%% =========== DISTRIBUTIONS, SEPARATE FOR EACH EXPERIMENT
figcount=1;
subplotrows=4;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];


numexpts = max(ExptCounter);
for i=1:numexpts
    
    % ############## 2) paired
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    inds = ExptCounter==i;
    
    bname = unique({Birdnames{inds}});
    ename = unique({Exptnames{inds}});
    assert(length([bname ename])==2, 'asfds');
    title([bname '-' ename]);
    
    x = [1 2];
    y = [Pitch_PBSraw(inds)', Pitch_PBSraw(inds)'+Pitch_MUSC_minus_PBS(inds)'];
    sylnames = {Syllables{inds}};
    
    plot(x, y, '-ok');
    plot(2, y(:,2), 'or');
    
    for j=1:length(y)
        lt_plot_text(2.1, y(j,2), sylnames{j}, 'r');
    end
    xlim([0 3]);
    xlabel('PBS -- MUSC');
    
end

%% ====== 1) distribution of shift caused by muscimol
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('distr of baseline musc shifts');

lt_plot_histogram(Pitch_MUSC_minus_PBS, '', 1, 1, 0, 1, 'k');
p=signrank(Pitch_MUSC_minus_PBS);
lt_plot_pvalue(p, 'signrank')

% --- use zscore
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('distr of baseline musc shifts');
ylabel('zscore of PBS across rends/days');
lt_plot_histogram(Pitch_MUSC_minus_PBS_zscore, '', 1, 1, 0, 1, 'k');
p=signrank(Pitch_MUSC_minus_PBS_zscore);
lt_plot_pvalue(p, 'signrank')


% --- does shift correlate with baseline std?
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('musc shifts vs. baseline std (across rends and days)');

plot(PitchSTD_PBS, Pitch_MUSC_minus_PBS, 'ok');
lt_regress(Pitch_MUSC_minus_PBS, PitchSTD_PBS, 1, 0, 1, 1, 'r', 0);

% --- does shift correlate with baseline std?
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('abs(musc shifts) vs. baseline std (across rends and days)');

plot(PitchSTD_PBS, abs(Pitch_MUSC_minus_PBS), 'ok');
lt_regress(abs(Pitch_MUSC_minus_PBS), PitchSTD_PBS, 1, 0, 1, 1, 'r', 0);





%% ############################# PAIRED SYL PLOTS ################
%% ###############################################################
%% ====== PLOT (scatter of pitch shifts due to musc)


figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];


% === plot all pairwise effect (45 degree)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[all pairs]');
color='k';
xlabel('syl1 shift (MUSC - PBS)')
ylabel('syl2 shift (MUSC - PBS)')
inds=find(AcousticDist_pairs>0);
X=[];
Y=[];
for j=inds;
    x=PitchShiftBothSyls_paired{j}(1);
    y=PitchShiftBothSyls_paired{j}(2);
    X=[X x];
    Y=[Y y];
end
% -- randomly switch positions of X and Y (since order doesn't matter)
XY = [X'  Y'];
indflip =rand(size(XY,1), 1)>0.5;
XY(indflip,:) = fliplr(XY(indflip,:));
lt_regress(XY(:,2), XY(:,1), 1, 0, 1, 1, color);
lt_plot_zeroline;lt_plot_zeroline_vert;

% === plot all pairwise effect (45 degree)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same');
color='b';
xlabel('syl1 shift (MUSC - PBS)')
ylabel('syl2 shift (MUSC - PBS)')
inds=find(AcousticDist_pairs<same_type_thr);
X=[];
Y=[];
for j=inds;
    x=PitchShiftBothSyls_paired{j}(1);
    y=PitchShiftBothSyls_paired{j}(2);
    X=[X x];
    Y=[Y y];
end
XY = [X'  Y'];
indflip =rand(size(XY,1), 1)>0.5;
XY(indflip,:) = fliplr(XY(indflip,:));
lt_regress(XY(:,2), XY(:,1), 1, 0, 1, 1, color);
lt_plot_zeroline;lt_plot_zeroline_vert;

% === plot all pairwise effect (45 degree)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('diff');
color='r';
xlabel('syl1 shift (MUSC - PBS)')
ylabel('syl2 shift (MUSC - PBS)')
inds=find(AcousticDist_pairs>=same_type_thr);
X=[];
Y=[];
for j=inds;
    x=PitchShiftBothSyls_paired{j}(1);
    y=PitchShiftBothSyls_paired{j}(2);
    X=[X x];
    Y=[Y y];
end
XY = [X'  Y'];
indflip =rand(size(XY,1), 1)>0.5;
XY(indflip,:) = fliplr(XY(indflip,:));
lt_regress(XY(:,2), XY(:,1), 1, 0, 1, 1, color);
lt_plot_zeroline;lt_plot_zeroline_vert;


% ################################################################################
% ######################################## same, but using hand label same type
% === plot all pairwise effect (45 degree)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same (handlab)');
color='b';
xlabel('syl1 shift (MUSC - PBS)')
ylabel('syl2 shift (MUSC - PBS)')
inds=find(IsSameSyl);
X=[];
Y=[];
for j=inds;
    x=PitchShiftBothSyls_paired{j}(1);
    y=PitchShiftBothSyls_paired{j}(2);
    X=[X x];
    Y=[Y y];
end
XY = [X'  Y'];
indflip =rand(size(XY,1), 1)>0.5;
XY(indflip,:) = fliplr(XY(indflip,:));
lt_regress(XY(:,2), XY(:,1), 1, 0, 1, 1, color);
lt_plot_zeroline;lt_plot_zeroline_vert;


% === plot all pairwise effect (45 degree)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('diff');
color='r';
xlabel('syl1 shift (MUSC - PBS)')
ylabel('syl2 shift (MUSC - PBS)')
inds=find(IsSameSyl==0);
X=[];
Y=[];
for j=inds;
    x=PitchShiftBothSyls_paired{j}(1);
    y=PitchShiftBothSyls_paired{j}(2);
    X=[X x];
    Y=[Y y];
end
XY = [X'  Y'];
indflip =rand(size(XY,1), 1)>0.5;
XY(indflip,:) = fliplr(XY(indflip,:));
lt_regress(XY(:,2), XY(:,1), 1, 0, 1, 1, color);
lt_plot_zeroline;lt_plot_zeroline_vert;


%% ################### SCATTER (SYL1 = starting higher pitch])

figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

% ######################################## same, but using hand label same type
% === plot all pairwise effect (45 degree) 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('same [syl1: higher hz]');
color='b';
xlabel('syl1 shift (MUSC - PBS)')
ylabel('syl2 shift (MUSC - PBS)')

inds=find(IsSameSyl);
XY = cell2mat(PitchShiftBothSyls_paired(inds)');

% -- sort so that syl1 is the one with higher hz
ffsyl1 = Pitch_PBSraw(Pairs_OriginalInds(inds, 1));
ffsyl2 = Pitch_PBSraw(Pairs_OriginalInds(inds, 2));
% -- for any cases where ffsyl2 > ffsyl1, flip them
indflip = ffsyl2>ffsyl1;
XY(indflip,:) = fliplr(XY(indflip,:));

% -- annotate cases where MUSC brought separation closer
indsMuscCloser = abs(separationMUSC_paired(inds)) < abs(separationPBS_paired(inds));

% --- PLOT
lt_regress(XY(:,2), XY(:,1), 1, 0, 1, 1, color);
plot(XY(indsMuscCloser, 1), XY(indsMuscCloser,2), 'mo');
lt_plot_zeroline;lt_plot_zeroline_vert;
line([-50 50], [-50 50], 'Color','k');

% --- also plot how much greater pitch of syl 2 is (rel syl 1)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('same [syl1: higher hz]');
color='b';
xlabel('syl1 shift (MUSC - PBS)')
ylabel('syl1 hz - syl2 hz (base)')

Ysep= abs(ffsyl2-ffsyl1);
plot(XY(:,1), Ysep, 'ok');

% ----------------------------
linkaxes(hsplots, 'x');

% ######################################## diff
% === plot all pairwise effect (45 degree) 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('diff [syl1: higher hz]');
color='r';
xlabel('syl1 shift (MUSC - PBS)')
ylabel('syl2 shift (MUSC - PBS)')

inds=find(IsSameSyl==0);
XY = cell2mat(PitchShiftBothSyls_paired(inds)');

% -- sort so that syl1 is the one with higher hz
ffsyl1 = Pitch_PBSraw(Pairs_OriginalInds(inds, 1));
ffsyl2 = Pitch_PBSraw(Pairs_OriginalInds(inds, 2));
% -- for any cases where ffsyl2 > ffsyl1, flip them
indflip = ffsyl2>ffsyl1;
XY(indflip,:) = fliplr(XY(indflip,:));

% -- annotate cases where MUSC brought separation closer
indsMuscCloser = abs(separationMUSC_paired(inds)) < abs(separationPBS_paired(inds));

% --- PLOT
lt_regress(XY(:,2), XY(:,1), 1, 0, 1, 1, color);
plot(XY(indsMuscCloser, 1), XY(indsMuscCloser,2), 'mo');
lt_plot_zeroline;lt_plot_zeroline_vert;
line([-50 50], [-50 50], 'Color','k');
xlim([-55 55]); ylim([-55 55]);

% --- also plot how much greater pitch of syl 2 is (rel syl 1)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('diff [syl1: higher hz]');
color='b';
xlabel('syl1 shift (MUSC - PBS)')
ylabel('syl1 hz - syl2 hz (base)')

Ysep= abs(ffsyl2-ffsyl1);
plot(XY(:,1), Ysep, 'ok');

% ----------------------------
linkaxes(hsplots, 'x');


%% ##### MORE LIKELY FOR AFP BIAS TO BE IN OPPOSITE DIRECTIONS IF START CLOSER?

% =========================== SAME
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('SAME');
xlabel('starting sepration (abs, hz)')
ylabel('diff in AFP bias (abs, hz)')

inds = IsSameSyl==1;

x = abs(separationPBS_paired(inds));
XY = cell2mat(PitchShiftBothSyls_paired(inds)');
DiffInAFPbias = abs(XY(:,2)-XY(:,1));

plot(x, DiffInAFPbias, 'ok');
lt_regress(DiffInAFPbias, x, 1)
% 
% y = abs(separationMUSC_paired(inds));
% plot(x, y, 'or');
% 
% 

%% ######## MORE OVERLAP IN PITCH DISTRIBUTION DURING MUSC?




%% ####################### compare distributions to shuffle distributions
% to shuffle: 1) null hypothesis is that AFP bias is independent of pair
% (i.e. each syl in a pair gets random draw from entire distribution of
% biases (across all datapoints across contexts and syls)

%% #################### DISTRIBUTION OF TYPES OF PAIRED EFFECTS
CategoriesStruct = struct;

% ====== SAME SYL
inds=find(IsSameSyl);

CategoriesStruct.same = fn_getAllShiftDistributions(inds, PitchShiftBothSyls_paired, ...
    separationPBS_paired, separationMUSC_paired, Pairs_OriginalInds, Pitch_PBSraw, 1);

% ====== DIFF SYL
inds=find(IsSameSyl==0);

CategoriesStruct.diff = fn_getAllShiftDistributions(inds, PitchShiftBothSyls_paired, ...
    separationPBS_paired, separationMUSC_paired, Pairs_OriginalInds, Pitch_PBSraw, 1);


% ######################## COMPARE SAME AND DIFF
goodclasses = [1 2 3]; % those that bring closer
% === the proportion of classes that brings closer vs. further
sylfield = 'same';
Ngoodclass = sum(CategoriesStruct.(sylfield).NAll(goodclasses));
Nbadclass = sum(CategoriesStruct.(sylfield).NAll) - Ngoodclass;
disp(['*** ' sylfield ': ' num2str(Ngoodclass) '/' num2str(Ngoodclass+Nbadclass)]);

sylfield = 'diff';
Ngoodclass = sum(CategoriesStruct.(sylfield).NAll(goodclasses));
Nbadclass = sum(CategoriesStruct.(sylfield).NAll) - Ngoodclass;
disp(['*** ' sylfield ': ' num2str(Ngoodclass) '/' num2str(Ngoodclass+Nbadclass)]);



%% SHUFFLE PROCEDURE (COMBINED SAME AND DIFF TO TEST DIFFERENCES BETWEEN THEM)
% ################################################ SAME/DIFF COMBINED
NAll_same = [];
NAll_diff = [];
SepPBS_same = [];
SepMUSC_same = [];

SepPBS_diff = [];
SepMUSC_diff = [];

Ncycles = 1000;
listofbirds = unique(DATSTRUCT.singlesyls.Birdnames);
for nn = 1:Ncycles
    
    % --- shuffle all biases (maintaining everything else)
    % NOTE: shuffle at level of syl, not at level of paired syls, since
    % a given syl contributes to multiple pairs
    
    if (1)
    % ============ DO THIS SEPARATELY FOR EACH BIRD [correct]
    indshuff = 1:length(Pitch_MUSC_minus_PBS);
        for bname = listofbirds
            idxtmp = find(strcmp(DATSTRUCT.singlesyls.Birdnames, bname)); % data for this bird
            
            indshuff(idxtmp) = indshuff(idxtmp(randperm(length(idxtmp)))); % shuffle just inds for this bird
            
        end
        PitchShifts_shuff = Pitch_MUSC_minus_PBS(indshuff);
    else
        % === pool across birds
        indshuff = randperm(length(Pitch_MUSC_minus_PBS));
        PitchShifts_shuff = Pitch_MUSC_minus_PBS(indshuff);
    end
    % ============================== recalculate Pitch shifts paired
    PitchShiftBothSyls_SHUFF = {};
    NumExpts=max(ExptCounter);
    
    for i=1:NumExpts
        Inds=find(ExptCounter==i);
        % ---- get all pairwise stuff
        for j=1:length(Inds)
            ind1=Inds(j);
            for jj=j+1:length(Inds)
                ind2=Inds(jj);
                % ==== collect pairwise stats
                PitchShiftBothSyls_SHUFF=[PitchShiftBothSyls_SHUFF [PitchShifts_shuff(ind1) PitchShifts_shuff(ind2)]];
            end
        end
    end
    
    % =================================== get distribution of classes
    % ----------- SAME
    inds = IsSameSyl==1;
    structtmp = fn_getAllShiftDistributions(inds, PitchShiftBothSyls_SHUFF, ...
    separationPBS_paired, separationMUSC_paired, Pairs_OriginalInds, Pitch_PBSraw, 0);
    
    NAll_same = [NAll_same; structtmp.NAll'];
    SepPBS_same = [SepPBS_same; structtmp.SeparationPBS'];
    SepMUSC_same = [SepMUSC_same; structtmp.SeparationMUSC'];

    % ----------- DIFF
    inds = IsSameSyl==0;
    structtmp = fn_getAllShiftDistributions(inds, PitchShiftBothSyls_SHUFF, ...
    separationPBS_paired, separationMUSC_paired, Pairs_OriginalInds, Pitch_PBSraw, 0);
    
    NAll_diff = [NAll_diff; structtmp.NAll'];
    SepPBS_diff = [SepPBS_diff; structtmp.SeparationPBS'];
    SepMUSC_diff = [SepMUSC_diff; structtmp.SeparationMUSC'];

   
end
CategoriesStruct.same_SHUFF.NAll_All = NAll_same;
CategoriesStruct.same_SHUFF.SepPBS = SepPBS_same;
CategoriesStruct.same_SHUFF.SepMUSC= SepMUSC_same;

CategoriesStruct.diff_SHUFF.NAll_All = NAll_diff;
CategoriesStruct.diff_SHUFF.SepPBS= SepPBS_diff;
CategoriesStruct.diff_SHUFF.SepMUSC= SepMUSC_diff;




%% ############################################ PLOT SHUFFLE RESULTS
    lt_figure; hold on;

% ====== COMPARE EACH POSITION (PROPORTION)
lt_subplot(2,1,1); hold on;

nsame = CategoriesStruct.same.NAll./sum(CategoriesStruct.same.NAll);
ndiff = CategoriesStruct.diff.NAll./sum(CategoriesStruct.diff.NAll);
ndat = nsame - ndiff;

nsame_null = CategoriesStruct.same_SHUFF.NAll_All./sum(CategoriesStruct.same.NAll);
ndiff_null = CategoriesStruct.diff_SHUFF.NAll_All./sum(CategoriesStruct.diff.NAll);
nnull = nsame_null - ndiff_null;

% - plot null
X = 1:length(ndat);
plot(X, nnull, '-r');

% - plot dat
lt_plot_bar(X, ndat);

% - calculate p
for j=1:length(ndat)
    
%     p = sum(nnull(:,j)>ndat(j))/size(nnull,1);
    p = sum(nnull(:,j)>ndat(j))/size(nnull,1);
    lt_plot_text(j-0.3, 1.2*ndat(j), ['p=' num2str(p)], 'b');
    lt_plot_annotation(1, 'one-sided', 'b');
    
end

% ==================== COMPARE "GOOD" CLASSES
nsame = sum(CategoriesStruct.same.NAll(goodclasses)./sum(CategoriesStruct.same.NAll));
ndiff = sum(CategoriesStruct.diff.NAll(goodclasses)./sum(CategoriesStruct.diff.NAll));
ndat = nsame - ndiff;

nsame_null = sum(CategoriesStruct.same_SHUFF.NAll_All(:,goodclasses)./sum(CategoriesStruct.same.NAll),2);
ndiff_null = sum(CategoriesStruct.diff_SHUFF.NAll_All(:,goodclasses)./sum(CategoriesStruct.diff.NAll), 2);
nnull = nsame_null - ndiff_null;

lt_subplot(2,2,3); hold on;
title('diff (good classes');
plot(1, nnull, 'o', 'Color', 'r');
lt_plot_bar(1, ndat);

% -- plot p
    p = sum(nnull>ndat)/size(nnull,1);
    lt_plot_text(1, 1.2*ndat, ['p=' num2str(p)], 'b');
    lt_plot_annotation(1, 'one-sided', 'b');




% ######################################## plot null distribution
lt_figure; hold on;

% ============================================ same
sylfield = 'same';
lt_subplot(3,1,1); hold on;

% ----------- RUN
title(sylfield);
NAll_All = CategoriesStruct.([sylfield '_SHUFF']).NAll_All;

 % --- shuffle
X = 1:length(CategoriesStruct.(sylfield).NAll);
plot(X, NAll_All', '-r');

% -- actual dat
lt_plot_bar(X, CategoriesStruct.(sylfield).NAll, {'Color', [0.6 0.6 0.6]});

% -- get p value for each categorie
for j=1:length(CategoriesStruct.(sylfield).NAll)
    nthis = CategoriesStruct.(sylfield).NAll(j);
    nnull = NAll_All(:,j);
    p = sum(nnull>nthis)/length(nnull);
    lt_plot_text(j-0.3, 1.2*nthis , ['p=' num2str(p)], 'b');
end


% ============================================ same
sylfield = 'diff';
lt_subplot(3,1,2); hold on;

% ----------- RUN
title(sylfield);
NAll_All = CategoriesStruct.([sylfield '_SHUFF']).NAll_All;

 % --- shuffle
X = 1:length(CategoriesStruct.(sylfield).NAll);
plot(X, NAll_All', '-r');

% -- actual dat
lt_plot_bar(X, CategoriesStruct.(sylfield).NAll, {'Color', [0.6 0.6 0.6]});

% -- get p value for each categorie
for j=1:length(CategoriesStruct.(sylfield).NAll)
    nthis = CategoriesStruct.(sylfield).NAll(j);
    nnull = NAll_All(:,j);
    p = sum(nnull>nthis)/length(nnull);
    lt_plot_text(j-0.3, 1.2*nthis , ['p=' num2str(p)], 'b');
end


% ====================================== [proportion in good class]
% [frequency]
lt_subplot(3,2,5); hold on;
xlabel('same -- diff');
ylabel('count');
% --------------------------------------- SAME
sylfield = 'same';
X = 1;

ndat = CategoriesStruct.(sylfield).NAll(goodclasses);
ndat = sum(ndat);
nnull = sum(CategoriesStruct.([sylfield '_SHUFF']).NAll_All(:, goodclasses),2);
% -- plot null
plot(X, nnull, 'o', 'Color', [0.6 0.6 0.6])
% -- plot dat
lt_plot_bar(X, ndat);
% -- p value
p = sum(nnull>ndat)/length(nnull);
lt_plot_text(X, 1.2*ndat, ['p=' num2str(p)], 'r');


% --------------------------------------- DIFF
sylfield = 'diff';
X = 2;

ndat = CategoriesStruct.(sylfield).NAll(goodclasses);
ndat = sum(ndat);
nnull = sum(CategoriesStruct.([sylfield '_SHUFF']).NAll_All(:, goodclasses),2);
% -- plot null
plot(X, nnull, 'o', 'Color', [0.6 0.6 0.6])
% -- plot dat
lt_plot_bar(X, ndat);
% -- p value
p = sum(nnull>ndat)/length(nnull);
lt_plot_text(X, 1.2*ndat, ['p=' num2str(p)], 'r');

lt_subplot(3,2,5); hold on;
xlabel('same -- diff');
ylabel('count');


% ====================================== [proportion in good class]
% [PROPORTION]
lt_subplot(3,2,6); hold on;
xlabel('same -- diff');
ylabel('fraction');

% ------------------------ SAME
sylfield = 'same';
X = 1;

ndat = sum(CategoriesStruct.(sylfield).NAll(goodclasses));
ndat_tot = sum(CategoriesStruct.(sylfield).NAll);
ndat = ndat/ndat_tot;

nnull = sum(CategoriesStruct.([sylfield '_SHUFF']).NAll_All(:, goodclasses),2);
nnull_tot =  sum(CategoriesStruct.([sylfield '_SHUFF']).NAll_All,2);
nnull = nnull./nnull_tot;

% -- plot null
plot(X, nnull, 'o', 'Color', [0.6 0.6 0.6])
% -- plot dat
lt_plot_bar(X, ndat);
% -- p value
p = sum(nnull>ndat)/length(nnull);
lt_plot_text(X, 1.2*ndat, ['p=' num2str(p)], 'r');


% ------------------------ SAME
sylfield = 'diff';
X = 2;

ndat = sum(CategoriesStruct.(sylfield).NAll(goodclasses));
ndat_tot = sum(CategoriesStruct.(sylfield).NAll);
ndat = ndat/ndat_tot;

nnull = sum(CategoriesStruct.([sylfield '_SHUFF']).NAll_All(:, goodclasses),2);
nnull_tot =  sum(CategoriesStruct.([sylfield '_SHUFF']).NAll_All,2);
nnull = nnull./nnull_tot;

% -- plot null
plot(X, nnull, 'o', 'Color', [0.6 0.6 0.6])
% -- plot dat
lt_plot_bar(X, ndat);
% -- p value
p = sum(nnull>ndat)/length(nnull);
lt_plot_text(X, 1.2*ndat, ['p=' num2str(p)], 'r');





    
    
%% ============= [SHUFFLE] - COMPARE CHANGE IN SEPARATION
lt_figure; hold on;
title('separation');

Y = {};
% ============= actual dat
sepChange_same = cell2mat(CategoriesStruct.same.SeparationMUSC') - ...
    cell2mat(CategoriesStruct.same.SeparationPBS');
Y{1} = sepChange_same;

sepChange_diff = cell2mat(CategoriesStruct.diff.SeparationMUSC') - ...
    cell2mat(CategoriesStruct.diff.SeparationPBS');
Y{2} = sepChange_diff;

ymean = [mean(Y{1}) mean(Y{2})];
ysem = [lt_sem(Y{1}) lt_sem(Y{2})];
% lt_plot_MultDist(Y, [1 2]);
lt_plot([1 2], ymean, {'Errors', ysem, 'Color', 'k'});

% --- overlay null data
numshuffs = size(CategoriesStruct.same_SHUFF.SepPBS,1);
Yshuff = [];
for j=1:numshuffs
sepchange_same_shuff = cell2mat(CategoriesStruct.same_SHUFF.SepMUSC(j,:)) - ...
    cell2mat(CategoriesStruct.same_SHUFF.SepPBS(j,:));

sepchange_diff_shuff = cell2mat(CategoriesStruct.diff_SHUFF.SepMUSC(j,:)) - ...
    cell2mat(CategoriesStruct.diff_SHUFF.SepPBS(j,:));

y = [mean(sepchange_same_shuff) mean(sepchange_diff_shuff)];
plot([1.1 1.9], y, '-', 'Color', [0.7 0.7 0.7]);

Yshuff = [Yshuff; y];
end

% ========== significance
% -- same
p = sum(Yshuff(:,1) > ymean(1))/size(Yshuff,1);
lt_plot_text(1, 1.2*max(Yshuff(:,1)), ['p=' num2str(p)], 'r');

% -- diff
p = sum(Yshuff(:,2) > ymean(2))/size(Yshuff,1);
lt_plot_text(2, 1.2*max(Yshuff(:,2)), ['p=' num2str(p)], 'r');


% -- same minus diff
p = sum((Yshuff(:,1) - Yshuff(:,2)) > (ymean(1) - ymean(2)))/size(Yshuff,1);
lt_plot_text(1.5, 1.2*max(Yshuff(:,2)), ['(vs) p=' num2str(p)], 'r');

    
    
    
%% ################### [SYL1 = starting higher pitch] [change as fraction of separation]
if (0)
figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

% ######################################## same, but using hand label same type
% === plot all pairwise effect (45 degree) 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('same [syl1: higher hz]');
color='b';
xlabel('syl1 shift (MUSC - PBS) [frac of separation]')
ylabel('syl2 shift (MUSC - PBS) [frac of seapration]')

inds=find(IsSameSyl);
XY = cell2mat(PitchShiftBothSyls_paired(inds)');

% -- sort so that syl1 is the one with higher hz
ffsyl1 = Pitch_PBSraw(Pairs_OriginalInds(inds, 1));
ffsyl2 = Pitch_PBSraw(Pairs_OriginalInds(inds, 2));
% -- for any cases where ffsyl2 > ffsyl1, flip them
indflip = ffsyl2>ffsyl1;
XY(indflip,:) = fliplr(XY(indflip,:));

% -- convert pitch shift to fraction of separation
Ysep= abs(ffsyl2-ffsyl1);
Ysep = repmat(Ysep', 1, 2);
XY = XY./Ysep;

% --- PLOT
lt_regress(XY(:,2), XY(:,1), 1, 0, 1, 1, color);
lt_plot_zeroline;lt_plot_zeroline_vert;
line([-50 50], [-50 50], 'Color','k');

% --- also plot how much greater pitch of syl 2 is (rel syl 1)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('same [syl1: higher hz]');
color='b';
xlabel('syl1 shift (MUSC - PBS)')
ylabel('syl1 hz - syl2 hz (base)')

Ysep= abs(ffsyl2-ffsyl1);
plot(XY(:,1), Ysep, 'ok');

% ----------------------------
linkaxes(hsplots, 'x');

end
%% OTHER PLOTS RELATIING MUSC EFFECT TO ACOUSTIC DISTANCE

% === 1) diff in effect of musc vs. acoustic distance
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[all pairs] diff in effect of musc vs. acoustic distance');
xlabel('acoustic dist')
ylabel('abs diff in effect of musc (hz)');

lt_regress(DiffInPitchShiftCausedByMusc_pairs, AcousticDist_pairs, 1, 0, 1, 1, 'k');


% === 2) similarity of sign of musc effect vs. acoustic distance
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[all pairs] similarity of sign of musc effect vs. acoustic distance');
xlabel('acoustic dist')
ylabel('similarity of sign of musc effect');

lt_regress(similarDirOfShift_pairs, AcousticDist_pairs, 1, 0, 1, 1, 'k');


% === 2) similarity of sign of musc effect vs. acoustic distance
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[all pairs] similarity of sign of musc effect vs. acoustic distance');
xlabel('acoustic dist')
ylabel('similarity of sign of musc effect');

[~, inds]=sort(AcousticDist_pairs);
X=AcousticDist_pairs(inds);
Y=similarDirOfShift_pairs(inds);

plot(X, Y, 'ok');
Yrun=lt_running_stats(Y, 50);
Xrun=lt_running_stats(X, 50);

shadedErrorBar(Xrun.Mean, Yrun.Mean, Yrun.SEM, {'Color','b'}, 1);




%% distributions of changes in separation

% ---
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[same-type pairs] ');
xlabel('change in abs seperation during MUSC vs. pbs')

inds=AcousticDist_pairs<same_type_thr;
Y=changeInAbsSeparationDuringMUSC_pairs(inds);

lt_plot_histogram(Y, '', 1, 1, 0, 1, 'k');
p=signrank(Y);
lt_plot_pvalue(p, 'signrank')


% ---
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[diff-type pairs] ');
xlabel('change in abs seperation during MUSC vs. pbs')

inds=AcousticDist_pairs>=same_type_thr;
Y=changeInAbsSeparationDuringMUSC_pairs(inds);

lt_plot_histogram(Y, '', 1, 1, 0, 1, 'k');
p=signrank(Y);
lt_plot_pvalue(p, 'signrank')


% ####################################### REPEAT, BUT USING HAND LAB
% ---
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[same-type pairs [HANDLAB] ');
xlabel('change in abs seperation during MUSC vs. pbs')

inds=IsSameSyl;
Y=changeInAbsSeparationDuringMUSC_pairs(inds);

lt_plot_histogram(Y, '', 1, 1, 0, 1, 'k');
p=signrank(Y);
lt_plot_pvalue(p, 'signrank')

% ---
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[diff-type pairs [HAND LAB] ');
xlabel('change in abs seperation during MUSC vs. pbs')

inds=IsSameSyl==0;
Y=changeInAbsSeparationDuringMUSC_pairs(inds);

lt_plot_histogram(Y, '', 1, 1, 0, 1, 'k');
p=signrank(Y);
lt_plot_pvalue(p, 'signrank')


%% COMPARE SEPARATION POST VS PRE MUSC

% ---
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[all pairs] ');
xlabel('change in abs seperation during MUSC vs. pbs')

Y=changeInAbsSeparationDuringMUSC_pairs;

lt_plot_histogram(Y, '', 1, 1, 0, 1, 'k');
p=signrank(Y);
lt_plot_pvalue(p, 'signrank')


% --- SAME TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[same-type] ');
xlabel('separation (PBS)')
ylabel('separation (MUSC)');

inds=AcousticDist_pairs<same_type_thr;

X=separationPBS_paired(inds);
Y=separationMUSC_paired(inds);

lt_plot_45degScatter(X, Y, 'b');

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[same-type] ');
xlabel('abs separation (PBS)')
ylabel('abs separation (MUSC)');

inds=AcousticDist_pairs<same_type_thr;

X=abs(separationPBS_paired(inds));
Y=abs(separationMUSC_paired(inds));

lt_plot_45degScatter(X, Y, 'b');


% --- DIFF TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[diff-type] ');
xlabel('separation (PBS)')
ylabel('separation (MUSC)');

inds=AcousticDist_pairs>=same_type_thr;

X=separationPBS_paired(inds);
Y=separationMUSC_paired(inds);

lt_plot_45degScatter(X, Y, 'r');

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[diff-type] ');
xlabel('abs separation (PBS)')
ylabel('abs separation (MUSC)');

inds=AcousticDist_pairs>=same_type_thr;

X=abs(separationPBS_paired(inds));
Y=abs(separationMUSC_paired(inds));

lt_plot_45degScatter(X, Y, 'r');



% ---
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[all pairs] ');
xlabel('acoustic dist')
ylabel('change in separation of raw pitch during musc vs pbs');

lt_regress(changeInAbsSeparationDuringMUSC_pairs, AcousticDist_pairs, 1, 0, 1, 1, 'k');


% ################################################### REPLOT, BUT USING
% HAND LABELED
% --- SAME TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[same-type HANDLAB] ');
xlabel('separation (PBS)')
ylabel('separation (MUSC)');

inds=IsSameSyl;

X=separationPBS_paired(inds);
Y=separationMUSC_paired(inds);

lt_plot_45degScatter(X, Y, 'b');

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[same-type HANDLAB] ');
xlabel('abs separation (PBS)')
ylabel('abs separation (MUSC)');

X=abs(separationPBS_paired(inds));
Y=abs(separationMUSC_paired(inds));

lt_plot_45degScatter(X, Y, 'b');


% --- DIFF TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[diff-type HANDLAB] ');
xlabel('separation (PBS)')
ylabel('separation (MUSC)');

inds=IsSameSyl==0;

X=separationPBS_paired(inds);
Y=separationMUSC_paired(inds);

lt_plot_45degScatter(X, Y, 'r');


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[diff-type HANDLAB] ');
xlabel('abs separation (PBS)')
ylabel('abs separation (MUSC)');

inds=IsSameSyl==0;

X=abs(separationPBS_paired(inds));
Y=abs(separationMUSC_paired(inds));

lt_plot_45degScatter(X, Y, 'r');

%% ====================== SEPRATAION (MUSC VS. PBS) DEP ON CONTEXTUAL SIMILIARTIY
lt_figure; hold on;

% ================= 0 shared
lt_subplot(3,2,1); hold on;
title('[same-type [0 shared syl] ');
xlabel('abs separation (PBS)')
ylabel('abs separation (MUSC)');

inds = IsSameSyl & DATSTRUCT.pairedsyls.NumSylsSharedInContext==0;
X=abs(separationPBS_paired(inds));
Y=abs(separationMUSC_paired(inds));

lt_plot_45degScatter(X, Y, 'b');


% ================= 1 shared
lt_subplot(3,2,2);  hold on;
title('[same-type [1 shared syl] ');
xlabel('abs separation (PBS)')
ylabel('abs separation (MUSC)');

inds = IsSameSyl & DATSTRUCT.pairedsyls.NumSylsSharedInContext==1;
X=abs(separationPBS_paired(inds));
Y=abs(separationMUSC_paired(inds));

lt_plot_45degScatter(X, Y, 'b');

% ================= 1 shared
lt_subplot(3,2,3);  hold on;
title('[same-type [2 shared syl] ');
xlabel('abs separation (PBS)')
ylabel('abs separation (MUSC)');

inds = IsSameSyl & DATSTRUCT.pairedsyls.NumSylsSharedInContext==2;
X=abs(separationPBS_paired(inds));
Y=abs(separationMUSC_paired(inds));

lt_plot_45degScatter(X, Y, 'b');


% ================== ALL IN ONE PLOT, JUST CHANGE IN SEPARATION
lt_subplot(3,2,4); hold on;
title('change in separation');
ylabel('change in sep (MUSC-PBS)');
Yall = {};
Yallmat = [];
Xallmat = [];
% - 0 shared
inds = IsSameSyl & DATSTRUCT.pairedsyls.NumSylsSharedInContext==0;
X=abs(separationPBS_paired(inds));
Y=abs(separationMUSC_paired(inds));
Yall{1} = Y-X;
Yallmat = [Yallmat Y-X];
Xallmat = [Xallmat 0*ones(size(Y-X))];

% - 1 shared
inds = IsSameSyl & DATSTRUCT.pairedsyls.NumSylsSharedInContext==1;
X=abs(separationPBS_paired(inds));
Y=abs(separationMUSC_paired(inds));
Yall{2} = Y-X;
Yallmat = [Yallmat Y-X];
Xallmat = [Xallmat 1*ones(size(Y-X))];

% - 2 shared
inds = IsSameSyl & DATSTRUCT.pairedsyls.NumSylsSharedInContext==2;
X=abs(separationPBS_paired(inds));
Y=abs(separationMUSC_paired(inds));
Yall{3} = Y-X;
Yallmat = [Yallmat Y-X];
Xallmat = [Xallmat 2*ones(size(Y-X))];

lt_plot_MultDist(Yall, [0 1 2]);

for j=1:length(Yall);
   p = signrank(Yall{j});
   lt_plot_text(j-1, max(Yall{j}), ['p=' num2str(p)], 'r');
end


% ======================
lt_subplot(3,2,5); hold on;
xlabel('num syls shared');
ylabel('change in sep');
lt_regress(Yallmat, Xallmat, 1, 0, 1, 1, 'k');
xlim([-1 3]);

tblinput = table(Xallmat', Yallmat', 'VariableNames', {'x', 'y'});
fitlme(tblinput, 'y ~ x')


%% ==== directly compare sametype vs. diff type [all datapoints]
Yallall = {}; % cell1 = same; cell2 = diff;

% ============ SAME
inds = IsSameSyl==1;
X=abs(separationPBS_paired(inds));
Y=abs(separationMUSC_paired(inds));

% Xallall{1} = X; % separation
Yallall{1} = Y-X; % change in separation
threshold = max(X)+20000;

% ================ DIFF
inds=IsSameSyl==0;

X=abs(separationPBS_paired(inds));
Y=abs(separationMUSC_paired(inds));

% Xallall{2} = X(X<=threshold);
Yallall{2} = Y(X<=threshold) - X(X<=threshold);

% ========================================== PLOIT
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

lt_plot_MultDist(Yallall, [1 2]);

%% ==== directly compare sametype vs. diff type, [constraining to those that start similar]
Yallall = {}; % cell1 = same; cell2 = diff;

% ============ SAME
inds = IsSameSyl==1;
X=abs(separationPBS_paired(inds));
Y=abs(separationMUSC_paired(inds));

% Xallall{1} = X; % separation
Yallall{1} = Y-X; % change in separation
threshold = max(X);

% ================ DIFF
inds=IsSameSyl==0;

X=abs(separationPBS_paired(inds));
Y=abs(separationMUSC_paired(inds));

% Xallall{2} = X(X<=threshold);
Yallall{2} = Y(X<=threshold) - X(X<=threshold);

% ========================================== PLOIT
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

lt_plot_MultDist(Yallall, [1 2]);


%% %%%%%%%%%%%%%%%%%%%%%%%% change in separation as a function of Motif position
figcount=1;
subplotrows=3;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];


% ====================== SAME MOTIF
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same motif');
Yall ={};

% ----- same pair
inds = IsSameMotif==1 & IsSameSyl==1;
y = abs(separationMUSC_paired(inds)) - abs(separationPBS_paired(inds));
Yall{1} = y;

% ----- diff pair
inds = IsSameMotif==1 & IsSameSyl==0;
y = abs(separationMUSC_paired(inds)) - abs(separationPBS_paired(inds));
Yall{2} = y;

lt_plot_MultDist(Yall, [1 2])



% ====================== DIFF MOTIF
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('diff motif');
Yall ={};

% ----- same motif, same pair
inds = IsSameMotif==0 & IsSameSyl==1;
y = abs(separationMUSC_paired(inds)) - abs(separationPBS_paired(inds));
Yall{1} = y;

% ----- same motif, diff pair
inds = IsSameMotif==0 & IsSameSyl==0;
y = abs(separationMUSC_paired(inds)) - abs(separationPBS_paired(inds));
Yall{2} = y;

lt_plot_MultDist(Yall, [1 2])


%% %%%%%%%%%%%%%%%%%%%%%%%% change in separation as a function of Motif position
figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];


% ====================== SAME SYL
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same SYL');
xlabel('same motif -- diff motif');
Yall ={};

% ----- same motif
inds = IsSameMotif==1 & IsSameSyl==1;
y = abs(separationMUSC_paired(inds)) - abs(separationPBS_paired(inds));
Yall{1} = y;

% ----- diff motif
inds = IsSameMotif==0 & IsSameSyl==1;
y = abs(separationMUSC_paired(inds)) - abs(separationPBS_paired(inds));
Yall{2} = y;

lt_plot_MultDist(Yall, [1 2])
ylabel('change in sep (musc - pbs)');

ylim([-120 120]);
lt_plot_zeroline;

% --- test each
for i=1:2
    p = signrank(Yall{i});
    if p<0.1
        lt_plot_text(i, max(Yall{i}), ['p=' num2str(p)], 'r');
    end
end
% --- compare
p = ranksum(Yall{1}, Yall{2});
lt_plot_text(1.5, 0.8, ['(vs)p=' num2str(p)], 'b')



% ====================== DIFF MOTIF
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('diff syl');
Yall ={};

% -----
inds = IsSameMotif==1 & IsSameSyl==0;
y = abs(separationMUSC_paired(inds)) - abs(separationPBS_paired(inds));
Yall{1} = y;

% -----
inds = IsSameMotif==0 & IsSameSyl==0;
y = abs(separationMUSC_paired(inds)) - abs(separationPBS_paired(inds));
Yall{2} = y;

lt_plot_MultDist(Yall, [1 2]);

ylim([-120 120]);
lt_plot_zeroline;

% --- test each
for i=1:2
    p = signrank(Yall{i});
    if p<0.1
        lt_plot_text(i, max(Yall{i}), ['p=' num2str(p)], 'r');
    end
end
% --- compare
p = ranksum(Yall{1}, Yall{2});
lt_plot_text(1.5, 0.8, ['(vs)p=' num2str(p)], 'b')


%% ####################### AFP BIAS MORE SIMILAR IF ON SAME MOTIF?

% ======================== SAME SYL, SAME MOTIF
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('SAME SYL, SAME MOTIF (bk = adjacent)');
inds = IsSameMotif==1 & IsSameSyl==1;
plotcol = 'b';

y = cell2mat(DATSTRUCT.pairedsyls.PitchShiftBothSyls_paired(inds)');

% -- plot syls in random order
indtmp = rand(size(y,1),1)>0.5;
y(indtmp,:) = fliplr(y(indtmp,:));

% % -- flip so that syl starting wtih higher pitch is on first column
% ffraw = DATSTRUCT.singlesyls.Pitch_PBSraw(DATSTRUCT.pairedsyls.Pairs_OriginalInds(inds,:));
% indtmp = ffraw(:,2)>ffraw(:,1);
% y(indtmp,:) = fliplr(y(indtmp,:));

lt_regress(y(:,2), y(:,1), 1, 0, 1, 1, plotcol);
% plot(y(:,1), y(:,2), 'o', 'Color', plotcol);

% -- note down those that are adjacent
nsyls = DATSTRUCT.pairedsyls.NumSylsInBetween(inds);
plot(y(nsyls==0,1), y(nsyls==0,2), 'ok');

AdjacentBiasSAME = [y(nsyls==0,1), y(nsyls==0,2)];
NonAdjacentBiasSAME = [y(nsyls>0,1), y(nsyls>0,2)];

xlim([-60 60]);
ylim([-60 60]);
lt_plot_zeroline;
lt_plot_zeroline_vert;

% ======================== SAME SYL, DIFF MOTIF
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('SAME SYL, DIFF MOTIF');
inds = IsSameMotif==0 & IsSameSyl==1;
plotcol = 'b';

y = cell2mat(DATSTRUCT.pairedsyls.PitchShiftBothSyls_paired(inds)');
% -- plot syls in random order
indtmp = rand(size(y,1),1)>0.5;
y(indtmp,:) = fliplr(y(indtmp,:));

lt_regress(y(:,2), y(:,1), 1, 0, 1, 1, plotcol);

NonAdjacentBiasSAME = [NonAdjacentBiasSAME; [y(:,1), y(:,2)]];

lt_plot_zeroline;
lt_plot_zeroline_vert;
xlim([-60 60]);
ylim([-60 60]);


% ################## quick binomial test for whether ajacent are more
% likely to have similar direction bias than chance
% [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% title('adjacent syls');
% lt_regress(AdjacentBiasSAME(:,2), AdjacentBiasSAME(:,1), 1);
% 
% lt_plot_zeroline;
% lt_plot_zeroline_vert;
% xlim([-60 60]);
% ylim([-60 60]);
% 
% [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% title('non-adjacent syls');
% lt_regress(NonAdjacentBiasSAME(:,2), NonAdjacentBiasSAME(:,1), 1);
% 
% lt_plot_zeroline;
% lt_plot_zeroline_vert;
% xlim([-60 60]);
% ylim([-60 60]);


% -------------- directly compare slopes
X = [AdjacentBiasSAME(:,1); NonAdjacentBiasSAME(:,1)];
Y = [AdjacentBiasSAME(:,2); NonAdjacentBiasSAME(:,2)];
Grp = [0*ones(size(AdjacentBiasSAME,1),1); 1*ones(size(NonAdjacentBiasSAME,1),1)];
aoctool(X, Y, Grp)


% ======================== DIFF SYL, SAME MOTIF
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('DIFF SYL, SAME MOTIF');
inds = IsSameMotif==1 & IsSameSyl==0;
plotcol = 'r';

y = cell2mat(DATSTRUCT.pairedsyls.PitchShiftBothSyls_paired(inds)');
% -- plot syls in random order
indtmp = rand(size(y,1),1)>0.5;
y(indtmp,:) = fliplr(y(indtmp,:));

lt_regress(y(:,2), y(:,1), 1, 0, 1, 1, plotcol);
% -- note down those that are adjacent
nsyls = DATSTRUCT.pairedsyls.NumSylsInBetween(inds);
plot(y(nsyls==0,1), y(nsyls==0,2), 'ok');

lt_plot_zeroline;
lt_plot_zeroline_vert;
xlim([-60 60]);
ylim([-60 60]);

AdjacentBiasSAME = [y(nsyls==0,1), y(nsyls==0,2)];
NonAdjacentBiasSAME = [y(nsyls>0,1), y(nsyls>0,2)];

% ======================== DIFF SYL, DIFF MOTIF
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('DIFF SYL, DIFF MOTIF');
inds = IsSameMotif==0 & IsSameSyl==0;
plotcol = 'r';

y = cell2mat(DATSTRUCT.pairedsyls.PitchShiftBothSyls_paired(inds)');
% -- plot syls in random order
indtmp = rand(size(y,1),1)>0.5;
y(indtmp,:) = fliplr(y(indtmp,:));

lt_regress(y(:,2), y(:,1), 1, 0, 1, 1, plotcol);
NonAdjacentBiasSAME = [NonAdjacentBiasSAME; [y(:,1), y(:,2)]];

lt_plot_zeroline;
lt_plot_zeroline_vert;
xlim([-60 60]);
ylim([-60 60]);

% ################## quick binomial test for whether ajacent are more
% likely to have similar direction bias than chance
% [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% title('adjacent syls');
% lt_regress(AdjacentBiasSAME(:,2), AdjacentBiasSAME(:,1), 1);
% 
% lt_plot_zeroline;
% lt_plot_zeroline_vert;
% xlim([-60 60]);
% ylim([-60 60]);
% 
% [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% title('non-adjacent syls');
% lt_regress(NonAdjacentBiasSAME(:,2), NonAdjacentBiasSAME(:,1), 1);
% 
% lt_plot_zeroline;
% lt_plot_zeroline_vert;
% xlim([-60 60]);
% ylim([-60 60]);


% -------------- directly compare slopes
X = [AdjacentBiasSAME(:,1); NonAdjacentBiasSAME(:,1)];
Y = [AdjacentBiasSAME(:,2); NonAdjacentBiasSAME(:,2)];
Grp = [0*ones(size(AdjacentBiasSAME,1),1); 1*ones(size(NonAdjacentBiasSAME,1),1)];
aoctool(X, Y, Grp)

%% ####################################### SCATTER PLOT [SEPARATION]
ChangeInSep = abs(separationMUSC_paired) - abs(separationPBS_paired);
% ======== SAME TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same (bk = samemotif; rd = diff)');
inds = IsSameSyl==1;

x = abs(separationPBS_paired(inds));
y = ChangeInSep(inds);
samemotif = IsSameMotif(inds);

% -- same motifs
plot(x(samemotif==1), y(samemotif==1), 'ok');

% -- diff motif
plot(x(samemotif==0), y(samemotif==0), 'or');

lt_plot_zeroline;
ylim([-100 100]);
aoctool(x, y, samemotif)

% ======== DIFF TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('diff (bk = samemotif; rd = diff)');
inds = IsSameSyl==0;

x = abs(separationPBS_paired(inds));
y = ChangeInSep(inds);
samemotif = IsSameMotif(inds);

% -- same motifs
plot(x(samemotif==1), y(samemotif==1), 'ok');

% -- diff motif
plot(x(samemotif==0), y(samemotif==0), 'or');
lt_plot_zeroline;
ylim([-100 100]);



%% ===== same diff motif (45 degree scatter) [SEPARATION]

% ####################################### SCATTER PLOT
ChangeInSep = abs(separationMUSC_paired) - abs(separationPBS_paired);
% ======== SAME TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same (bk = samemotif; rd = diff)');
inds = IsSameSyl==1;

x = abs(separationPBS_paired(inds));
y = abs(separationMUSC_paired(inds));
samemotif = IsSameMotif(inds);

% -- same motifs
plot(x(samemotif==1), y(samemotif==1), 'ok');

% -- diff motif
plot(x(samemotif==0), y(samemotif==0), 'or');

lt_plot_zeroline;

% ylim([-100 100]);

aoctool(x, y, samemotif)



% ======== DIFF TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same (bk = samemotif; rd = diff)');
inds = IsSameSyl==0;

x = abs(separationPBS_paired(inds));
y = abs(separationMUSC_paired(inds));
samemotif = IsSameMotif(inds);

% -- same motifs
plot(x(samemotif==1), y(samemotif==1), 'ok');

% -- diff motif
plot(x(samemotif==0), y(samemotif==0), 'or');

lt_plot_zeroline;
% ylim([-100 100]);

aoctool(x, y, samemotif)


%% ######################33 [SEAPRATION], DEPENDENCE ON  POSOTIION WITHIN MOTIF
figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];


ChangeInSep = abs(separationMUSC_paired) - abs(separationPBS_paired);

inds = IsSameMotif==1 & ~isnan(NumSylsInBetween);

nsyls = NumSylsInBetween(inds);
syltype = IsSameSyl(inds);
Ycorrchange = ChangeInSep(inds);

% ========================== SAME TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same type');
ylabel('separationFF (MUSC minus PBS');
xlabel('num syls in between');
x = nsyls(syltype==1);
y = Ycorrchange(syltype==1);
plot(x, y, 'ok');
% - get mean
[ymean, ysem] = grpstats(y, x, {'mean', 'sem'});
lt_plot(unique(x)+0.1, ymean, {'Errors', ysem, 'Color', 'r'});
lt_plot_zeroline;
xlim([-1 max(x)+1]);
% ylim([-1 1]);
mean(y(x==0)), mean(y(x~=0));

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same type');
% -- OLS regression
tblinput = table(x', y', 'VariableNames', {'x', 'y'});
lt_regress(y, x, 1, 0, 1, 1, 'k')
fitlme(tblinput, 'y ~ x')
lt_plot_zeroline;
xlim([-1 max(x)+1]);
% ylim([-1 1]);

% ============================= DIFF TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('diff type');

x = nsyls(syltype==0);
y = Ycorrchange(syltype==0);
plot(x, y, 'ok');
% - get mean
[ymean, ysem] = grpstats(y, x, {'mean', 'sem'});
lt_plot(unique(x)+0.1, ymean, {'Errors', ysem, 'Color', 'r'});
lt_plot_zeroline;
xlim([-1 max(x)+1]);
% ylim([-1 1]);

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('diff type');
% -- OLS regression
tblinput = table(x', y', 'VariableNames', {'x', 'y'});
lt_regress(y, x, 1, 0, 1, 1, 'k')
fitlme(tblinput, 'y ~ x')
lt_plot_zeroline;
xlim([-1 max(x)+1]);
% ylim([-1 1]);



%% ==== effect of musc on acoustic separation,

figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];


% === ALL PAIRS
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[all pairs]');
color='k';
xlabel('Acoustic Sep (PBS)')
ylabel('Acoustic Sep (MUSC)')
inds=find(AcousticDist_pairs>0);

X=PairedAcousticDist(inds);
Y=PairedAcousticDist_MUSC(inds);

lt_plot_45degScatter(X, Y, color)


% === SAME TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[same]');
color='b';
xlabel('Acoustic Sep (PBS)')
ylabel('Acoustic Sep (MUSC)')
if useHandLabForSametype==1
    inds = IsSameSyl==1;
else
    inds=AcousticDist_pairs<same_type_thr;
end

X=PairedAcousticDist(inds);
Y=PairedAcousticDist_MUSC(inds);

lt_plot_45degScatter(X, Y, color)


% === plot all pairwise effect (45 degree)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[diff]');
color='r';
xlabel('Acoustic Sep (PBS)')
ylabel('Acoustic Sep (MUSC)')
if useHandLabForSametype==1
    inds = IsSameSyl==0;
else
    inds=AcousticDist_pairs>=same_type_thr;
end

X=PairedAcousticDist(inds);
Y=PairedAcousticDist_MUSC(inds);

lt_plot_45degScatter(X, Y, color)


%% ==== change in acoustic separation (vs. MUSC acoustic sep)
figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];


% === ALL
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[all pairs]');
color='k';
xlabel('Acoustic Sep (MUSC)')
ylabel('Change in sep (PBS - MUSC)')
inds=find(AcousticDist_pairs>0);

X=PairedAcousticDist_MUSC(inds);
Y=PairedAcousticDist(inds) - PairedAcousticDist_MUSC(inds);

lt_regress(Y, X, 1, 0, 1, 1, color);
lt_plot_zeroline;

% === SAME TYPE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[same]');
color='b';
xlabel('Acoustic Sep (MUSC)')
ylabel('Change in sep (PBS - MUSC)')
if useHandLabForSametype==1
    inds = IsSameSyl==1;
else
    inds=AcousticDist_pairs<same_type_thr;
end

X=PairedAcousticDist_MUSC(inds);
Y=PairedAcousticDist(inds) - PairedAcousticDist_MUSC(inds);

lt_regress(Y, X, 1, 0, 1, 1, color);
lt_plot_zeroline;

% === plot all pairwise effect (45 degree)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[diff]');
color='r';
xlabel('Acoustic Sep (MUSC)')
ylabel('Change in sep (PBS - MUSC)')
if useHandLabForSametype==1
    inds = IsSameSyl==0;
else
    inds=AcousticDist_pairs>=same_type_thr;
end

X=PairedAcousticDist_MUSC(inds);
Y=PairedAcousticDist(inds) - PairedAcousticDist_MUSC(inds);

lt_regress(Y, X, 1, 0, 1, 1, color);
lt_plot_zeroline;

end




function structout = fn_getAllShiftDistributions(inds, PitchShiftBothSyls_paired, ...
    separationPBS_paired, separationMUSC_paired, Pairs_OriginalInds, Pitch_PBSraw, ...
    plotON)

NAll = [];
NlessSepAll = [];
PitchShiftsAll = {};
PitchPBSAll = {};
SeparationPBS = {};
SeparationMUSC = {};

% -------------------- EXTRACT THINGS
XY = cell2mat(PitchShiftBothSyls_paired(inds)');
YsepPBS = abs(separationPBS_paired(inds));
if (0)
YsepMUSC = abs(separationMUSC_paired(inds));
else
YsepMUSC = Pitch_PBSraw(Pairs_OriginalInds(inds, :)) + XY;
YsepMUSC = abs(diff(YsepMUSC,1,2))';
end

% -- sort so that syl1 is the one with higher hz
ffsyl1 = Pitch_PBSraw(Pairs_OriginalInds(inds, 1));
ffsyl2 = Pitch_PBSraw(Pairs_OriginalInds(inds, 2));

% -- for any cases where ffsyl2 > ffsyl1, flip them
indflip = ffsyl2>ffsyl1;
XY(indflip,:) = fliplr(XY(indflip,:));
ffsyl2COPY = ffsyl2;
ffsyl2(indflip) = ffsyl1(indflip);
ffsyl1(indflip) = ffsyl2COPY(indflip);
clear ffsyl2COPY


% -- annotate cases where MUSC brought separation closer
indsMuscCloser = YsepMUSC < YsepPBS;

% ================ GO THRU CATEGORIES ONE BY ONE AND EXTRACT 1) N, 2) HOW
% MANY DECREASED SEPARATION, 3) PITCH SHIFTS
% ########## 1) down, up [i.e. syl 1 down, syl 2 up]
indtmp = XY(:,1) < 0 & XY(:,2) > 0 ;
NAll = [NAll; sum(indtmp)];
NlessSepAll = [NlessSepAll; sum(indsMuscCloser(indtmp))];
if ~any(indtmp)
PitchShiftsAll = [PitchShiftsAll; {[]}];
PitchPBSAll = [PitchPBSAll; {[]}];
SeparationPBS = [SeparationPBS; {[]}];
SeparationMUSC = [SeparationMUSC; {[]}];
else
PitchShiftsAll = [PitchShiftsAll; XY(indtmp,:)];
PitchPBSAll = [PitchPBSAll; [ffsyl1(indtmp)' ffsyl2(indtmp)']];
SeparationPBS = [SeparationPBS; YsepPBS(indtmp)];
SeparationMUSC = [SeparationMUSC; YsepMUSC(indtmp)];
end

% ########## 2) up(less), up(more)
indtmp = XY(:,1) > 0 & XY(:,2) > 0 & XY(:,1)<XY(:,2);
NAll = [NAll; sum(indtmp)];
NlessSepAll = [NlessSepAll; sum(indsMuscCloser(indtmp))];
if ~any(indtmp)
PitchShiftsAll = [PitchShiftsAll; {[]}];
PitchPBSAll = [PitchPBSAll; {[]}];
SeparationPBS = [SeparationPBS; {[]}];
SeparationMUSC = [SeparationMUSC; {[]}];
else
PitchShiftsAll = [PitchShiftsAll; XY(indtmp,:)];
PitchPBSAll = [PitchPBSAll; [ffsyl1(indtmp)' ffsyl2(indtmp)']];
SeparationPBS = [SeparationPBS; YsepPBS(indtmp)];
SeparationMUSC = [SeparationMUSC; YsepMUSC(indtmp)];
end

% ########## 3) dn(more), dn(less)
indtmp = XY(:,1) < 0 & XY(:,2) < 0 & XY(:,1)<XY(:,2);
NAll = [NAll; sum(indtmp)];
NlessSepAll = [NlessSepAll; sum(indsMuscCloser(indtmp))];
if ~any(indtmp)
PitchShiftsAll = [PitchShiftsAll; {[]}];
PitchPBSAll = [PitchPBSAll; {[]}];
SeparationPBS = [SeparationPBS; {[]}];
SeparationMUSC = [SeparationMUSC; {[]}];
else
PitchShiftsAll = [PitchShiftsAll; XY(indtmp,:)];
PitchPBSAll = [PitchPBSAll; [ffsyl1(indtmp)' ffsyl2(indtmp)']];
SeparationPBS = [SeparationPBS; YsepPBS(indtmp)];
SeparationMUSC = [SeparationMUSC; YsepMUSC(indtmp)];
end

% ########## 4) up, dn
indtmp = XY(:,1) > 0 & XY(:,2) < 0 ;
NAll = [NAll; sum(indtmp)];
NlessSepAll = [NlessSepAll; sum(indsMuscCloser(indtmp))];
if ~any(indtmp)
PitchShiftsAll = [PitchShiftsAll; {[]}];
PitchPBSAll = [PitchPBSAll; {[]}];
SeparationPBS = [SeparationPBS; {[]}];
SeparationMUSC = [SeparationMUSC; {[]}];
else
PitchShiftsAll = [PitchShiftsAll; XY(indtmp,:)];
PitchPBSAll = [PitchPBSAll; [ffsyl1(indtmp)' ffsyl2(indtmp)']];
SeparationPBS = [SeparationPBS; YsepPBS(indtmp)];
SeparationMUSC = [SeparationMUSC; YsepMUSC(indtmp)];
end

% ########## 2) up(more), up(less)
indtmp = XY(:,1) > 0 & XY(:,2) > 0 & XY(:,1)>XY(:,2);
NAll = [NAll; sum(indtmp)];
NlessSepAll = [NlessSepAll; sum(indsMuscCloser(indtmp))];
if ~any(indtmp)
PitchShiftsAll = [PitchShiftsAll; {[]}];
PitchPBSAll = [PitchPBSAll; {[]}];
SeparationPBS = [SeparationPBS; {[]}];
SeparationMUSC = [SeparationMUSC; {[]}];
else
PitchShiftsAll = [PitchShiftsAll; XY(indtmp,:)];
PitchPBSAll = [PitchPBSAll; [ffsyl1(indtmp)' ffsyl2(indtmp)']];
SeparationPBS = [SeparationPBS; YsepPBS(indtmp)];
SeparationMUSC = [SeparationMUSC; YsepMUSC(indtmp)];
end

% ########## 3) dn(less), dn(more)
indtmp = XY(:,1) < 0 & XY(:,2) < 0 & XY(:,1)>XY(:,2);
NAll = [NAll; sum(indtmp)];
NlessSepAll = [NlessSepAll; sum(indsMuscCloser(indtmp))];
if ~any(indtmp)
PitchShiftsAll = [PitchShiftsAll; {[]}];
PitchPBSAll = [PitchPBSAll; {[]}];
SeparationPBS = [SeparationPBS; {[]}];
SeparationMUSC = [SeparationMUSC; {[]}];
else
PitchShiftsAll = [PitchShiftsAll; XY(indtmp,:)];
PitchPBSAll = [PitchPBSAll; [ffsyl1(indtmp)' ffsyl2(indtmp)']];
SeparationPBS = [SeparationPBS; YsepPBS(indtmp)];
SeparationMUSC = [SeparationMUSC; YsepMUSC(indtmp)];
end





% ################################### PLOT
if plotON==1
lt_figure; hold on;

% =========== 1) bar plot of N(with fill showing how many reduced
% separation)
lt_subplot(4,1,1); hold on;
ylabel('Num cases (red = sep decrease)');
X = 1:length(NAll);
lt_plot_bar(X, NAll, {'Color', [0.6 0.6 0.6]});
lt_plot_bar(X, NlessSepAll, {'Color', 'r'});

% =========== 2) pitch shifts
lt_subplot(4,1,2); hold on
ylabel('pitch shift (hz');
for j=1:length(PitchShiftsAll)
   
    plot([j-0.2 j+0.2], PitchShiftsAll{j}','-ok');

end
lt_plot_zeroline

% =========== 3) pitch (PBS) -- pitch (MUSC)
lt_subplot(4,1,3); hold on
% ylabel('pitch shift (hz');
for j=1:length(PitchShiftsAll)
    
    numcases = size(PitchShiftsAll{j},1);
    plotcols = lt_make_plot_colors(numcases, 0, 0);
    for jj=1:numcases
    
    ffPBS = PitchPBSAll{j}(jj,:)';
    ffPBS = ffPBS-mean(ffPBS); % subtract mean
    ffMUSC = ffPBS + PitchShiftsAll{j}(jj,:)';
    
    % plot 2 contexts separate x places
%      plot([j-0.2 j+0.2], [ffPBS ffMUSC]', '-', 'Color', plotcols{jj});
    plot([j-0.3 j-0.2], [ffPBS(1) ffMUSC(1)]', '-', 'Color', plotcols{jj}, 'LineWidth', 2);
    plot([j+0.2 j+0.3], [ffPBS(2) ffMUSC(2)]', '-', 'Color', plotcols{jj}, 'LineWidth', 2);
    plot([j-0.2 j+0.2], [ffMUSC(1) ffPBS(2)], ':', 'Color', [0.7 0.7 0.7]);
        
    end
    
%     
%     plot([j-0.2 j+0.2], PitchShiftsAll{j}','-ok');

end
lt_plot_zeroline


% =========== 4) separtiaon (PBS) -- sep (MUSC)
lt_subplot(4,1,4); hold on
title('seprataiopn');
for j=1:length(PitchShiftsAll)
   
    yPBS = SeparationPBS{j};
    yMUSC = SeparationMUSC{j};
    plot([j-0.2 j+0.2], [yPBS' yMUSC']','-ok');

end
lt_plot_zeroline
end


% ================== OUTPUT
structout.NAll = NAll;
structout.NlessSepAll = NlessSepAll;
structout.PitchShiftsAll = PitchShiftsAll;
structout.PitchPBSAll = PitchPBSAll;
structout.SeparationPBS =SeparationPBS;
structout.SeparationMUSC = SeparationMUSC;


end