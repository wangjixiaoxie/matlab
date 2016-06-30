%% LT 4/26/16 - plot leraning distributions, acounting for hand vs. computer classified

function lt_seq_dep_pitch_ACROSSBIRDS_LearDistribv2(SeqDepPitch_AcrossBirds, Params, acoustThresh)



figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];



%% 
NumBirds=length(SeqDepPitch_AcrossBirds.birds);


%% not including targets
%             Targetstatus=[];
            Similar_com=[];
            Similar_hand=[];
            Generalization=[];
            AcousticDist=[];
            Shift_targ_dir=[];
            Targshift_targdir=[];
            Corr_Motif=[];
            
            Shift_targ_dir_WNend=[];
            Generalization_WNend=[];
            
            BirdnameAll={};
            SylAll={};
            ExptnameAll={};
            
            
for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    for ii=1:numexpts
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        for j=1:length(SylsUnique)
            syl=SylsUnique{j};
            
            istarg=SeqDepPitch_AcrossBirds.birds{i}.experiment...
                {ii}.Syl_ID_Dimensions.(syl).is_target;

            if istarg==1
                continue
            end
            
            BirdnameAll=[BirdnameAll birdname];
            SylAll=[SylAll syl];
            ExptnameAll=[ExptnameAll exptname];
            
            sim_com=SeqDepPitch_AcrossBirds.birds{i}.experiment...
                {ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
            Similar_com=[Similar_com sim_com];

            sim_hand=SeqDepPitch_AcrossBirds.birds{i}.experiment...
                {ii}.Syl_ID_Dimensions.(syl).similar_to_targ_HandLab;
            Similar_hand=[Similar_hand sim_hand];
            
            gener=SeqDepPitch_AcrossBirds.birds{i}.experiment...
                {ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean_rel_targ;
            Generalization=[Generalization gener];
            
            acoust=SeqDepPitch_AcrossBirds.birds{i}.experiment...
                {ii}.Syl_ID_Dimensions.(syl).eucldist_from_targ_zscore;
            AcousticDist=[AcousticDist acoust];
            
            shift=SeqDepPitch_AcrossBirds.birds{i}.experiment...
                {ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean;
            targlearndir=SeqDepPitch_AcrossBirds.birds{i}.experiment...
                {ii}.INFORMATION.targ_learn_dir;
            shift=shift*targlearndir;
            Shift_targ_dir=[Shift_targ_dir shift];
            
            targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
            try
                if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment...
                        {ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.corrcoeff_vs, targsyl);
                    corr=SeqDepPitch_AcrossBirds.birds{i}.experiment...
                        {ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.corrcoeff_vs.(targsyl);
                else
                    corr=nan;
                end
            catch err
                corr = nan;
                disp([birdname '-' exptname '; NO CORR!!']);
            end
            Corr_Motif=[Corr_Motif corr];
            
            targshift=SeqDepPitch_AcrossBirds.birds{i}.experiment...
                {ii}.Syl_ID_Dimensions.(targsyl).LEARNING.learning_metric.mean;
            targshift=targshift*targlearndir;
            Targshift_targdir=[Targshift_targdir targshift];
            
            % ==== end of WN shift
            shift_end=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_ZSCORE.Epoch.EndOfSingleTargEpoch.DATA.(syl).mean_Zscore;
            shift_end=shift_end*targlearndir;
            shift_gen=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_ZSCORE.Epoch.EndOfSingleTargEpoch.DATA.(syl).Generalization_zscore;
            
            Shift_targ_dir_WNend=[Shift_targ_dir_WNend shift_end];
            Generalization_WNend=[Generalization_WNend shift_gen];
            
            
            disp([num2str(i) '-' num2str(ii) '-' syl])
        end
    
    end
end


%% 1) ===== plot computer labeled [generalization, learn days]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
ylabel('genearlization');
xlabel('S(com) - D(com) - S(hand) - D(hand) - D(com)+S(hand) - D(com)+D(hand)');


% -- same type (computer)
x=1;
inds=Similar_com==1;

Y=Generalization(inds);
plot(x+0.1, Y, 'ok');
lt_plot_bar(x, mean(Y), {'Errors', lt_sem(Y)});

p=signrank(Y);
if p<0.1
    lt_plot_text(x, max(Y)*1.1, ['p=' num2str(p)], 'r');
end

% -- diff type (computer)
x=2;
inds=Similar_com==0;

Y=Generalization(inds);
plot(x+0.1, Y, 'ok');
lt_plot_bar(x, mean(Y), {'Errors', lt_sem(Y)});
p=signrank(Y);
if p<0.1
    lt_plot_text(x, max(Y)*1.1, ['p=' num2str(p)], 'r');
end

% -- same type (hand)
x=3;
inds=Similar_hand==1;

Y=Generalization(inds);
plot(x+0.1, Y, 'ok');
lt_plot_bar(x, mean(Y), {'Errors', lt_sem(Y)});
p=signrank(Y);
if p<0.1
    lt_plot_text(x, max(Y)*1.1, ['p=' num2str(p)], 'r');
end


% -- diff type (hand)
x=4;
inds=Similar_hand==0;

Y=Generalization(inds);
plot(x+0.1, Y, 'ok');
lt_plot_bar(x, mean(Y), {'Errors', lt_sem(Y)});
p=signrank(Y);
if p<0.1
    lt_plot_text(x, max(Y)*1.1, ['p=' num2str(p)], 'r');
end


% -- 
x=7;
inds=Similar_com==0 & Similar_hand==1;

Y=Generalization(inds);
plot(x+0.1, Y, 'ok');
lt_plot_bar(x, mean(Y), {'Errors', lt_sem(Y)});
p=signrank(Y);
if p<0.1
    lt_plot_text(x, max(Y)*1.1, ['p=' num2str(p)], 'r');
end


% -- 
x=8;
inds=Similar_com==0 & Similar_hand==0;

Y=Generalization(inds);
plot(x+0.1, Y, 'ok');
lt_plot_bar(x, mean(Y), {'Errors', lt_sem(Y)});
p=signrank(Y);
if p<0.1
    lt_plot_text(x, max(Y)*1.1, ['p=' num2str(p)], 'r');
end


%% 1) ===== plot computer labeled [gener, end of WN]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
ylabel('genearlization (end of WN training)');
xlabel('S(com) - D(com) - S(hand) - D(hand) - D(com)+S(hand) - D(com)+D(hand)');


% -- same type (computer)
x=1;
inds=Similar_com==1;

Y=Generalization_WNend(inds);
plot(x+0.1, Y, 'ok');
lt_plot_bar(x, mean(Y), {'Errors', lt_sem(Y)});

p=signrank(Y);
if p<0.1
    lt_plot_text(x, max(Y)*1.1, ['p=' num2str(p)], 'r');
end

% -- diff type (computer)
x=2;
inds=Similar_com==0;

Y=Generalization_WNend(inds);
plot(x+0.1, Y, 'ok');
lt_plot_bar(x, mean(Y), {'Errors', lt_sem(Y)});
p=signrank(Y);
if p<0.1
    lt_plot_text(x, max(Y)*1.1, ['p=' num2str(p)], 'r');
end

% -- same type (hand)
x=3;
inds=Similar_hand==1;

Y=Generalization_WNend(inds);
plot(x+0.1, Y, 'ok');
lt_plot_bar(x, mean(Y), {'Errors', lt_sem(Y)});
p=signrank(Y);
if p<0.1
    lt_plot_text(x, max(Y)*1.1, ['p=' num2str(p)], 'r');
end


% -- diff type (hand)
x=4;
inds=Similar_hand==0;

Y=Generalization_WNend(inds);
plot(x+0.1, Y, 'ok');
lt_plot_bar(x, mean(Y), {'Errors', lt_sem(Y)});
p=signrank(Y);
if p<0.1
    lt_plot_text(x, max(Y)*1.1, ['p=' num2str(p)], 'r');
end


% -- 
x=7;
inds=Similar_com==0 & Similar_hand==1;

Y=Generalization_WNend(inds);
plot(x+0.1, Y, 'ok');
lt_plot_bar(x, mean(Y), {'Errors', lt_sem(Y)});
p=signrank(Y);
if p<0.1
    lt_plot_text(x, max(Y)*1.1, ['p=' num2str(p)], 'r');
end


% -- 
x=8;
inds=Similar_com==0 & Similar_hand==0;

Y=Generalization_WNend(inds);
plot(x+0.1, Y, 'ok');
lt_plot_bar(x, mean(Y), {'Errors', lt_sem(Y)});
p=signrank(Y);
if p<0.1
    lt_plot_text(x, max(Y)*1.1, ['p=' num2str(p)], 'r');
end


%% 1) ===== plot computer labeled [shift, learn days]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
ylabel('shift (z-score)');
xlabel('S(com) - D(com) - S(hand) - D(hand) - D(com)+S(hand) - D(com)+D(hand) - TARG');


% -- same type (computer)
x=1;
inds=Similar_com==1;

Y=Shift_targ_dir(inds);
plot(x+0.1, Y, '.k');
lt_plot_bar(x, mean(Y), {'Errors', lt_sem(Y)});

p=signrank(Y);
if p<0.1
    lt_plot_text(x, max(Y)*1.1, ['p=' num2str(p)], 'r');
end

% -- diff type (computer)
x=2;
inds=Similar_com==0;

Y=Shift_targ_dir(inds);
plot(x+0.1, Y, '.k');
lt_plot_bar(x, mean(Y), {'Errors', lt_sem(Y)});
p=signrank(Y);
if p<0.1
    lt_plot_text(x, max(Y)*1.1, ['p=' num2str(p)], 'r');
end

% -- same type (hand)
x=3;
inds=Similar_hand==1;

Y=Shift_targ_dir(inds);
plot(x+0.1, Y, '.k');
lt_plot_bar(x, mean(Y), {'Errors', lt_sem(Y)});
p=signrank(Y);
if p<0.1
    lt_plot_text(x, max(Y)*1.1, ['p=' num2str(p)], 'r');
end


% -- diff type (hand)
x=4;
inds=Similar_hand==0;

Y=Shift_targ_dir(inds);
plot(x+0.1, Y, '.k');
lt_plot_bar(x, mean(Y), {'Errors', lt_sem(Y)});
p=signrank(Y);
if p<0.1
    lt_plot_text(x, max(Y)*1.1, ['p=' num2str(p)], 'r');
end


% -- 
x=7;
inds=Similar_com==0 & Similar_hand==1;

Y=Shift_targ_dir(inds);
plot(x+0.1, Y, '.k');
lt_plot_bar(x, mean(Y), {'Errors', lt_sem(Y)});
p=signrank(Y);
if p<0.1
    lt_plot_text(x, max(Y)*1.1, ['p=' num2str(p)], 'r');
end


% -- 
x=8;
inds=Similar_com==0 & Similar_hand==0;

Y=Shift_targ_dir(inds);
plot(x+0.1, Y, '.k');
lt_plot_bar(x, mean(Y), {'Errors', lt_sem(Y)});
p=signrank(Y);
if p<0.1
    lt_plot_text(x, max(Y)*1.1, ['p=' num2str(p)], 'r');
end

% --- targ
x=10;
Y=unique(Targshift_targdir);

plot(x+0.1, Y, '.k');
lt_plot_bar(x, mean(Y), {'Errors', lt_sem(Y)});
p=signrank(Y);
if p<0.1
    lt_plot_text(x, max(Y)*1.1, ['p=' num2str(p)], 'r');
end




%% 1) ===== plot computer labeled [shift, WN end]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
ylabel('shift (z-score) [WN end]');
xlabel('S(com) - D(com) - S(hand) - D(hand) - D(com)+S(hand) - D(com)+D(hand)');


% -- same type (computer)
x=1;
inds=Similar_com==1;

Y=Shift_targ_dir_WNend(inds);
plot(x+0.1, Y, '.k');
lt_plot_bar(x, mean(Y), {'Errors', lt_sem(Y)});

p=signrank(Y);
if p<0.1
    lt_plot_text(x, max(Y)*1.1, ['p=' num2str(p)], 'r');
end

% -- diff type (computer)
x=2;
inds=Similar_com==0;

Y=Shift_targ_dir_WNend(inds);
plot(x+0.1, Y, '.k');
lt_plot_bar(x, mean(Y), {'Errors', lt_sem(Y)});
p=signrank(Y);
if p<0.1
    lt_plot_text(x, max(Y)*1.1, ['p=' num2str(p)], 'r');
end

% -- same type (hand)
x=3;
inds=Similar_hand==1;

Y=Shift_targ_dir_WNend(inds);
plot(x+0.1, Y, '.k');
lt_plot_bar(x, mean(Y), {'Errors', lt_sem(Y)});
p=signrank(Y);
if p<0.1
    lt_plot_text(x, max(Y)*1.1, ['p=' num2str(p)], 'r');
end


% -- diff type (hand)
x=4;
inds=Similar_hand==0;

Y=Shift_targ_dir_WNend(inds);
plot(x+0.1, Y, '.k');
lt_plot_bar(x, mean(Y), {'Errors', lt_sem(Y)});
p=signrank(Y);
if p<0.1
    lt_plot_text(x, max(Y)*1.1, ['p=' num2str(p)], 'r');
end


% -- 
x=7;
inds=Similar_com==0 & Similar_hand==1;

Y=Shift_targ_dir_WNend(inds);
plot(x+0.1, Y, '.k');
lt_plot_bar(x, mean(Y), {'Errors', lt_sem(Y)});
p=signrank(Y);
if p<0.1
    lt_plot_text(x, max(Y)*1.1, ['p=' num2str(p)], 'r');
end


% -- 
x=8;
inds=Similar_com==0 & Similar_hand==0;

Y=Shift_targ_dir_WNend(inds);
plot(x+0.1, Y, '.k');
lt_plot_bar(x, mean(Y), {'Errors', lt_sem(Y)});
p=signrank(Y);
if p<0.1
    lt_plot_text(x, max(Y)*1.1, ['p=' num2str(p)], 'r');
end


%% ==== gen vs. acoustic dist (colored by hand label)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('acoustic dist');
ylabel('gen [learn days]');
title('colored by hand labeled. vert line is comp threshold')

% -- hand same
inds=Similar_hand==1;

X=AcousticDist(inds);
Y=Generalization(inds);

plot(X, Y, '.b');


% --- hand diff
inds=Similar_hand==0;

X=AcousticDist(inds);
Y=Generalization(inds);

plot(X, Y, '.r');


lt_plot_zeroline;
line([acoustThresh acoustThresh], ylim);




%% [generalization, WN end]

figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
plottext=1; % birdnames, etc

% 1) ==== gen vs. acoustic dist (colored by hand label) [WN end]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('acoustic dist');
ylabel('gen [end WN]');
title('colored by hand labeled. vert line is comp threshold')

% -- hand same
inds=Similar_hand==1;

X=AcousticDist(inds);
Y=Generalization_WNend(inds);

plot(X, Y, '.b');
if plottext==1
    for ind=find(inds);
        
    bname=BirdnameAll{ind};
    ename=ExptnameAll{ind};
    sname=SylAll{ind};
    x=AcousticDist(ind);
    y=Generalization_WNend(ind);
   stringgg=[bname(1:4) '-' ename(end-2:end) '-' sname];
    lt_plot_text(x, y, stringgg, 'k', 8);

    end
end

% --- hand diff
inds=Similar_hand==0;

X=AcousticDist(inds);
Y=Generalization_WNend(inds);

plot(X, Y, '.r');

if plottext==1
    for ind=find(inds);
        
    bname=BirdnameAll{ind};
    ename=ExptnameAll{ind};
    sname=SylAll{ind};
    x=AcousticDist(ind);
    y=Generalization_WNend(ind);
   stringgg=[bname(1:4) '-' ename(end-2:end) '-' sname];
    lt_plot_text(x, y, stringgg, 'k', 8);

    end
end


lt_plot_zeroline;
line([acoustThresh acoustThresh], ylim);

% ========== 2) ==== corr between learning and acoustic dist? [separate all types
% -- same (com)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('acoustic dist');
ylabel('gen [end WN]');
title('same (com)')

inds=Similar_com==1;

X=AcousticDist(inds);
Y=Generalization_WNend(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'b');
% -- same (hand), diff(com)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('acoustic dist');
ylabel('gen [end WN]');
title('same (hand), diff(com)')

inds=Similar_com==0 & Similar_hand==1;

X=AcousticDist(inds);
Y=Generalization_WNend(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'b');


% --  diff(hand)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('acoustic dist');
ylabel('gen [end WN]');
title('diff (hand)')

inds=Similar_hand==0;

X=AcousticDist(inds);
Y=Generalization_WNend(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'r');


%% [raw shift, targ dir, WN end]

figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

% 1) ==== gen vs. acoustic dist (colored by hand label) [WN end]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('acoustic dist');
ylabel('shift (targ dir) [end WN]');
title('colored by hand labeled. vert line is comp threshold')

% -- hand same
inds=Similar_hand==1;

X=AcousticDist(inds);
Y=Shift_targ_dir_WNend(inds);

plot(X, Y, '.b');
if plottext==1
    for ind=find(inds);
        
    bname=BirdnameAll{ind};
    ename=ExptnameAll{ind};
    sname=SylAll{ind};
    x=AcousticDist(ind);
    y=Shift_targ_dir_WNend(ind);
   stringgg=[bname(1:4) '-' ename(end-2:end) '-' sname];
    lt_plot_text(x, y, stringgg, 'k', 8);

    end
end

% --- hand diff
inds=Similar_hand==0;

X=AcousticDist(inds);
Y=Shift_targ_dir_WNend(inds);

plot(X, Y, '.r');

if plottext==1
    for ind=find(inds);
        
    bname=BirdnameAll{ind};
    ename=ExptnameAll{ind};
    sname=SylAll{ind};
    x=AcousticDist(ind);
    y=Shift_targ_dir_WNend(ind);
   stringgg=[bname(1:4) '-' ename(end-2:end) '-' sname];
    lt_plot_text(x, y, stringgg, 'k', 8);

    end
end


lt_plot_zeroline;
line([acoustThresh acoustThresh], ylim);

% ========== 2) ==== corr between learning and acoustic dist? [separate all types
% -- same (com)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('acoustic dist');
ylabel('shift (targ dir) [end WN]');
title('same (com)')

inds=Similar_com==1;

X=AcousticDist(inds);
Y=Shift_targ_dir_WNend(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'b');
% -- same (hand), diff(com)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('acoustic dist');
ylabel('shift (targ dir) [end WN]');
title('same (hand), diff(com)')

inds=Similar_com==0 & Similar_hand==1;

X=AcousticDist(inds);
Y=Shift_targ_dir_WNend(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'b');


% --  diff(hand)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('acoustic dist');
ylabel('shift (targ dir) [end WN]');
title('diff (hand)')

inds=Similar_hand==0;

X=AcousticDist(inds);
Y=Shift_targ_dir_WNend(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'r');


%% DISTRIBUTIONS OF WN END
% [generalization, WN end]

figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
plottext=1; % birdnames, etc

% 1) ==== gen vs. acoustic dist (colored by hand label) [WN end]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('generalization (WN end)');
title('all');

Y=Generalization_WNend;
[~, Xcenters]=lt_plot_histogram(Y, '', 1, 0, 1, 1, 'k');

Xcenters=linspace(Xcenters(1), Xcenters(end), round(length(Xcenters)/3));

% === same (com)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('generalization (WN end)');
title('same (com)')

inds=Similar_com==1;

Y=Generalization_WNend(inds);
[~, ~]=lt_plot_histogram(Y, Xcenters, 1, 0, 1, 1, 'b');



% -- same (hand), diff(com)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('generalization (WN end)');
title('same (hand), diff(com)')

inds=Similar_com==0 & Similar_hand==1;

Y=Generalization_WNend(inds);
[~, ~]=lt_plot_histogram(Y, Xcenters, 1, 0, 1, 1, 'b');



% --  diff(hand)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('generalization (WN end)');
title('diff(hand)')

inds=Similar_hand==0;

Y=Generalization_WNend(inds);
[~, ~]=lt_plot_histogram(Y, Xcenters, 1, 0, 1, 1, 'r');




%% ==== for diff types (computer labeled) that are hand labeled same type, do they generalize more strongly than dif ftype?
disp('KEYBOARD!!!');
keyboard

inds=Similar_com==0;

X=AcousticDist(inds);
Y=Generalization(inds);
Groups=Similar_hand(inds); % 1 = sim(hand lab)

[~, ~, ~, stats]=aoctool(X, Y, Groups, '', 'acoustic', 'gen', '1=same(hand)', 'on', 'parallel lines');

multcompare(stats, 'estimate', 'pmm')

disp('marginal mean (hand sim vs hand diff) not significantly different');
disp('slopes not different');


% ========== TAKE RANGE AND COMPARE HAND SIM VS. HAND DIFF
maxAcoustic=4.0; % should contain all overlap of diff(hand) and same(hand);

inds=Similar_com==0 & AcousticDist<maxAcoustic;

X=AcousticDist(inds);
Y=Generalization(inds);
Groups=Similar_hand(inds); % 1 = sim(hand lab)

[~, ~, ~, stats]=aoctool(X, Y, Groups, '', 'acoustic', 'gen', '1=same(hand)', 'on', 'parallel lines');

multcompare(stats, 'estimate', 'pmm')

% === significant difference between generalization?
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('diff (computer) separated into hand labeled sim and diff');

plot(1, Y(Groups==1), 'ok'); % same(hand)
lt_plot_bar(1, mean(Y(Groups==1)), {'Errors', lt_sem(Y(Groups==1))});
plot(2, Y(Groups==0), 'ok');
lt_plot_bar(2, mean(Y(Groups==0)), {'Errors', lt_sem(Y(Groups==1))});


p=ranksum(Y(Groups==1), Y(Groups==0));
lt_plot_pvalue(p, 'ranksum', 1);



%% ==== for diff types (computer labeled) that are hand labeled same type, do have different correlation with pitch corr?

inds=Similar_com==0 & ~isnan(Corr_Motif);

X=Corr_Motif(inds);
Y=Generalization(inds);
Groups=Similar_hand(inds); % 1 = sim(hand lab)

[~, ~, ~, stats]=aoctool(X, Y, Groups, '', 'corr(motif)', 'gen', '1=same(hand)', 'on', 'parallel lines');

multcompare(stats, 'estimate', 'pmm')


% ==== plot regressions
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% -- same(hand)
inds=Similar_com==0 & ~isnan(Corr_Motif) & Similar_hand==1;
X=Corr_Motif(inds);
Y=Generalization(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'b');


% -- diff(hand)
inds=Similar_com==0 & ~isnan(Corr_Motif) & Similar_hand==0;
X=Corr_Motif(inds);
Y=Generalization(inds);

lt_regress(Y, X, 1, 0, 1, 1, 'r');

line(xlim, [0 0], 'Color','k');
line([0 0], ylim, 'Color','k');


%% ======= greater antigeneralization for lower learning experiments?

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
ylabel('shift (targ dir)');
xlabel('targ shift (targ dir)');
title('hand labeled');

% -- hand labeled (diff)
inds=Similar_hand==0;
color='r';

X=Targshift_targdir(inds);
Y=Shift_targ_dir(inds);

lt_regress(Y, X, 1, 0, 1, 1, color);

% -- hand labeled (sim)
inds=Similar_hand==1;
color='b';

X=Targshift_targdir(inds);
Y=Shift_targ_dir(inds);

lt_regress(Y, X, 1, 0, 1, 1, color);

lt_plot_zeroline



% ------
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
ylabel('shift (targ dir)');
xlabel('targ shift (targ dir)');
title('computer labeled');

% -- hand labeled (diff)
inds=Similar_com==0;
color='r';

X=Targshift_targdir(inds);
Y=Shift_targ_dir(inds);

lt_regress(Y, X, 1, 0, 1, 1, color);

% -- hand labeled (sim)
inds=Similar_com==1;
color='b';

X=Targshift_targdir(inds);
Y=Shift_targ_dir(inds);

lt_regress(Y, X, 1, 0, 1, 1, color);

lt_plot_zeroline





