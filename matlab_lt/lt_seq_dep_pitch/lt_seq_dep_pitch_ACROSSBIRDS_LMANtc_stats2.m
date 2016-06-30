%% LT 2/29/16 - using output structures from the multidir and samedir analyses
function lt_seq_dep_pitch_ACROSSBIRDS_LMANtc_stats2(OUTPUT_multidir, OUTPUT_samedir, DATSTRUCT_multidir, DATSTRUCT_samedir, DayBinToUse);

lt_figure; hold on;
title('note, dealing with mult pairs for bird by taking 1st pair, should take average instead');
lt_plot_annotation(1, 'Note: need to redo these plots, but including the 2 bidir expts for wh25', 'k');

%%
figcount=1;
subplotrows=4;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];

% DayBinToUse=3;


%% ==== UNPAIRED


X1=OUTPUT_samedir.separation(DayBinToUse).MUSC./OUTPUT_samedir.separation(DayBinToUse).PBS;

X2=OUTPUT_multidir.separation(DayBinToUse).MUSC./OUTPUT_multidir.separation(DayBinToUse).PBS;



% ==== PLOT
% ======================== all [same type], (i.e excluding pu11-2 which is diff type
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
ylabel('consolidation');
title('unpaired, only same type(excluding pu11-2)');
[~, p]=ttest2(X1, X2);

plot(1.1, X1, 'ok');
plot(2.1, X2, 'ok');

lt_plot_bar([1 2], [mean(X1) mean(X2)], {'Errors', [lt_sem(X1) lt_sem(X2)]});
lt_plot_text(1.5, 1.1*max([X1 X2]), num2str(p, '%3.2g'), 'r')
xlim([0 3]);
ylim([0 1.5])
line(xlim, [1 1], 'LineStyle', '--', 'Color','k')


%% ==== PAIRED
X1=[];
X2=[];
Paired_Samedir_PBS=[];
Paired_Samedir_MUSC=[];
Paired_Multdir_PBS=[];
Paired_Multdir_MUSC=[];

for i=1:length(OUTPUT_samedir.INFORMATION(DayBinToUse).experimentNum);
    
    exptnumSD=OUTPUT_samedir.INFORMATION(DayBinToUse).experimentNum(i);
    birdnameSD=DATSTRUCT_samedir.INFORMATION(exptnumSD).birdname;
    %     exptnameSD=DATSTRUCT_samedir.INFORMATION(exptnumSD).exptname;
    
    
    % figure out if it has a multidir counterpart
    for ii=1:length(OUTPUT_multidir.INFORMATION(DayBinToUse).experimentNum);
        
        
        exptnumMD=OUTPUT_multidir.INFORMATION(DayBinToUse).experimentNum(ii);
        birdnameMD=DATSTRUCT_multidir.INFORMATION(exptnumMD).birdname;
        %         exptnameMD=DATSTRUCT_multidir.INFORMATION(exptnumMD).exptname;
        
        % ----- check if birdname and exptname are same
        if strcmp(birdnameSD, birdnameMD)
            
            % ===== COLLECT DATA FOR BOTH
            Paired_Samedir_PBS=[Paired_Samedir_PBS OUTPUT_samedir.separation(DayBinToUse).PBS(i)];
            Paired_Samedir_MUSC=[Paired_Samedir_MUSC OUTPUT_samedir.separation(DayBinToUse).MUSC(i)];
            Paired_Multdir_PBS=[Paired_Multdir_PBS OUTPUT_multidir.separation(DayBinToUse).PBS(ii)];
            Paired_Multdir_MUSC=[Paired_Multdir_MUSC OUTPUT_multidir.separation(DayBinToUse).MUSC(ii)];
            
            
            
            % x1=OUTPUT_samedir.separation(DayBinToUse).MUSC(i)./OUTPUT_samedir.separation(DayBinToUse).PBS(i);
            % x2=OUTPUT_multidir.separation(DayBinToUse).MUSC(ii)./OUTPUT_multidir.separation(DayBinToUse).PBS(ii);
            %
            % X1=[X1 x1];
            % X2=[X2 x2];
            
            disp(birdnameSD)
            break
        end
    end
end

X1=Paired_Samedir_MUSC./Paired_Samedir_PBS;
X2=Paired_Multdir_MUSC./Paired_Multdir_PBS;



            
            
% ============================= all, rd23 one syl same, one syl diff, between phases
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
ylabel('consolidation');
xlabel('samedir -- diffdir');
title('paired (one rd23 syl not overlapped)');
[~, p]=ttest(X1, X2);

plot([1.1 2.1], [X1; X2], '-ok');

lt_plot_bar([1 2], [mean(X1) mean(X2)], {'Errors', [lt_sem(X1) lt_sem(X2)]});
lt_plot_text(1.5, 1.1*max([X1 X2]), ['ttest ' num2str(p, '%3.2g')], 'r')
xlim([0 3]);
ylim([0 1.5])
line(xlim, [1 1], 'LineStyle', '--', 'Color','k')
            
        


%% ====== plot consolidation across days for both samedir and diff dir on one plot
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('day bin');
ylabel('consolidation');
title('[left: samedir; right: multidir]');

NumDayBins=3;
for i=1:NumDayBins
X=[3*i-1 3*i];
% ==== samedir
consolVals=OUTPUT_samedir.separation(i).MUSC./OUTPUT_samedir.separation(i).PBS;

Y1=mean(consolVals);
Y1err=lt_sem(consolVals);
plot(X(1)+0.1, consolVals, 'ob');

% ==== diffdir
consolVals=OUTPUT_multidir.separation(i).MUSC./OUTPUT_multidir.separation(i).PBS;

Y2=mean(consolVals);
Y2err=lt_sem(consolVals);
plot(X(2)+0.1, consolVals, 'or');

%=== PLOT
lt_plot_bar(X(1), [Y1], {'Errors', [Y1err], 'Color','b'});
lt_plot_bar(X(2), [Y2], {'Errors',  [Y2err], 'Color','r'});

end



%% ====== plot consolidation across days for both samedir and diff dir on one plot [PAIRED]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('day bin');
ylabel('consolidation');
title('[left: samedir; right: multidir]');

NumDayBins=3;
for i=1:NumDayBins
X=[3*i-1 3*i];
% ==== COLLECT SAMEDIR AND DIFFDIR
Paired_Samedir_PBS=[];
Paired_Samedir_MUSC=[];
Paired_Multdir_PBS=[];
Paired_Multdir_MUSC=[];

for j=1:length(OUTPUT_samedir.INFORMATION(i).experimentNum);
    
    exptnumSD=OUTPUT_samedir.INFORMATION(i).experimentNum(j);
    birdnameSD=DATSTRUCT_samedir.INFORMATION(exptnumSD).birdname;
    
    % figure out if it has a multidir counterpart
    for jj=1:length(OUTPUT_multidir.INFORMATION(i).experimentNum);
        
        
        exptnumMD=OUTPUT_multidir.INFORMATION(i).experimentNum(jj);
        birdnameMD=DATSTRUCT_multidir.INFORMATION(exptnumMD).birdname;
        
        % ----- check if birdname and exptname are same
        if strcmp(birdnameSD, birdnameMD)
            
            % ===== COLLECT DATA FOR BOTH
Paired_Samedir_PBS=[Paired_Samedir_PBS OUTPUT_samedir.separation(i).PBS(j)];
Paired_Samedir_MUSC=[Paired_Samedir_MUSC OUTPUT_samedir.separation(i).MUSC(j)];
Paired_Multdir_PBS=[Paired_Multdir_PBS OUTPUT_multidir.separation(i).PBS(jj)];
Paired_Multdir_MUSC=[Paired_Multdir_MUSC OUTPUT_multidir.separation(i).MUSC(jj)];

            disp(birdnameSD)
            break

        end
    end
end


% === plot lines
Y=[Paired_Samedir_MUSC./Paired_Samedir_PBS; Paired_Multdir_MUSC./Paired_Multdir_PBS];
plot(X+0.1, Y', '-k');

% === samedir
consolVals=Paired_Samedir_MUSC./Paired_Samedir_PBS;

Y1=mean(consolVals);
Y1err=lt_sem(consolVals);
plot(X(1)+0.1, consolVals, 'ob');
consolSame=consolVals;

% ==== diffdir
consolVals=Paired_Multdir_MUSC./Paired_Multdir_PBS;

Y2=mean(consolVals);
Y2err=lt_sem(consolVals);
plot(X(2)+0.1, consolVals, 'or');
consolDiff=consolVals;

%=== PLOT
lt_plot_bar(X(1), [Y1], {'Errors', [Y1err], 'Color','b'});
lt_plot_bar(X(2), [Y2], {'Errors',  [Y2err], 'Color','r'});


% === stats (paired)
    % ==== STATS (on unpaired)
    p=signrank(consolSame, consolDiff);
%     [~, p]=ttest(consolSame, consolDiff);
    
    if p<0.1;
        lt_plot_text(X(1), 1.1, ['paired p=' num2str(p)],'k');
    end
    

end


%% ====== plot consolidation across days for both samedir and diff dir on one plot [ALL, but put lines for paired]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('day bin');
ylabel('consolidation');
title('[left: samedir; right: multidir]');

NumDayBins=3;
for i=1:NumDayBins
    X=[3*i-1 3*i];
    % ==== COLLECT SAMEDIR AND DIFFDIR
    Paired_Samedir_PBS=[];
    Paired_Samedir_MUSC=[];
    Paired_Multdir_PBS=[];
    Paired_Multdir_MUSC=[];
    
    for j=1:length(OUTPUT_samedir.INFORMATION(i).experimentNum);
        
        exptnumSD=OUTPUT_samedir.INFORMATION(i).experimentNum(j);
        birdnameSD=DATSTRUCT_samedir.INFORMATION(exptnumSD).birdname;
        
        % figure out if it has a multidir counterpart
        for jj=1:length(OUTPUT_multidir.INFORMATION(i).experimentNum);
            
            
            exptnumMD=OUTPUT_multidir.INFORMATION(i).experimentNum(jj);
            birdnameMD=DATSTRUCT_multidir.INFORMATION(exptnumMD).birdname;
            
            % ----- check if birdname and exptname are same
            if strcmp(birdnameSD, birdnameMD)
                
                % ===== COLLECT DATA FOR BOTH
                Paired_Samedir_PBS=[Paired_Samedir_PBS OUTPUT_samedir.separation(i).PBS(j)];
                Paired_Samedir_MUSC=[Paired_Samedir_MUSC OUTPUT_samedir.separation(i).MUSC(j)];
                Paired_Multdir_PBS=[Paired_Multdir_PBS OUTPUT_multidir.separation(i).PBS(jj)];
                Paired_Multdir_MUSC=[Paired_Multdir_MUSC OUTPUT_multidir.separation(i).MUSC(jj)];
                
                            disp(birdnameSD)
            break

            end
        end
    end
    
    
    % === plot lines [PAIRED
    Y=[Paired_Samedir_MUSC./Paired_Samedir_PBS; Paired_Multdir_MUSC./Paired_Multdir_PBS];
    plot(X+0.1, Y', '-k');
    
    % === PLOT UNPAIRED MEANS, ETC
    % ==== samedir
    consolVals=OUTPUT_samedir.separation(i).MUSC./OUTPUT_samedir.separation(i).PBS;
    
    Y1=mean(consolVals);
    Y1err=lt_sem(consolVals);
    plot(X(1)+0.1, consolVals, '.b');
    
    consolSame=consolVals;
    
    % ==== diffdir
    consolVals=OUTPUT_multidir.separation(i).MUSC./OUTPUT_multidir.separation(i).PBS;
    
    Y2=mean(consolVals);
    Y2err=lt_sem(consolVals);
    plot(X(2)+0.1, consolVals, '.r');
     consolDiff=consolVals;
   
    %=== PLOT
    lt_plot_bar(X(1), [Y1], {'Errors', [Y1err], 'Color','b'});
    lt_plot_bar(X(2), [Y2], {'Errors',  [Y2err], 'Color','r'});
    
    
    % ==== STATS (on unpaired)
    p=ranksum(consolSame, consolDiff);
%     [~, p]=ttest2(consolSame, consolDiff);
    
    if p<0.1;
        lt_plot_text(X(1), 1.1, ['unpaired p=' num2str(p)],'k');
    end
    
    
    
end

line(xlim, [1 1], 'Color', 'k', 'LineStyle', '--');





%% ====== plot AFP BIAS across days for both samedir and diff dir on one plot [ALL, but put lines for paired]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('day bin');
ylabel('AFP BIAS');
title('[left: samedir; right: multidir]');

NumDayBins=3;
for i=1:NumDayBins
    X=[3*i-1 3*i];
    % ==== COLLECT SAMEDIR AND DIFFDIR
    Paired_Samedir_PBS=[];
    Paired_Samedir_MUSC=[];
    Paired_Multdir_PBS=[];
    Paired_Multdir_MUSC=[];
    
    for j=1:length(OUTPUT_samedir.INFORMATION(i).experimentNum);
        
        exptnumSD=OUTPUT_samedir.INFORMATION(i).experimentNum(j);
        birdnameSD=DATSTRUCT_samedir.INFORMATION(exptnumSD).birdname;
        
        % figure out if it has a multidir counterpart
        for jj=1:length(OUTPUT_multidir.INFORMATION(i).experimentNum);
            
            
            exptnumMD=OUTPUT_multidir.INFORMATION(i).experimentNum(jj);
            birdnameMD=DATSTRUCT_multidir.INFORMATION(exptnumMD).birdname;
            
            % ----- check if birdname and exptname are same
            if strcmp(birdnameSD, birdnameMD)
                
                % ===== COLLECT DATA FOR BOTH
                Paired_Samedir_PBS=[Paired_Samedir_PBS OUTPUT_samedir.separation(i).PBS(j)];
                Paired_Samedir_MUSC=[Paired_Samedir_MUSC OUTPUT_samedir.separation(i).MUSC(j)];
                Paired_Multdir_PBS=[Paired_Multdir_PBS OUTPUT_multidir.separation(i).PBS(jj)];
                Paired_Multdir_MUSC=[Paired_Multdir_MUSC OUTPUT_multidir.separation(i).MUSC(jj)];
                
                            disp(birdnameSD)
            break

            end
        end
    end
    
    
    % === plot lines [PAIRED
    Y=[Paired_Samedir_PBS-Paired_Samedir_MUSC; Paired_Multdir_PBS - Paired_Multdir_MUSC];
%     Y=[Paired_Samedir_MUSC./Paired_Samedir_PBS; Paired_Multdir_MUSC./Paired_Multdir_PBS];
    plot(X+0.1, Y', '-k');
    
    % === PLOT UNPAIRED MEANS, ETC
    % ==== samedir
    afpVals=OUTPUT_samedir.separation(i).PBS-OUTPUT_samedir.separation(i).MUSC;
%     consolVals=OUTPUT_samedir.separation(i).MUSC./OUTPUT_samedir.separation(i).PBS;
    
    Y1=mean(afpVals);
    Y1err=lt_sem(afpVals);
    plot(X(1)+0.1, afpVals, '.b');
    
    consolSame=afpVals;
    
    % ==== diffdir
        afpVals=OUTPUT_multidir.separation(i).PBS-OUTPUT_multidir.separation(i).MUSC;
% consolVals=OUTPUT_multidir.separation(i).MUSC./OUTPUT_multidir.separation(i).PBS;
    
    Y2=mean(afpVals);
    Y2err=lt_sem(afpVals);
    plot(X(2)+0.1, afpVals, '.r');
     consolDiff=afpVals;
   
    %=== PLOT
    lt_plot_bar(X(1), [Y1], {'Errors', [Y1err], 'Color','b'});
    lt_plot_bar(X(2), [Y2], {'Errors',  [Y2err], 'Color','r'});
    
    
    % ==== STATS (on unpaired)
    p=ranksum(consolSame, consolDiff);
%     [~, p]=ttest2(consolSame, consolDiff);
    
    if p<0.2;
        lt_plot_text(X(1), 1.1*max(Y1), ['unpaired p=' num2str(p)],'k');
    end
    
    % === stats (paired)
    yy1=[Paired_Samedir_PBS-Paired_Samedir_MUSC];
    yy2=[Paired_Multdir_PBS-Paired_Multdir_MUSC];
    
    p=signrank(yy1, yy2);
    
    if p<0.2;
        lt_plot_text(X(1), 1.2*max(Y1), ['paired p=' num2str(p)],'m');
    end
    
end

line(xlim, [1 1], 'Color', 'k', 'LineStyle', '--');


%% ====== plot AFP BIAS [norm laerning] across days for both samedir and diff dir on one plot [ALL, but put lines for paired]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('day bin');
ylabel('AFP BIAS/LEARNING');
title('[left: samedir; right: multidir]');

NumDayBins=3;
for i=1:NumDayBins
    X=[3*i-1 3*i];
    % ==== COLLECT SAMEDIR AND DIFFDIR
    Paired_Samedir_PBS=[];
    Paired_Samedir_MUSC=[];
    Paired_Multdir_PBS=[];
    Paired_Multdir_MUSC=[];
    
    for j=1:length(OUTPUT_samedir.INFORMATION(i).experimentNum);
        
        exptnumSD=OUTPUT_samedir.INFORMATION(i).experimentNum(j);
        birdnameSD=DATSTRUCT_samedir.INFORMATION(exptnumSD).birdname;
        
        % figure out if it has a multidir counterpart
        for jj=1:length(OUTPUT_multidir.INFORMATION(i).experimentNum);
            
            
            exptnumMD=OUTPUT_multidir.INFORMATION(i).experimentNum(jj);
            birdnameMD=DATSTRUCT_multidir.INFORMATION(exptnumMD).birdname;
            
            % ----- check if birdname and exptname are same
            if strcmp(birdnameSD, birdnameMD)
                
                % ===== COLLECT DATA FOR BOTH
                Paired_Samedir_PBS=[Paired_Samedir_PBS OUTPUT_samedir.separation(i).PBS(j)];
                Paired_Samedir_MUSC=[Paired_Samedir_MUSC OUTPUT_samedir.separation(i).MUSC(j)];
                Paired_Multdir_PBS=[Paired_Multdir_PBS OUTPUT_multidir.separation(i).PBS(jj)];
                Paired_Multdir_MUSC=[Paired_Multdir_MUSC OUTPUT_multidir.separation(i).MUSC(jj)];
                
                            disp(birdnameSD)
            break

            end
        end
    end
    afpVals
    
    % === plot lines [PAIRED
    Y=[(Paired_Samedir_PBS-Paired_Samedir_MUSC)./Paired_Samedir_PBS; (Paired_Multdir_PBS - Paired_Multdir_MUSC)./Paired_Multdir_PBS];
    plot(X+0.1, Y', '-k');
    
    % === PLOT UNPAIRED MEANS, ETC
    % ==== samedir
    afpVals=OUTPUT_samedir.separation(i).PBS-OUTPUT_samedir.separation(i).MUSC;
    learningVals=OUTPUT_samedir.separation(i).PBS;
    Yvals=afpVals./learningVals;
    
    Y1=mean(Yvals);
    Y1err=lt_sem(Yvals);
    plot(X(1)+0.1, Yvals, '.b');
    
    consolSame=Yvals;
    
    % ==== diffdir
    afpVals=OUTPUT_multidir.separation(i).PBS-OUTPUT_multidir.separation(i).MUSC;
    learningVals=OUTPUT_multidir.separation(i).PBS;
    Yvals=afpVals./learningVals;
    
    Y2=mean(Yvals);
    Y2err=lt_sem(Yvals);
    plot(X(2)+0.1, Yvals, '.r');
     consolDiff=Yvals;
   
    %=== PLOT
    lt_plot_bar(X(1), [Y1], {'Errors', [Y1err], 'Color','b'});
    lt_plot_bar(X(2), [Y2], {'Errors',  [Y2err], 'Color','r'});
    
    
    % ==== STATS (on unpaired)
    p=ranksum(consolSame, consolDiff);
    
    if p<0.1;
        lt_plot_text(X(1), 1.1, ['unpaired p=' num2str(p)],'k');
    end
    

    % === stats (paired)
    yy1=(Paired_Samedir_PBS-Paired_Samedir_MUSC)./Paired_Samedir_PBS;
    yy2=(Paired_Multdir_PBS - Paired_Multdir_MUSC)./Paired_Multdir_PBS;
    
    p=signrank(yy1, yy2);
    
    if p<0.2;
        lt_plot_text(X(1), 1.2*max(Y1), ['paired p=' num2str(p)],'m');
    end

    
    
end

line(xlim, [1 1], 'Color', 'k', 'LineStyle', '--');



%% ====== plot consolidation across days for both samedir and diff dir on one plot [ALL] [ONLY THE BIDIR EXPERIMENTS IN WHICH BOTH PAST BASELINE]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('day bin');
ylabel('consolidation');
title('[left: samedir; right: multidir] [BIDIR: only expts where past baseline [could be diff expt for each bin]');

NumDayBins=3;
for i=1:NumDayBins
    X=[3*i-1 3*i];
    
    % === PLOT UNPAIRED MEANS, ETC
    % ==== samedir
    consolVals=OUTPUT_samedir.separation(i).MUSC./OUTPUT_samedir.separation(i).PBS;
    
    Y1=mean(consolVals);
    Y1err=lt_sem(consolVals);
    plot(X(1)+0.1, consolVals, '.b');
    
    consolSame=consolVals;
    
    % ==== diffdir [first filter out expts that have opposite sign for
    % each targ
    IndsToKeep=find(sign(OUTPUT_multidir.secondtarg.FFrelBaseInBin(i).PBS)==-1); % for these second targ is neg learning
    consolVals=OUTPUT_multidir.separation(i).MUSC(IndsToKeep)./OUTPUT_multidir.separation(i).PBS(IndsToKeep);
    
    Y2=mean(consolVals);
    Y2err=lt_sem(consolVals);
    plot(X(2)+0.1, consolVals, '.r');
     consolDiff=consolVals;
   
    %=== PLOT
    lt_plot_bar(X(1), [Y1], {'Errors', [Y1err], 'Color','b'});
    lt_plot_bar(X(2), [Y2], {'Errors',  [Y2err], 'Color','r'});
    
    
    % ==== STATS (on unpaired)
    p=ranksum(consolSame, consolDiff);
%     [~, p]=ttest2(consolSame, consolDiff);
    
    if p<0.1;
        lt_plot_text(X(1), 1.1, ['unpaired p=' num2str(p)],'k');
    end
    
    
    
end

line(xlim, [1 1], 'Color', 'k', 'LineStyle', '--');


%% ====== plot consolidation across days for both samedir and diff dir on one plot [ALL, but paired lines] [1ST TARG ONLY]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('day bin');
ylabel('consolidation');
title('[left: samedir; right: multidir] !!1st targ');

NumDayBins=3;
for i=1:NumDayBins
    X=[3*i-1 3*i];

        % ==== COLLECT SAMEDIR AND DIFFDIR
    Paired_Samedir_PBS=[];
    Paired_Samedir_MUSC=[];
    Paired_Multdir_PBS=[];
    Paired_Multdir_MUSC=[];
    
    for j=1:length(OUTPUT_samedir.INFORMATION(i).experimentNum);
        
        exptnumSD=OUTPUT_samedir.INFORMATION(i).experimentNum(j);
        birdnameSD=DATSTRUCT_samedir.INFORMATION(exptnumSD).birdname;
        
        % figure out if it has a multidir counterpart
        for jj=1:length(OUTPUT_multidir.INFORMATION(i).experimentNum);
            
            
            exptnumMD=OUTPUT_multidir.INFORMATION(i).experimentNum(jj);
            birdnameMD=DATSTRUCT_multidir.INFORMATION(exptnumMD).birdname;
            
            % ----- check if birdname and exptname are same
            if strcmp(birdnameSD, birdnameMD)
                
%                 % ===== COLLECT DATA FOR BOTH
%                 Paired_Samedir_PBS=[Paired_Samedir_PBS OUTPUT_samedir.separation(i).PBS(j)];
%                 Paired_Samedir_MUSC=[Paired_Samedir_MUSC OUTPUT_samedir.separation(i).MUSC(j)];
%                 Paired_Multdir_PBS=[Paired_Multdir_PBS OUTPUT_multidir.separation(i).PBS(jj)];
%                 Paired_Multdir_MUSC=[Paired_Multdir_MUSC OUTPUT_multidir.separation(i).MUSC(jj)];
                
                % ===== COLLECT DATA FOR BOTH
                Paired_Samedir_PBS=[Paired_Samedir_PBS OUTPUT_samedir.firsttarget.FFrelBaseInBin(i).PBS(j)];
                Paired_Samedir_MUSC=[Paired_Samedir_MUSC OUTPUT_samedir.firsttarget.FFrelBaseInBin(i).MUSC(j)];
                Paired_Multdir_PBS=[Paired_Multdir_PBS OUTPUT_multidir.firsttarget.FFrelBaseInBin(i).PBS(jj)];
                Paired_Multdir_MUSC=[Paired_Multdir_MUSC OUTPUT_multidir.firsttarget.FFrelBaseInBin(i).MUSC(jj)];

                            disp(birdnameSD)
            break

            end
        end
    end

    % === plot lines [PAIRED
    Y=[Paired_Samedir_MUSC./Paired_Samedir_PBS; Paired_Multdir_MUSC./Paired_Multdir_PBS];
    plot(X+0.1, Y', '-k');

    
    % === PLOT UNPAIRED MEANS, ETC
    % ==== samedir
    consolVals=OUTPUT_samedir.firsttarget.FFrelBaseInBin(i).MUSC./OUTPUT_samedir.firsttarget.FFrelBaseInBin(i).PBS;
    
    Y1=mean(consolVals);
    Y1err=lt_sem(consolVals);
    plot(X(1)+0.1, consolVals, '.b');
    
    consolSame=consolVals;
    
    % ==== diffdir
    consolVals=OUTPUT_multidir.firsttarget.FFrelBaseInBin(i).MUSC./OUTPUT_multidir.firsttarget.FFrelBaseInBin(i).PBS;
    
    Y2=mean(consolVals);
    Y2err=lt_sem(consolVals);
    plot(X(2)+0.1, consolVals, '.r');
     consolDiff=consolVals;
   
    %=== PLOT
    lt_plot_bar(X(1), [Y1], {'Errors', [Y1err], 'Color','b'});
    lt_plot_bar(X(2), [Y2], {'Errors',  [Y2err], 'Color','r'});
    
    
    % ==== STATS (on unpaired)
    p=ranksum(consolSame, consolDiff);
%     [~, p]=ttest2(consolSame, consolDiff);
    
    if p<0.1;
        lt_plot_text(X(1), 1.1, ['unpaired p=' num2str(p)],'k');
    end
    
    % === stats (paired)
    yy1=Paired_Samedir_MUSC./Paired_Samedir_PBS;
    yy2=Paired_Multdir_MUSC./Paired_Multdir_PBS;
    
    p=signrank(yy1, yy2);
    
    if p<0.2;
        lt_plot_text(X(1), 1.2*max(Y1), ['paired p=' num2str(p)],'m');
    end

end

line(xlim, [1 1], 'Color', 'k', 'LineStyle', '--');


%% ====== [1ST TARG] plot AFP BIAS across days for both samedir and diff dir on one plot [ALL, but put lines for paired]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('day bin');
ylabel('AFP BIAS');
title('[left: samedir; right: multidir] [1st targ]');

NumDayBins=3;
for i=1:NumDayBins
    X=[3*i-1 3*i];
    % ==== COLLECT SAMEDIR AND DIFFDIR
    Paired_Samedir_PBS=[];
    Paired_Samedir_MUSC=[];
    Paired_Multdir_PBS=[];
    Paired_Multdir_MUSC=[];
    
    for j=1:length(OUTPUT_samedir.INFORMATION(i).experimentNum);
        
        exptnumSD=OUTPUT_samedir.INFORMATION(i).experimentNum(j);
        birdnameSD=DATSTRUCT_samedir.INFORMATION(exptnumSD).birdname;
        
        % figure out if it has a multidir counterpart
        for jj=1:length(OUTPUT_multidir.INFORMATION(i).experimentNum);
            
            
            exptnumMD=OUTPUT_multidir.INFORMATION(i).experimentNum(jj);
            birdnameMD=DATSTRUCT_multidir.INFORMATION(exptnumMD).birdname;
            
            % ----- check if birdname and exptname are same
            if strcmp(birdnameSD, birdnameMD)
                
                % ===== COLLECT DATA FOR BOTH
                Paired_Samedir_PBS=[Paired_Samedir_PBS OUTPUT_samedir.firsttarget.FFrelBaseInBin(i).PBS(j)];
                Paired_Samedir_MUSC=[Paired_Samedir_MUSC OUTPUT_samedir.firsttarget.FFrelBaseInBin(i).MUSC(j)];
                Paired_Multdir_PBS=[Paired_Multdir_PBS OUTPUT_multidir.firsttarget.FFrelBaseInBin(i).PBS(jj)];
                Paired_Multdir_MUSC=[Paired_Multdir_MUSC OUTPUT_multidir.firsttarget.FFrelBaseInBin(i).MUSC(jj)];
                            disp(birdnameSD)
            break

            end
        end
    end
    
    
    % === plot lines [PAIRED
    Y=[Paired_Samedir_PBS-Paired_Samedir_MUSC; Paired_Multdir_PBS - Paired_Multdir_MUSC];
%     Y=[Paired_Samedir_MUSC./Paired_Samedir_PBS; Paired_Multdir_MUSC./Paired_Multdir_PBS];
    plot(X+0.1, Y', '-k');
    
    % === PLOT UNPAIRED MEANS, ETC
    % ==== samedir
    afpVals=OUTPUT_samedir.firsttarget.FFrelBaseInBin(i).PBS-OUTPUT_samedir.firsttarget.FFrelBaseInBin(i).MUSC;
%     consolVals=OUTPUT_samedir.separation(i).MUSC./OUTPUT_samedir.separation(i).PBS;
    
    Y1=mean(afpVals);
    Y1err=lt_sem(afpVals);
    plot(X(1)+0.1, afpVals, '.b');
    
    consolSame=afpVals;
    
    % ==== diffdir
        afpVals=OUTPUT_multidir.firsttarget.FFrelBaseInBin(i).PBS-OUTPUT_multidir.firsttarget.FFrelBaseInBin(i).MUSC;
% consolVals=OUTPUT_multidir.separation(i).MUSC./OUTPUT_multidir.separation(i).PBS;
    
    Y2=mean(afpVals);
    Y2err=lt_sem(afpVals);
    plot(X(2)+0.1, afpVals, '.r');
     consolDiff=afpVals;
   
    %=== PLOT
    lt_plot_bar(X(1), [Y1], {'Errors', [Y1err], 'Color','b'});
    lt_plot_bar(X(2), [Y2], {'Errors',  [Y2err], 'Color','r'});
    
    
    % ==== STATS (on unpaired)
    p=ranksum(consolSame, consolDiff);
%     [~, p]=ttest2(consolSame, consolDiff);
    
    if p<0.2;
        lt_plot_text(X(1), 1.1*max(Y1), ['unpaired p=' num2str(p)],'k');
    end
    
    % === stats (paired)
    yy1=[Paired_Samedir_PBS-Paired_Samedir_MUSC];
    yy2=[Paired_Multdir_PBS-Paired_Multdir_MUSC];
    
    p=signrank(yy1, yy2);
    
    if p<0.2;
        lt_plot_text(X(1), 1.2*max(Y1), ['paired p=' num2str(p)],'m');
    end
    
end

line(xlim, [1 1], 'Color', 'k', 'LineStyle', '--');


%% ====== [1st targ] plot AFP BIAS [norm laerning] across days for both samedir and diff dir on one plot [ALL, but put lines for paired]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('day bin');
ylabel('AFP BIAS/LEARNING');
title('[1st targ] [left: samedir; right: multidir]');

NumDayBins=3;
for i=1:NumDayBins
    X=[3*i-1 3*i];
    % ==== COLLECT SAMEDIR AND DIFFDIR
    Paired_Samedir_PBS=[];
    Paired_Samedir_MUSC=[];
    Paired_Multdir_PBS=[];
    Paired_Multdir_MUSC=[];
    
    for j=1:length(OUTPUT_samedir.INFORMATION(i).experimentNum);
        
        exptnumSD=OUTPUT_samedir.INFORMATION(i).experimentNum(j);
        birdnameSD=DATSTRUCT_samedir.INFORMATION(exptnumSD).birdname;
        
        % figure out if it has a multidir counterpart
        for jj=1:length(OUTPUT_multidir.INFORMATION(i).experimentNum);
            
            
            exptnumMD=OUTPUT_multidir.INFORMATION(i).experimentNum(jj);
            birdnameMD=DATSTRUCT_multidir.INFORMATION(exptnumMD).birdname;
            
            % ----- check if birdname and exptname are same
            if strcmp(birdnameSD, birdnameMD)
                
                % ===== COLLECT DATA FOR BOTH
                Paired_Samedir_PBS=[Paired_Samedir_PBS OUTPUT_samedir.firsttarget.FFrelBaseInBin(i).PBS(j)];
                Paired_Samedir_MUSC=[Paired_Samedir_MUSC OUTPUT_samedir.firsttarget.FFrelBaseInBin(i).MUSC(j)];
                Paired_Multdir_PBS=[Paired_Multdir_PBS OUTPUT_multidir.firsttarget.FFrelBaseInBin(i).PBS(jj)];
                Paired_Multdir_MUSC=[Paired_Multdir_MUSC OUTPUT_multidir.firsttarget.FFrelBaseInBin(i).MUSC(jj)];
            disp(birdnameSD)
            break

            end
        end
    end
    
    % === plot lines [PAIRED
    Y=[(Paired_Samedir_PBS-Paired_Samedir_MUSC)./Paired_Samedir_PBS; (Paired_Multdir_PBS - Paired_Multdir_MUSC)./Paired_Multdir_PBS];
    plot(X+0.1, Y', '-k');
    
    % === PLOT UNPAIRED MEANS, ETC
    % ==== samedir
    afpVals=OUTPUT_samedir.firsttarget.FFrelBaseInBin(i).PBS-OUTPUT_samedir.firsttarget.FFrelBaseInBin(i).MUSC;
    learningVals=OUTPUT_samedir.firsttarget.FFrelBaseInBin(i).PBS;
    Yvals=afpVals./learningVals;
    
    Y1=mean(Yvals);
    Y1err=lt_sem(Yvals);
    plot(X(1)+0.1, Yvals, '.b');
    
    consolSame=Yvals;
    
    % ==== diffdir
    afpVals=OUTPUT_multidir.firsttarget.FFrelBaseInBin(i).PBS-OUTPUT_multidir.firsttarget.FFrelBaseInBin(i).MUSC;
    learningVals=OUTPUT_multidir.firsttarget.FFrelBaseInBin(i).PBS;
    Yvals=afpVals./learningVals;
    
    Y2=mean(Yvals);
    Y2err=lt_sem(Yvals);
    plot(X(2)+0.1, Yvals, '.r');
     consolDiff=Yvals;
   
    %=== PLOT
    lt_plot_bar(X(1), [Y1], {'Errors', [Y1err], 'Color','b'});
    lt_plot_bar(X(2), [Y2], {'Errors',  [Y2err], 'Color','r'});
    
    
    % ==== STATS (on unpaired)
    p=ranksum(consolSame, consolDiff);
    
    if p<0.1;
        lt_plot_text(X(1), 1.1, ['unpaired p=' num2str(p)],'k');
    end
    

    % === stats (paired)
    yy1=(Paired_Samedir_PBS-Paired_Samedir_MUSC)./Paired_Samedir_PBS;
    yy2=(Paired_Multdir_PBS - Paired_Multdir_MUSC)./Paired_Multdir_PBS;
    
    p=signrank(yy1, yy2);
    
    if p<0.2;
        lt_plot_text(X(1), 1.2*max(Y1), ['paired p=' num2str(p)],'m');
    end

    
    
end

line(xlim, [1 1], 'Color', 'k', 'LineStyle', '--');


%% ====== plot consolidation across days for both samedir and diff dir on one plot [ALL] [1ST TARG ONLY] [LEARNING VS. CONSOLIDATION]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('learning (PBS)');
ylabel('consolidation');
title('1st target, on day bin to plot');

i=DayBinToUse;

    % === PLOT UNPAIRED MEANS, ETC
    % ==== samedir
    consolVals=OUTPUT_samedir.firsttarget.FFrelBaseInBin(i).MUSC./OUTPUT_samedir.firsttarget.FFrelBaseInBin(i).PBS;
    learningVals=OUTPUT_samedir.firsttarget.FFrelBaseInBin(i).PBS;
    
    plot(learningVals, consolVals, 'ob');
        
    consolSame=consolVals;
    
    % ==== diffdir
    consolVals=OUTPUT_multidir.firsttarget.FFrelBaseInBin(i).MUSC./OUTPUT_multidir.firsttarget.FFrelBaseInBin(i).PBS;
    learningVals=OUTPUT_multidir.firsttarget.FFrelBaseInBin(i).PBS;
    
    plot(learningVals, consolVals, 'or');
        
    consolDiff=consolVals;
   
%     % ==== STATS (on unpaired)
%     p=ranksum(consolSame, consolDiff);
% %     [~, p]=ttest2(consolSame, consolDiff);
%     
%     if p<0.1;
%         lt_plot_text(X(1), 1.1, ['unpaired p=' num2str(p)],'k');
%     end

%% ====== [ALL] [1ST TARG ONLY]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('learning (PBS)');
ylabel('MUSC pitch');
title('1st target, on day bin to plot');
i=DayBinToUse;

    % === PLOT UNPAIRED MEANS, ETC
    % ==== samedir
    muscVals=OUTPUT_samedir.firsttarget.FFrelBaseInBin(i).MUSC;
    learningVals=OUTPUT_samedir.firsttarget.FFrelBaseInBin(i).PBS;
    
    plot(learningVals, muscVals, 'ob');
        
    consolSame=consolVals;
    
    % ==== diffdir
    muscVals=OUTPUT_multidir.firsttarget.FFrelBaseInBin(i).MUSC;
    learningVals=OUTPUT_multidir.firsttarget.FFrelBaseInBin(i).PBS;
    
    plot(learningVals, muscVals, 'or');
        
    consolDiff=consolVals;
   
%     % ==== STATS (on unpaired)
%     p=ranksum(consolSame, consolDiff);
% %     [~, p]=ttest2(consolSame, consolDiff);
%     
%     if p<0.1;
%         lt_plot_text(X(1), 1.1, ['unpaired p=' num2str(p)],'k');
%     end

%% ======  [ALL] [1ST TARG ONLY] [LEARNING VS. CONSOLIDATION]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('learning (PBS)');
ylabel('AFP bias');
title('1st target, on day bin to plot');

% === PLOT UNPAIRED MEANS, ETC
% ==== samedir
muscVals=OUTPUT_samedir.firsttarget.FFrelBaseInBin(i).MUSC;
learningVals=OUTPUT_samedir.firsttarget.FFrelBaseInBin(i).PBS;
afpBias=learningVals-muscVals;

plot(learningVals, afpBias, 'ob');

consolSame=consolVals;

% ==== diffdir
muscVals=OUTPUT_multidir.firsttarget.FFrelBaseInBin(i).MUSC;
learningVals=OUTPUT_multidir.firsttarget.FFrelBaseInBin(i).PBS;
afpBias=learningVals-muscVals;

plot(learningVals, afpBias, 'or');

consolDiff=consolVals;

%     % ==== STATS (on unpaired)
%     p=ranksum(consolSame, consolDiff);
% %     [~, p]=ttest2(consolSame, consolDiff);
%
%     if p<0.1;
%         lt_plot_text(X(1), 1.1, ['unpaired p=' num2str(p)],'k');
%     end

%% ============= COLLECT PAIRED STUFF FOR ANALYSES BELOW
% ==== COLLECT SAMEDIR AND DIFFDIR
Paired_Samedir_PBS=[];
Paired_Samedir_MUSC=[];
Paired_Multdir_PBS=[];
Paired_Multdir_MUSC=[];

i=DayBinToUse;
for j=1:length(OUTPUT_samedir.INFORMATION(i).experimentNum);
    
    exptnumSD=OUTPUT_samedir.INFORMATION(i).experimentNum(j);
    birdnameSD=DATSTRUCT_samedir.INFORMATION(exptnumSD).birdname;
    
    % figure out if it has a multidir counterpart
    for jj=1:length(OUTPUT_multidir.INFORMATION(i).experimentNum);
        
        
        exptnumMD=OUTPUT_multidir.INFORMATION(i).experimentNum(jj);
        birdnameMD=DATSTRUCT_multidir.INFORMATION(exptnumMD).birdname;
        
        % ----- check if birdname and exptname are same
        if strcmp(birdnameSD, birdnameMD)
            
            % ===== COLLECT DATA FOR BOTH
            Paired_Samedir_PBS=[Paired_Samedir_PBS OUTPUT_samedir.firsttarget.FFrelBaseInBin(i).PBS(j)];
            Paired_Samedir_MUSC=[Paired_Samedir_MUSC OUTPUT_samedir.firsttarget.FFrelBaseInBin(i).MUSC(j)];
            Paired_Multdir_PBS=[Paired_Multdir_PBS OUTPUT_multidir.firsttarget.FFrelBaseInBin(i).PBS(jj)];
            Paired_Multdir_MUSC=[Paired_Multdir_MUSC OUTPUT_multidir.firsttarget.FFrelBaseInBin(i).MUSC(jj)];
            disp(birdnameSD)
            break
            
        end
    end
end

%% SAMEDIR AND MULTIDIR, SIGNIFICANT EFFECTS WITHIN THEM? (i.e effect of musc)
% % ==== all paired
% % --- same dir
% [~, p]=ttest(Paired_Samedir_MUSC, Paired_Samedir_PBS)
% [p]=signrank(Paired_Samedir_MUSC, Paired_Samedir_PBS)
% 
% % --- diff dir
% [~, p]=ttest(Paired_Multdir_MUSC, Paired_Multdir_PBS)
% [p]=signrank(Paired_Multdir_MUSC, Paired_Multdir_PBS)
% 
% 
% % === all unpaired
% % -- same dir
% [~, p]=ttest(OUTPUT_samedir.separation(DayBinToUse).MUSC, OUTPUT_samedir.separation(DayBinToUse).PBS)
% [p]=signrank(OUTPUT_samedir.separation(DayBinToUse).MUSC, OUTPUT_samedir.separation(DayBinToUse).PBS)
% 
% % -- diff dir
% [~, p]=ttest(OUTPUT_multidir.separation(DayBinToUse).MUSC, OUTPUT_multidir.separation(DayBinToUse).PBS)
% [p]=signrank(OUTPUT_multidir.separation(DayBinToUse).MUSC, OUTPUT_multidir.separation(DayBinToUse).PBS)

%% == consolidation effect is not explained by difference in magnitude of shift

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('magnitude of shift (mean of absolute value)');
ylabel('consolidation');

% --- samedir experiments
color='b';
Consolidation=OUTPUT_samedir.separation(DayBinToUse).MUSC./OUTPUT_samedir.separation(DayBinToUse).PBS;
ShiftMagnitude=OUTPUT_samedir.separation(DayBinToUse).PBS/2;

plot(ShiftMagnitude, Consolidation, 'o','Color',color)

% -- diff dir 
color='r';
Consolidation=OUTPUT_multidir.separation(DayBinToUse).MUSC./OUTPUT_multidir.separation(DayBinToUse).PBS;
ShiftMagnitude=OUTPUT_multidir.separation(DayBinToUse).PBS/2;

plot(ShiftMagnitude, Consolidation, 'o','Color',color)



[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('paired expts')
xlabel('magnitude of shift (mean of absolute value)');
ylabel('consolidation');

% --- samedir experiments
color='b';
Consolidation=Paired_Samedir_MUSC./Paired_Samedir_PBS;
ShiftMagnitude=Paired_Samedir_PBS/2;

plot(ShiftMagnitude, Consolidation, 'o','Color',color)

% -- diff dir 
color='r';
Consolidation=Paired_Multdir_MUSC./Paired_Multdir_PBS;
ShiftMagnitude=Paired_Multdir_PBS/2;

plot(ShiftMagnitude, Consolidation, 'o','Color',color)


%% is there a difference in magnitude of shift?

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
ylabel('magnitude of shift (mean of absolute value)');
xlabel('samedir -- multidir');
% --- samedir experiments
ShiftMagnitude1=OUTPUT_samedir.separation(DayBinToUse).PBS/2;
% -- diff dir
ShiftMagnitude2=OUTPUT_multidir.separation(DayBinToUse).PBS/2;

X=[1 2];

plot(1.1, ShiftMagnitude1, 'o');
plot(2.1, ShiftMagnitude2, 'o');

lt_plot_bar(X,[mean(ShiftMagnitude1) mean(ShiftMagnitude2)], {'Errors',[lt_sem(ShiftMagnitude1) lt_sem(ShiftMagnitude2)]});

[p]=ranksum(ShiftMagnitude1, ShiftMagnitude2);
lt_plot_text(1.5, max(ShiftMagnitude1)*1.1, num2str(p, '%3.2g'), 'r');


%% is there a difference in magnitude of shift? [first targ]

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('first targ');
ylabel('magnitude of shift (mean of absolute value)');
xlabel('samedir -- multidir');
% --- samedir experiments
ShiftMagnitude1=OUTPUT_samedir.firsttarget.FFrelBaseInBin(DayBinToUse).PBS;
% -- diff dir
ShiftMagnitude2=OUTPUT_multidir.firsttarget.FFrelBaseInBin(DayBinToUse).PBS;

X=[1 2];

plot(1.1, ShiftMagnitude1, 'ok');
plot(2.1, ShiftMagnitude2, 'ok');

lt_plot_bar(X,[mean(ShiftMagnitude1) mean(ShiftMagnitude2)], {'Errors',[lt_sem(ShiftMagnitude1) lt_sem(ShiftMagnitude2)]});

[p]=ranksum(ShiftMagnitude1, ShiftMagnitude2);
lt_plot_text(1.5, max(ShiftMagnitude1)*1.1, num2str(p, '%3.2g'), 'r');



%% ==== does difference in variance during holding (across days) account for difference in consolidation?) - e.g. instability impairs consoldiation

% ++++++++++++++++++++++++++++++++++++++++++ FIRST TARG ONLY
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
ylabel('consolidation');
xlabel('STD of PBS during consol');
title('only first targ');

% === multidir
Consolidation = OUTPUT_multidir.firsttarget.FFrelBaseInBin(DayBinToUse).MUSC./OUTPUT_multidir.firsttarget.FFrelBaseInBin(DayBinToUse).PBS;
ExptNums= OUTPUT_multidir.INFORMATION(DayBinToUse).experimentNum;

STDpbsAll=[];
for i=1:length(ExptNums)
    exptnum=ExptNums(i);
    
    % for this expt, figure out the variance of the pitch for each syllable
    % over the consolidation period
    % -- first targ
    FFpbs=DATSTRUCT_multidir.firsttarget(exptnum).FFminusBase_Mean_PBS;
    consolPeriod=DATSTRUCT_multidir.INFORMATION(exptnum).consolPeriod;
    
    FFpbs=FFpbs(consolPeriod(1):consolPeriod(2));
    
    FF_std=std(FFpbs(~isnan(FFpbs)));
    STDpbsAll=[STDpbsAll FF_std];
end

plot(STDpbsAll, Consolidation, 'or');

% === same dir
Consolidation = OUTPUT_samedir.firsttarget.FFrelBaseInBin(DayBinToUse).MUSC./OUTPUT_samedir.firsttarget.FFrelBaseInBin(DayBinToUse).PBS;
ExptNums= OUTPUT_samedir.INFORMATION(DayBinToUse).experimentNum;

STDpbsAll=[];
for i=1:length(ExptNums)
    exptnum=ExptNums(i);
    
    % for this expt, figure out the variance of the pitch for each syllable
    % over the consolidation period
    % -- first targ
    FFpbs=DATSTRUCT_samedir.firsttarget(exptnum).FFminusBase_Mean_PBS;
    consolPeriod=DATSTRUCT_samedir.INFORMATION(exptnum).consolPeriod;
    
    FFpbs=FFpbs(consolPeriod(1):consolPeriod(2));
    
    FF_std=std(FFpbs(~isnan(FFpbs)));
    STDpbsAll=[STDpbsAll FF_std];
end

plot(STDpbsAll, Consolidation, 'ob');


    
% +++++++++++++++++++++++++++++++ SECOND TARG ONLY
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
ylabel('consolidation');
xlabel('STD of PBS during consol');
title('only second targ');

% === multidir
Consolidation = OUTPUT_multidir.secondtarg.FFrelBaseInBin(DayBinToUse).MUSC./OUTPUT_multidir.secondtarg.FFrelBaseInBin(DayBinToUse).PBS;
ExptNums= OUTPUT_multidir.INFORMATION(DayBinToUse).experimentNum;

STDpbsAll=[];
for i=1:length(ExptNums)
    exptnum=ExptNums(i);
    
    % for this expt, figure out the variance of the pitch for each syllable
    % over the consolidation period
    % -- first targ
    FFpbs=DATSTRUCT_multidir.secondtarget(exptnum).FFminusBase_Mean_PBS;
    consolPeriod=DATSTRUCT_multidir.INFORMATION(exptnum).consolPeriod;
    
    FFpbs=FFpbs(consolPeriod(1):consolPeriod(2));
    
    FF_std=std(FFpbs(~isnan(FFpbs)));
    STDpbsAll=[STDpbsAll FF_std];
end

plot(STDpbsAll, Consolidation, 'or');

% === same dir
Consolidation = OUTPUT_samedir.secondtarg.FFrelBaseInBin(DayBinToUse).MUSC./OUTPUT_samedir.secondtarg.FFrelBaseInBin(DayBinToUse).PBS;
ExptNums= OUTPUT_samedir.INFORMATION(DayBinToUse).experimentNum;

STDpbsAll=[];
for i=1:length(ExptNums)
    exptnum=ExptNums(i);
    
    % for this expt, figure out the variance of the pitch for each syllable
    % over the consolidation period
    % -- first targ
    FFpbs=DATSTRUCT_samedir.secondtarget(exptnum).FFminusBase_Mean_PBS;
    consolPeriod=DATSTRUCT_samedir.INFORMATION(exptnum).consolPeriod;
    
    FFpbs=FFpbs(consolPeriod(1):consolPeriod(2));
    
    FF_std=std(FFpbs(~isnan(FFpbs)));
    STDpbsAll=[STDpbsAll FF_std];
end

plot(STDpbsAll, Consolidation, 'ob');


% +++++++++++++++++++++++++++++++ MEAN OF BOTH TARGETS (FIRST GET STD, THEN
% TAKE MEAN)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
ylabel('consolidation');
xlabel('STD of PBS during consol');
title('both targs (first take std, then mean)');

% === multidir
Consolidation = OUTPUT_multidir.separation(DayBinToUse).MUSC./OUTPUT_multidir.separation(DayBinToUse).PBS;
ExptNums= OUTPUT_multidir.INFORMATION(DayBinToUse).experimentNum;

STDpbsAll=[];
for i=1:length(ExptNums)
    exptnum=ExptNums(i);
    
    % for this expt, figure out the variance of the pitch for each syllable
    % over the consolidation period
    consolPeriod=DATSTRUCT_multidir.INFORMATION(exptnum).consolPeriod;
    % -- first targ
    FFpbs=DATSTRUCT_multidir.firsttarget(exptnum).FFminusBase_Mean_PBS;
    FFpbs=FFpbs(consolPeriod(1):consolPeriod(2));
    FF_std1=std(FFpbs(~isnan(FFpbs)));

    % --- second targ
    FFpbs=DATSTRUCT_multidir.secondtarget(exptnum).FFminusBase_Mean_PBS;
    FFpbs=FFpbs(consolPeriod(1):consolPeriod(2));
    FF_std2=std(FFpbs(~isnan(FFpbs)));
    
    % --- meanlt_plot_text(STDpbsAll, Consolidation+0.05, BirdnameAll, 'r');

    FF_std=mean([FF_std1, FF_std2]);
    STDpbsAll=[STDpbsAll FF_std];
    
end

plot(STDpbsAll, Consolidation, 'or');


% === samedir
Consolidation = OUTPUT_samedir.separation(DayBinToUse).MUSC./OUTPUT_samedir.separation(DayBinToUse).PBS;
ExptNums= OUTPUT_samedir.INFORMATION(DayBinToUse).experimentNum;

STDpbsAll=[];
for i=1:length(ExptNums)
    exptnum=ExptNums(i);
    
    % for this expt, figure out the variance of the pitch for each syllable
    % over the consolidation period
    consolPeriod=DATSTRUCT_samedir.INFORMATION(exptnum).consolPeriod;
    % -- first targ
    FFpbs=DATSTRUCT_samedir.firsttarget(exptnum).FFminusBase_Mean_PBS;
    FFpbs=FFpbs(consolPeriod(1):consolPeriod(2));
    FF_std1=std(FFpbs(~isnan(FFpbs)));

    % --- second targ
    FFpbs=DATSTRUCT_samedir.secondtarget(exptnum).FFminusBase_Mean_PBS;
    FFpbs=FFpbs(consolPeriod(1):consolPeriod(2));
    FF_std2=std(FFpbs(~isnan(FFpbs)));
    
    % --- mean
    FF_std=mean([FF_std1, FF_std2]);
    STDpbsAll=[STDpbsAll FF_std];
    
end

plot(STDpbsAll, Consolidation, 'ob');


% +++++++++++++++++++++++++++++++ MEAN OF BOTH TARGETS (FIRST GET STD, THEN
% TAKE MEAN)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
ylabel('consolidation');
xlabel('STD of PBS during consol');
title('both targs (first take std, then mean)');

% === multidir
Consolidation = OUTPUT_multidir.separation(DayBinToUse).MUSC./OUTPUT_multidir.separation(DayBinToUse).PBS;
ExptNums= OUTPUT_multidir.INFORMATION(DayBinToUse).experimentNum;

STDpbsAll=[];
BirdnameAll={};
for i=1:length(ExptNums)
    exptnum=ExptNums(i);
    
    % for this expt, figure out the variance of the pitch for each syllable
    % over the consolidation period
    consolPeriod=DATSTRUCT_multidir.INFORMATION(exptnum).consolPeriod;
    % -- first targ
    FFpbs=DATSTRUCT_multidir.firsttarget(exptnum).FFminusBase_Mean_PBS;
    FFpbs=FFpbs(consolPeriod(1):consolPeriod(2));
    FF_std1=std(FFpbs(~isnan(FFpbs)));

    % --- second targ
    FFpbs=DATSTRUCT_multidir.secondtarget(exptnum).FFminusBase_Mean_PBS;
    FFpbs=FFpbs(consolPeriod(1):consolPeriod(2));
    FF_std2=std(FFpbs(~isnan(FFpbs)));
    
    % --- meanlt_plot_text(STDpbsAll, Consolidation+0.05, BirdnameAll, 'r');

    FF_std=mean([FF_std1, FF_std2]);
    STDpbsAll=[STDpbsAll FF_std];
    
    % --- collect exptname and birdname
    birdname=DATSTRUCT_multidir.INFORMATION(exptnum).birdname;
    BirdnameAll=[BirdnameAll birdname];
end

plot(STDpbsAll, Consolidation, 'or');
lt_plot_text(STDpbsAll, Consolidation+0.05, BirdnameAll, 'r');


% === samedir
Consolidation = OUTPUT_samedir.separation(DayBinToUse).MUSC./OUTPUT_samedir.separation(DayBinToUse).PBS;
ExptNums= OUTPUT_samedir.INFORMATION(DayBinToUse).experimentNum;

STDpbsAll=[];
BirdnameAll={};
for i=1:length(ExptNums)
    exptnum=ExptNums(i);
    
    % for this expt, figure out the variance of the pitch for each syllable
    % over the consolidation period
    consolPeriod=DATSTRUCT_samedir.INFORMATION(exptnum).consolPeriod;
    % -- first targ
    FFpbs=DATSTRUCT_samedir.firsttarget(exptnum).FFminusBase_Mean_PBS;
    FFpbs=FFpbs(consolPeriod(1):consolPeriod(2));
    FF_std1=std(FFpbs(~isnan(FFpbs)));

    % --- second targ
    FFpbs=DATSTRUCT_samedir.secondtarget(exptnum).FFminusBase_Mean_PBS;
    FFpbs=FFpbs(consolPeriod(1):consolPeriod(2));
    FF_std2=std(FFpbs(~isnan(FFpbs)));
    
    % --- mean
    FF_std=mean([FF_std1, FF_std2]);
    STDpbsAll=[STDpbsAll FF_std];
    
    % --- collect exptname and birdname
    birdname=DATSTRUCT_samedir.INFORMATION(exptnum).birdname;
    BirdnameAll=[BirdnameAll birdname];
    
end

plot(STDpbsAll, Consolidation, 'ob');
lt_plot_text(STDpbsAll, Consolidation+0.05, BirdnameAll, 'b');


%% ==== Plot end consol vs. start consol in consolidation vs. starting consolidation
EndDayBin=DayBinToUse;
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
ylabel('end consolid');
xlabel('start consolid');

% ==== MULTIDIR
StartConsolid=OUTPUT_multidir.separation(1).MUSC./OUTPUT_multidir.separation(1).PBS;
StartExptnums=OUTPUT_multidir.INFORMATION(1).experimentNum;

EndConsolid=OUTPUT_multidir.separation(EndDayBin).MUSC./OUTPUT_multidir.separation(EndDayBin).PBS;
EndExptnums=OUTPUT_multidir.INFORMATION(EndDayBin).experimentNum;

% --- plot, go through each potential experiment number
for i=1:length(StartExptnums);
    exptnum=StartExptnums(i);
        
    ii=find(EndExptnums==exptnum);
    
    if ~isempty(ii);
    % then plot one dot
    plot(StartConsolid(i), EndConsolid(ii), 'or');
    
    end
       
    
end


% ==== SAMEDIR
StartConsolid=OUTPUT_samedir.separation(1).MUSC./OUTPUT_samedir.separation(1).PBS;
StartExptnums=OUTPUT_samedir.INFORMATION(1).experimentNum;

EndConsolid=OUTPUT_samedir.separation(EndDayBin).MUSC./OUTPUT_samedir.separation(EndDayBin).PBS;
EndExptnums=OUTPUT_samedir.INFORMATION(EndDayBin).experimentNum;

% --- plot, go through each potential experiment number
for i=1:length(StartExptnums);
    exptnum=StartExptnums(i);
        
    ii=find(EndExptnums==exptnum);
    
    if ~isempty(ii);
    % then plot one dot
    plot(StartConsolid(i), EndConsolid(ii), 'ob');
    
    end
       
    
end

line([0 1], [0 1]);


%% plot change in consol vs. start consol.
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('multidir');
ylabel('end minus start consolid');
xlabel('start consolid');

% ==== MULTIDIR
StartConsolid=OUTPUT_multidir.separation(1).MUSC./OUTPUT_multidir.separation(1).PBS;
StartExptnums=OUTPUT_multidir.INFORMATION(1).experimentNum;

EndConsolid=OUTPUT_multidir.separation(EndDayBin).MUSC./OUTPUT_multidir.separation(EndDayBin).PBS;
EndExptnums=OUTPUT_multidir.INFORMATION(EndDayBin).experimentNum;

% --- plot, go through each potential experiment number
X=[];
Y=[];
for i=1:length(StartExptnums);
    exptnum=StartExptnums(i);
        
    ii=find(EndExptnums==exptnum);
    
    if ~isempty(ii);
    % then plot one dot
    plot(StartConsolid(i), EndConsolid(ii)-StartConsolid(i), 'or');
    
    % save for correlation analysis
    X=[X StartConsolid(i)];
    Y=[Y EndConsolid(ii)-StartConsolid(i)];
    end
    
end
lt_regress(Y, X, 1, 0, 1, 1, 'r');


% ==== SAMEDIR
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('samedir');
ylabel('end minus start consolid');
xlabel('start consolid');

StartConsolid=OUTPUT_samedir.separation(1).MUSC./OUTPUT_samedir.separation(1).PBS;
StartExptnums=OUTPUT_samedir.INFORMATION(1).experimentNum;

EndConsolid=OUTPUT_samedir.separation(EndDayBin).MUSC./OUTPUT_samedir.separation(EndDayBin).PBS;
EndExptnums=OUTPUT_samedir.INFORMATION(EndDayBin).experimentNum;

% --- plot, go through each potential experiment number
X=[];
Y=[];
for i=1:length(StartExptnums);
    exptnum=StartExptnums(i);
        
    ii=find(EndExptnums==exptnum);
    
    if ~isempty(ii);
    % then plot one dot
    plot(StartConsolid(i), EndConsolid(ii)-StartConsolid(i), 'ob');
    
        % save for correlation analysis
    X=[X StartConsolid(i)];
    Y=[Y EndConsolid(ii)-StartConsolid(i)];

    end
       
    
end
lt_regress(Y, X, 1, 0, 1, 1, 'b');


line(xlim, [0 0]);
line([0 0], ylim);


%% ==== is there an increase in consolidation over days? (for same-dir and multidir)

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('change in consolidation');

% ==== MULTIDIR
StartConsolid=OUTPUT_multidir.separation(1).MUSC./OUTPUT_multidir.separation(1).PBS;
StartExptnums=OUTPUT_multidir.INFORMATION(1).experimentNum;

EndConsolid=OUTPUT_multidir.separation(EndDayBin).MUSC./OUTPUT_multidir.separation(EndDayBin).PBS;
EndExptnums=OUTPUT_multidir.INFORMATION(EndDayBin).experimentNum;

% --- plot, go through each potential experiment number
X=[];
Y=[];
for i=1:length(StartExptnums);
    exptnum=StartExptnums(i);
        
    ii=find(EndExptnums==exptnum);
    
    if ~isempty(ii);
    
    % save for correlation analysis
    X=[X StartConsolid(i)];
    Y=[Y EndConsolid(ii)];
    end
    
end

% -- plot
plot([1.1 2.1], [X; Y], '-or')
lt_plot_bar([1 2], [mean(X), mean(Y)], {'Errors',[lt_sem(X) lt_sem(Y)], 'Color', 'r'});

p=signrank(X, Y);
lt_plot_pvalue(p, 'multidir', 2);

% ==== SAMEDIR
StartConsolid=OUTPUT_samedir.separation(1).MUSC./OUTPUT_samedir.separation(1).PBS;
StartExptnums=OUTPUT_samedir.INFORMATION(1).experimentNum;

EndConsolid=OUTPUT_samedir.separation(EndDayBin).MUSC./OUTPUT_samedir.separation(EndDayBin).PBS;
EndExptnums=OUTPUT_samedir.INFORMATION(EndDayBin).experimentNum;

% --- plot, go through each potential experiment number
X=[];
Y=[];
for i=1:length(StartExptnums);
    exptnum=StartExptnums(i);
        
    ii=find(EndExptnums==exptnum);
    
    if ~isempty(ii);
    
    % save for correlation analysis
    X=[X StartConsolid(i)];
    Y=[Y EndConsolid(ii)];
    end
    
end

% -- plot
plot([4.1 5.1], [X; Y], '-ob')
lt_plot_bar([4.1 5.1], [mean(X), mean(Y)], {'Errors',[lt_sem(X) lt_sem(Y)], 'Color', 'b'});

p=signrank(X, Y);
lt_plot_pvalue(p, 'samedir', 1);


%% ==== is there an increase in consolidation over days? (for same-dir and multidir)

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
ylabel('change in consolidation');
title('first targ');

% ==== MULTIDIR
StartConsolid=OUTPUT_multidir.firsttarget.FFrelBaseInBin(1).MUSC./OUTPUT_multidir.firsttarget.FFrelBaseInBin(1).PBS;
StartExptnums=OUTPUT_multidir.INFORMATION(1).experimentNum;

EndConsolid=OUTPUT_multidir.firsttarget.FFrelBaseInBin(EndDayBin).MUSC./OUTPUT_multidir.firsttarget.FFrelBaseInBin(EndDayBin).PBS;
EndExptnums=OUTPUT_multidir.INFORMATION(EndDayBin).experimentNum;

% --- plot, go through each potential experiment number
X=[];
Y=[];
for i=1:length(StartExptnums);
    exptnum=StartExptnums(i);
        
    ii=find(EndExptnums==exptnum);
    
    if ~isempty(ii);
    
    % save for correlation analysis
    X=[X StartConsolid(i)];
    Y=[Y EndConsolid(ii)];
    end
    
end

% -- plot
plot([1.1 2.1], [X; Y], '-or')
lt_plot_bar([1 2], [mean(X), mean(Y)], {'Errors',[lt_sem(X) lt_sem(Y)], 'Color', 'r'});

p=signrank(X, Y);
lt_plot_pvalue(p, 'multidir', 2);

% ==== SAMEDIR
StartConsolid=OUTPUT_samedir.firsttarget.FFrelBaseInBin(1).MUSC./OUTPUT_samedir.firsttarget.FFrelBaseInBin(1).PBS;
StartExptnums=OUTPUT_samedir.INFORMATION(1).experimentNum;

EndConsolid=OUTPUT_samedir.firsttarget.FFrelBaseInBin(EndDayBin).MUSC./OUTPUT_samedir.firsttarget.FFrelBaseInBin(EndDayBin).PBS;
EndExptnums=OUTPUT_samedir.INFORMATION(EndDayBin).experimentNum;

% --- plot, go through each potential experiment number
X=[];
Y=[];
for i=1:length(StartExptnums);
    exptnum=StartExptnums(i);
        
    ii=find(EndExptnums==exptnum);
    
    if ~isempty(ii);
    
    % save for correlation analysis
    X=[X StartConsolid(i)];
    Y=[Y EndConsolid(ii)];
    end
    
end

% -- plot
plot([4.1 5.1], [X; Y], '-ob')
lt_plot_bar([4.1 5.1], [mean(X), mean(Y)], {'Errors',[lt_sem(X) lt_sem(Y)], 'Color', 'b'});

p=signrank(X, Y);
lt_plot_pvalue(p, 'samedir', 1);


%% ===== is there correlation between time froms tart of entire experiments and magnitude of consolidation?
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('days from entire expt start (1st consol day)');
ylabel('consolidation (end)');

% ====== DIFF DIR
Consolidation=OUTPUT_multidir.separation(DayBinToUse).MUSC./OUTPUT_multidir.separation(DayBinToUse).PBS;
ExptNums=OUTPUT_multidir.INFORMATION(DayBinToUse).experimentNum;

% figure out the time from onset of entire experiment to first
% consoldiation day for each of these experiments
DaysFromExptStart=[];
for i=1:length(ExptNums);
   
   exptnum=ExptNums(i);
   
   day=DATSTRUCT_multidir.INFORMATION(exptnum).day1_FromStartWN + DATSTRUCT_multidir.INFORMATION(exptnum).consolPeriod(1) -1;
           
    DaysFromExptStart=[DaysFromExptStart day];
end

plot(DaysFromExptStart, Consolidation, 'or');


% ====== SAME DIR
Consolidation=OUTPUT_samedir.separation(DayBinToUse).MUSC./OUTPUT_samedir.separation(DayBinToUse).PBS;
ExptNums=OUTPUT_samedir.INFORMATION(DayBinToUse).experimentNum;

% figure out the time from onset of entire experiment to first
% consoldiation day for each of these experiments
DaysFromExptStart=[];
for i=1:length(ExptNums);
   
   exptnum=ExptNums(i);
   
   day=DATSTRUCT_samedir.INFORMATION(exptnum).day1_FromStartWN + DATSTRUCT_samedir.INFORMATION(exptnum).consolPeriod(1) -1;
           
    DaysFromExptStart=[DaysFromExptStart day];
end

plot(DaysFromExptStart, Consolidation, 'ob');



lt_plot_zeroline;





%% ===== is there correlation between time froms tart of entire experiments and magnitude of consolidation?
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('days from entire expt start (1st consol day)');
ylabel('consolidation (end)');
title('first targ only');

% ====== DIFF DIR
Consolidation=OUTPUT_multidir.firsttarget.FFrelBaseInBin(DayBinToUse).MUSC./OUTPUT_multidir.firsttarget.FFrelBaseInBin(DayBinToUse).PBS;
ExptNums=OUTPUT_multidir.INFORMATION(DayBinToUse).experimentNum;

% figure out the time from onset of entire experiment to first
% consoldiation day for each of these experiments
DaysFromExptStart=[];
for i=1:length(ExptNums);
   
   exptnum=ExptNums(i);
   
   day=DATSTRUCT_multidir.INFORMATION(exptnum).day1_FromStartWN + DATSTRUCT_multidir.INFORMATION(exptnum).consolPeriod(1) -1;
           
    DaysFromExptStart=[DaysFromExptStart day];
end

plot(DaysFromExptStart, Consolidation, 'or');


% ====== SAME DIR
Consolidation=OUTPUT_samedir.firsttarget.FFrelBaseInBin(DayBinToUse).MUSC./OUTPUT_samedir.firsttarget.FFrelBaseInBin(DayBinToUse).PBS;
ExptNums=OUTPUT_samedir.INFORMATION(DayBinToUse).experimentNum;

% figure out the time from onset of entire experiment to first
% consoldiation day for each of these experiments
DaysFromExptStart=[];
for i=1:length(ExptNums);
   
   exptnum=ExptNums(i);
   
   day=DATSTRUCT_samedir.INFORMATION(exptnum).day1_FromStartWN + DATSTRUCT_samedir.INFORMATION(exptnum).consolPeriod(1) -1;
           
    DaysFromExptStart=[DaysFromExptStart day];
end

plot(DaysFromExptStart, Consolidation, 'ob');



lt_plot_zeroline;

%% ==== duration from start of expt diff?

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
ylabel('days from entire expt start (1st consol day)');
title('first targ only');

% ====== DIFF DIR
ExptNums=OUTPUT_multidir.INFORMATION(DayBinToUse).experimentNum;

% figure out the time from onset of entire experiment to first
% consoldiation day for each of these experiments
DaysFromExptStart=[];
for i=1:length(ExptNums);
   
   exptnum=ExptNums(i);
   
   day=DATSTRUCT_multidir.INFORMATION(exptnum).day1_FromStartWN + DATSTRUCT_multidir.INFORMATION(exptnum).consolPeriod(1) -1;
           
    DaysFromExptStart=[DaysFromExptStart day];
end

DaysFromExptStart_diff=DaysFromExptStart;

% ====== SAME DIR
ExptNums=OUTPUT_samedir.INFORMATION(DayBinToUse).experimentNum;

% figure out the time from onset of entire experiment to first
% consoldiation day for each of these experiments
DaysFromExptStart=[];
for i=1:length(ExptNums);
   
   exptnum=ExptNums(i);
   
   day=DATSTRUCT_samedir.INFORMATION(exptnum).day1_FromStartWN + DATSTRUCT_samedir.INFORMATION(exptnum).consolPeriod(1) -1;
           
    DaysFromExptStart=[DaysFromExptStart day];
end
DaysFromExptStart_same=DaysFromExptStart;


%$ ===== plot
plot(1, DaysFromExptStart_diff, 'or');
plot(2, DaysFromExptStart_same, 'ob');


lt_plot_bar([1 2], [mean(DaysFromExptStart_diff), mean(DaysFromExptStart_same)], {'Errors', [lt_sem(DaysFromExptStart_diff) lt_sem(DaysFromExptStart_same)]});

% -- sign test
p=ranksum(DaysFromExptStart_diff, DaysFromExptStart_same);
lt_plot_pvalue(p, 'ranksum',1);

xlim([0 3]);




%%  ======================== ANCOVA, account for potential covariation of consolidation with time from expt start
% FIRST TARG
DaysALL=[];
ConsolidationALL=[];
SameDirALL=[];


% ====== DIFF DIR
Consolidation=OUTPUT_multidir.firsttarget.FFrelBaseInBin(DayBinToUse).MUSC./OUTPUT_multidir.firsttarget.FFrelBaseInBin(DayBinToUse).PBS;
ExptNums=OUTPUT_multidir.INFORMATION(DayBinToUse).experimentNum;

% figure out the time from onset of entire experiment to first
% consoldiation day for each of these experiments
DaysFromExptStart=[];
for i=1:length(ExptNums);
   
   exptnum=ExptNums(i);
   
   day=DATSTRUCT_multidir.INFORMATION(exptnum).day1_FromStartWN + DATSTRUCT_multidir.INFORMATION(exptnum).consolPeriod(1) -1;
           
    DaysFromExptStart=[DaysFromExptStart day];
end

plot(DaysFromExptStart, Consolidation, 'or');

DaysALL=[DaysALL DaysFromExptStart];
ConsolidationALL=[ConsolidationALL Consolidation];
SameDirALL=[SameDirALL 0*ones(1, length(Consolidation))];


% ====== SAME DIR
Consolidation=OUTPUT_samedir.firsttarget.FFrelBaseInBin(DayBinToUse).MUSC./OUTPUT_samedir.firsttarget.FFrelBaseInBin(DayBinToUse).PBS;
ExptNums=OUTPUT_samedir.INFORMATION(DayBinToUse).experimentNum;

% figure out the time from onset of entire experiment to first
% consoldiation day for each of these experiments
DaysFromExptStart=[];
for i=1:length(ExptNums);
   
   exptnum=ExptNums(i);
   
   day=DATSTRUCT_samedir.INFORMATION(exptnum).day1_FromStartWN + DATSTRUCT_samedir.INFORMATION(exptnum).consolPeriod(1) -1;
           
    DaysFromExptStart=[DaysFromExptStart day];
end

plot(DaysFromExptStart, Consolidation, 'ob');

DaysALL=[DaysALL DaysFromExptStart];
ConsolidationALL=[ConsolidationALL Consolidation];
SameDirALL=[SameDirALL 1*ones(1, length(Consolidation))];


aoctool(DaysALL, ConsolidationALL, SameDirALL, 0.05, 'days', 'consol', 'congruent?', 'on', 'parallel lines');


%% ===== is there correlation between time froms tart of that specific phase and magnitude of consolidation?
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('days from start of phase (to 1st consol day)');
ylabel('consolidation (end)');

% ====== DIFF DIR
Consolidation=OUTPUT_multidir.separation(DayBinToUse).MUSC./OUTPUT_multidir.separation(DayBinToUse).PBS;
ExptNums=OUTPUT_multidir.INFORMATION(DayBinToUse).experimentNum;

% figure out the time from onset of entire experiment to first
% consoldiation day for each of these experiments
DaysFromExptStart=[];
for i=1:length(ExptNums);
   
   exptnum=ExptNums(i);
   
   day=DATSTRUCT_multidir.INFORMATION(exptnum).consolPeriod(1);
           
    DaysFromExptStart=[DaysFromExptStart day];
end

plot(DaysFromExptStart, Consolidation, 'or');


% ====== SAME DIR
Consolidation=OUTPUT_samedir.separation(DayBinToUse).MUSC./OUTPUT_samedir.separation(DayBinToUse).PBS;
ExptNums=OUTPUT_samedir.INFORMATION(DayBinToUse).experimentNum;

% figure out the time from onset of entire experiment to first
% consoldiation day for each of these experiments
DaysFromExptStart=[];
for i=1:length(ExptNums);
   
   exptnum=ExptNums(i);
   
   day=DATSTRUCT_samedir.INFORMATION(exptnum).consolPeriod(1);
           
    DaysFromExptStart=[DaysFromExptStart day];
end

plot(DaysFromExptStart, Consolidation, 'ob');



lt_plot_zeroline;





%% ===== is there correlation between time froms tart of that specific phase and magnitude of consolidation?
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('days from start of phase (to 1st consol day)');
ylabel('consolidation (end)');
title('first targ');

% ====== DIFF DIR
Consolidation=OUTPUT_multidir.firsttarget.FFrelBaseInBin(DayBinToUse).MUSC./OUTPUT_multidir.firsttarget.FFrelBaseInBin(DayBinToUse).PBS;
ExptNums=OUTPUT_multidir.INFORMATION(DayBinToUse).experimentNum;

% figure out the time from onset of entire experiment to first
% consoldiation day for each of these experiments
DaysFromExptStart=[];
for i=1:length(ExptNums);
   
   exptnum=ExptNums(i);
   
   day=DATSTRUCT_multidir.INFORMATION(exptnum).consolPeriod(1);
           
    DaysFromExptStart=[DaysFromExptStart day];
end

plot(DaysFromExptStart, Consolidation, 'or');


% ====== SAME DIR
Consolidation=OUTPUT_samedir.firsttarget.FFrelBaseInBin(DayBinToUse).MUSC./OUTPUT_samedir.firsttarget.FFrelBaseInBin(DayBinToUse).PBS;
ExptNums=OUTPUT_samedir.INFORMATION(DayBinToUse).experimentNum;

% figure out the time from onset of entire experiment to first
% consoldiation day for each of these experiments
DaysFromExptStart=[];
for i=1:length(ExptNums);
   
   exptnum=ExptNums(i);
   
   day=DATSTRUCT_samedir.INFORMATION(exptnum).consolPeriod(1);
           
    DaysFromExptStart=[DaysFromExptStart day];
end

plot(DaysFromExptStart, Consolidation, 'ob');



lt_plot_zeroline;






