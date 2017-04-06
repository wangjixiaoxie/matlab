function lt_neural_v2_ANALY_VocModel_anova(VOCALSTRUCTall)
%% for each time bin, predicts fr based on sequence information
% 

do2by2 = 0; % ad hoc, to pull out 2x2 matrix. add interaction term if do this. (omnly for bk7)
logtransform = 1; % log trqansform fr (add 1 before transform)

%% ad hoc, pull out a 2x2 matrix for anova.


if do2by2 ==1
    
    PresylPool = {'p','h'};
    SylPool = {'p','r'};
    DesiredBird = 'bk7';
    useinteraction=1;
    
    for z = 1:Numbirds
        
        if ~strcmp(SummaryStruct.birds(z).birdname, DesiredBird);
            continue
        end
        
        Numneurons = length(SummaryStruct.birds(z).neurons);
        for zz = 1:Numneurons
            
            Vocalstruct = VOCALSTRUCTall.birds(z).neurons(zz).dat;
            
            inds = [];
            for j=1:length(Vocalstruct.data_vocalrend)
                
                if any(strcmp(Vocalstruct.data_vocalrend(j).syl, SylPool)) ...
                        & any(strcmp(Vocalstruct.data_vocalrend(j).presyl, PresylPool));
                    inds = [inds j];
                end
            end
            
            % resave output
            VOCALSTRUCTall.birds(z).neurons(zz).dat.data_vocalrend = Vocalstruct.data_vocalrend(inds);
        end
    end
else
    useinteraction=0;
end



%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++ ANOVA
%% ==================== EXTRACT DATA
% useinteraction =1; use for if do 2 x 2 above.
if do2by2==1
    Groupnames = {'Syl', 'Presyl'}; % adds interaction, and pulls out stats assuming there's only 2 factors
else
    Groupnames = {'Syl', 'Presyl', 'Postsyl'};
end
numtimebins = length(VOCALSTRUCTall.global.bincenters);
ANOVAOUTPUT = struct;

% lt_figure; hold on;
% title('anova rel to syl onset');
Numbirds = length(VOCALSTRUCTall.birds);

for z=1:Numbirds
    Numneurons = length(VOCALSTRUCTall.birds(z).neurons);
    %     z=1
    for zz=1:Numneurons
        VOCALSTRUCT = VOCALSTRUCTall.birds(z).neurons(zz).dat;
        %         zz=12
        for i=1:numtimebins
            
            % ==========
            tmp = {VOCALSTRUCT.data_vocalrend.Spk_bincounts};
            SpkCntMat = cell2mat(tmp');
            SpkCnts = SpkCntMat(:, i);
            
            if logtransform==1
                SpkCnts = log10(SpkCnts+1);
            end
            
            Syl = [VOCALSTRUCT.data_vocalrend.syl]';
            Presyl = [VOCALSTRUCT.data_vocalrend.presyl]';
            Postsyl = [VOCALSTRUCT.data_vocalrend.postsyl]';
            
            % ------ 1) LME model
            if (0)
                modelform = 'SpkCnts ~ 1 + Syl + Presyl';
                X = table(SpkCnts, Syl, Presyl);
                lme = fitlme(X, modelform);
            end
            
            % ----- 2) ANOVA
            Group = {};
            for j=1:length(Groupnames)
                Group{j} = eval(Groupnames{j});
            end
            if useinteraction ==1;
                [p, tabletmp] = anovan(SpkCnts, Group, 'display', 'off', ...
                    'model', 'interaction');
            else
                [p, tabletmp] = anovan(SpkCnts, Group, 'display', 'off');
            end
            
            numgroups = length(Group)+useinteraction;
            ss_total = tabletmp{3+numgroups,2};
            ms_error = tabletmp{2+numgroups,5};
            
            Eta2_all = [];
            Omega2_all = [];
            
            for j=1:numgroups
                ss_effect = tabletmp{1+j,2};
                df_effect = tabletmp{1+j,3};
                
                numerator = ss_effect - df_effect*ms_error;
                denominator = ss_total + ms_error;
                omega2 = numerator/denominator;
                eta2 = ss_effect/ss_total;
                
                Eta2_all = [Eta2_all eta2];
                Omega2_all = [Omega2_all omega2];
            end
            
            ANOVAOUTPUT.timebin(i).Eta2_all = Eta2_all;
            ANOVAOUTPUT.timebin(i).Omega2_all = Omega2_all;
        end
        ANOVAOUTPUT.global.Groupnames = Groupnames;
        if useinteraction ==1
            ANOVAOUTPUT.global.Groupnames = [Groupnames 'interaction'];
        end
        VOCALSTRUCTall.birds(z).neurons(zz).anova.output = ANOVAOUTPUT;
        
        if (0)
            % =================== PLOT ANOVA RESULTS
            tmp = {ANOVAOUTPUT.timebin.Omega2_all};
            OmegaAcrossBins = cell2mat(tmp');
            
            numgroups = length(ANOVAOUTPUT.global.Groupnames);
            plotcols = lt_make_plot_colors(length(Groupnames), 0,0);
            for i=1:numgroups
                groupname = ANOVAOUTPUT.global.Groupnames{i};
                
                omegavals = OmegaAcrossBins(:,i);
                xvals = VOCALSTRUCTall.global.bincenters;
                
                plot(xvals, omegavals, '-', 'Color' ,plotcols{i});
            end
        end
    end
end



%% ================== PLOT ANOVA RESULTS
figcount=1;
subplotrows=3;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

numgroups = length(ANOVAOUTPUT.global.Groupnames);
plotcols = lt_make_plot_colors(length(ANOVAOUTPUT.global.Groupnames), 0,0);


% ============ FOR A SINGLE NEURON
if (0)
    z=1;
    zz=12;
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    birdname = SummaryStruct.birds(z).birdname;
    xlabel('timebin');
    ylabel('omega^2');
    
    ANOVAOUTPUT = VOCALSTRUCTall.birds(z).neurons(zz).anova.output;
    for g =1:numgroups
        groupname = Groupnames{g};
        
        % =================== PLOT ANOVA RESULTS
        tmp = {ANOVAOUTPUT.timebin.Omega2_all};
        OmegaAcrossBins = cell2mat(tmp');
        
        omegavals = OmegaAcrossBins(:,g);
        xvals = VOCALSTRUCTall.global.bincenters;
        
        plot(xvals, omegavals, '-', 'Color' ,plotcols{g});
    end
    lt_plot_zeroline;
    line([VOCALSTRUCTall.global.pretime VOCALSTRUCTall.global.pretime], ylim, 'Color', 'k');
    ylim([-0.05 0.2]);
end


% ======== ONE FOR EACH BIRD
for z=1:Numbirds
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    birdname = SummaryStruct.birds(z).birdname;
    title(['bird: ' birdname]);
    xlabel('timebin');
    ylabel('omega^2');
    
    Numneurons = length(VOCALSTRUCTall.birds(z).neurons);
    for zz=1:Numneurons
        ANOVAOUTPUT = VOCALSTRUCTall.birds(z).neurons(zz).anova.output;
        for g =1:numgroups
            groupname = ANOVAOUTPUT.global.Groupnames{g};
            
            % =================== PLOT ANOVA RESULTS
            tmp = {ANOVAOUTPUT.timebin.Omega2_all};
            OmegaAcrossBins = cell2mat(tmp');
            
            omegavals = OmegaAcrossBins(:,g);
            xvals = VOCALSTRUCTall.global.bincenters;
            
            plot(xvals, omegavals, '-', 'Color' ,plotcols{g});
        end
    end
    lt_plot_zeroline;
    line([VOCALSTRUCTall.global.pretime VOCALSTRUCTall.global.pretime], ylim, 'Color', 'k');
    ylim([-0.05 0.2]);
end


% ==== 1) ONE FOR EACH GROUP, EACH SYL ONE LINE

for g =1 :numgroups
    groupname = ANOVAOUTPUT.global.Groupnames{g};
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(groupname);
    xlabel('timebin');
    ylabel('omega^2');
    
    for z=1:Numbirds
        Numneurons = length(VOCALSTRUCTall.birds(z).neurons);
        
        for zz=1:Numneurons
            ANOVAOUTPUT = VOCALSTRUCTall.birds(z).neurons(zz).anova.output;
            
            % =================== PLOT ANOVA RESULTS
            tmp = {ANOVAOUTPUT.timebin.Omega2_all};
            OmegaAcrossBins = cell2mat(tmp');
            
            omegavals = OmegaAcrossBins(:,g);
            xvals = VOCALSTRUCTall.global.bincenters;
            
            plot(xvals, omegavals, '-', 'Color' ,plotcols{g});
        end
    end
    lt_plot_zeroline;
    line([VOCALSTRUCTall.global.pretime VOCALSTRUCTall.global.pretime], ylim, 'Color', 'k');
    ylim([-0.05 0.2]);
end


% ======= 2) combine all in one plot
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title(groupname);
xlabel('timebin');
ylabel('omega^2');

for g =1 :numgroups
    groupname = ANOVAOUTPUT.global.Groupnames{g};
    
    Omega_all = [];
    for z=1:Numbirds
        Numneurons = length(VOCALSTRUCTall.birds(z).neurons);
        
        for zz=1:Numneurons
            ANOVAOUTPUT = VOCALSTRUCTall.birds(z).neurons(zz).anova.output;
            
            % =================== PLOT ANOVA RESULTS
            tmp = {ANOVAOUTPUT.timebin.Omega2_all};
            OmegaAcrossBins = cell2mat(tmp');
            
            omegavals = OmegaAcrossBins(:,g);
            
            Omega_all = [Omega_all omegavals];
        end
    end
    
    x = VOCALSTRUCTall.global.bincenters;
    ymean = mean(Omega_all');
    ystd = std(Omega_all', 0,1);
    ysem = lt_sem(Omega_all');
    shadedErrorBar(x, ymean, ysem, {'Color', plotcols{g}},1);
end

lt_plot_zeroline;
% legend(gca, Groupnames);
line([VOCALSTRUCTall.global.pretime VOCALSTRUCTall.global.pretime], ylim, 'Color', 'k');
axis tight

