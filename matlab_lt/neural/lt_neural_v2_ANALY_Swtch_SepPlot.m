function lt_neural_v2_ANALY_Swtch_SepPlot(Datstruct, MOTIFSTATS_Compiled, ...
    SwitchStruct, getdprime, expttypewanted)
%% lt 9/17/17 - plots 
plotanova =1;

% expttypewanted = 'one targ context';
% % expttypewanted = 'mult targ context - samedir';
% % expttypewanted = 'mult targ context - diff dir';
% % expttypewanted=''; % get all

%% ======================================== PLOT FOR EACH TYPE OF EXPERIMENT

numexpts = length(Datstruct.exptcount);


AllExptcount = [];
AllNeuronnum = [];
AllPairtype = [];
AllDprime_base = {};
AllDprime_train = {};
AllDprime_change = {};
AllXbin_dprime = {};

AllAnova_change = {};
AllAnova_xbin = {};

AllBirdExptSw = {};
AllSylPairs = {};

for ee = 1:numexpts
    if ~isempty(expttypewanted)
        if ~strcmp(Datstruct.exptcount(ee).expttype, expttypewanted)
            continue
        end
    end
    
    birdname = SwitchStruct.bird(Datstruct.exptcount(ee).birdnum_old).birdname;
    exptname = SwitchStruct.bird(Datstruct.exptcount(ee).birdnum_old).exptnum(Datstruct.exptcount(ee).exptnum_old).exptname;
    swnum = Datstruct.exptcount(ee).switchnum_old;
    
    
    numneurons = length(Datstruct.exptcount(ee).neuronnum);
    
    for nn=1:numneurons
        
        nummotifpairs = length(Datstruct.exptcount(ee).neuronnum(nn).motifpair);
        
        %         lt_figure; hold on;
        %         subplot(221); hold on;
        %         title('');
        %         subplot(222); hold on;
        %         title('');
        %         subplot(223); hold on;
        %         title('');
        %         subplot(224); hold on;
        %         title('');
        %
        %
        targdirs = Datstruct.exptcount(ee).learndirs;
        disp(['---------------- ' targdirs]);
        
        for mm = 1:nummotifpairs
            
            
            issamesyl = Datstruct.exptcount(ee).neuronnum(nn).motifpair(mm).issamesyl;
            istargsyl = Datstruct.exptcount(ee).neuronnum(nn).motifpair(mm).istargsyl;
            istargcontext = Datstruct.exptcount(ee).neuronnum(nn).motifpair(mm).istargcontext;
            learndir = Datstruct.exptcount(ee).neuronnum(nn).motifpair(mm).learndir;
            
            
            % ============================ plot all kinds of pairs
            
            if getdprime==1
            dprime_base = Datstruct.exptcount(ee).neuronnum(nn).motifpair(mm).DAT_dprime_base;
            dprime_train = Datstruct.exptcount(ee).neuronnum(nn).motifpair(mm).DAT_dprime_train;
            xbin = Datstruct.exptcount(ee).neuronnum(nn).motifpair(mm).DAT_xbins;
            dprime_change = dprime_train - dprime_base;
            end
            
            anova_change = Datstruct.exptcount(ee).neuronnum(nn).motifpair(mm).DAT_omega_change;
            xbin_anova = Datstruct.exptcount(ee).neuronnum(nn).motifpair(mm).DAT_anovaxbins;
            
            sylpairs = Datstruct.exptcount(ee).neuronnum(nn).motifpair(mm).motifs;
            
            % 1) diff syls
            if issamesyl == 0
                %                 subplot(221);
                %                 plot(xbin, dprime_change, '-r');
                
                AllPairtype = [AllPairtype 1];
                disp(['diff syl: ' sylpairs]);
            else
                assert(length(unique(istargsyl)) ==1, 'asdf');
                
                % 2) same syl (target syl), both contexts nontarget
                if all(istargsyl==1) & all(istargcontext==0)
                    %                     subplot(222);
                    %                     plot(xbin, dprime_change, '-r');
                    
                    AllPairtype = [AllPairtype 2];
                    
                    disp(['same syl (targ syl, not targ): ' sylpairs]);
                    
                    
                    % 3) same syl (target syl), one target context
                elseif all(istargsyl==1) & sum(istargcontext)==1
                    %                     subplot(223);
                    %
                    %                     plot(xbin, dprime_change, '-r');
                    
                    AllPairtype = [AllPairtype 3];
                    
                    disp(['same syl (targ syl, one context targeted): ' sylpairs]);
                    
                    
                    % 4) same syl (not targ syl)
                elseif all(istargsyl ==0)
                    assert(all(istargcontext ==0));
                    
                    %                     subplot(224);
                    %
                    %                     plot(xbin, dprime_change, '-r');
                    %
                    AllPairtype = [AllPairtype 4];
                    disp(['same syl (not targ syl)' sylpairs]);
                    
                    % 5) both targ context (same dir)
                elseif all(istargcontext==1) & length(unique(learndir))==1
                    AllPairtype = [AllPairtype 5];
                    
                    disp(['same syl (both targ context, same dir): ' sylpairs]);
                    
                    
                    
                    % 6) both targ context (opposite dir)
                elseif all(istargcontext==1) & length(unique(learndir))>1
                    AllPairtype = [AllPairtype 6];
                    disp(['same syl (both targ context, opposite dir): ' sylpairs]);
                    
                    
                end
                
            end
            
            % =================== OUTPUT
            AllExptcount = [AllExptcount ee];
            AllNeuronnum = [AllNeuronnum nn];
            
            AllAnova_change = [AllAnova_change anova_change];
            AllAnova_xbin = [AllAnova_xbin xbin_anova];
            
            if getdprime==1
            AllDprime_base = [AllDprime_base dprime_base];
            AllDprime_train = [AllDprime_train dprime_train];
            AllDprime_change = [AllDprime_change dprime_change];
            AllXbin_dprime = [AllXbin_dprime xbin];
            end
            
            AllBirdExptSw = [AllBirdExptSw [birdname '-' exptname '-sw' num2str(swnum)]];
            AllSylPairs = [AllSylPairs [sylpairs{1} '-' sylpairs{2}]];
        end
        if (0)
            % if want to visually inspect the syls that make up each pair.
            pause
        end
    end
    
end

% ==== plot types
PlotTypeList = {'diff syls', 'same (targ) syl, no targ C', 'same (targ) syl, one targ C', ...
    'same (not targ) syl', 'same syl, 2 targ C (samedir)', 'same syl, 2 targ C (bidir)'};


%% ============================== PLOT [each expt]
premotorWind = MOTIFSTATS_Compiled.birds(1).exptnum(1).MOTIFSTATS.params.motif_predur;
NumExpts = max(AllExptcount);
for ee =1:NumExpts
    lt_figure;
    hsplots = [];
    NumPairTypes = max(AllPairtype);
    
        plotname = unique(AllBirdExptSw(AllExptcount==ee));
        assert(length(plotname)<2, 'asfds');
        
    for pp=1:NumPairTypes
        hsplot = lt_subplot(3,2,pp); hold on;
        hsplots = [hsplots hsplot];
        inds = AllPairtype==pp & AllExptcount==ee;
        
        
        
        % --- runs
        if plotanova==1
        NeurSimChangeCell = AllAnova_change(inds);
        XbinCell = AllAnova_xbin(inds);
        else
        NeurSimChangeCell = AllDprime_change(inds);
                XbinCell = AllXbin_dprime(inds);
        end
        
        func = @(X) length(X);
        tmp = min(cellfun(func, NeurSimChangeCell));
        Ymat = nan(length(NeurSimChangeCell), tmp); % trial x bins
        for i=1:length(NeurSimChangeCell)
            Ymat(i,:) = NeurSimChangeCell{i}(1:tmp);
        end
        
       % -- plot
        if ~isempty(Ymat)
        x = XbinCell{1}(1:tmp);
            % raw
            plot(x, Ymat', '-', 'Color', [0.6 0.6 0.6]);
        end
        if size(Ymat,1)>1
            % mean
            ymean = mean(Ymat,1);
            ysem = lt_sem(Ymat);
            shadedErrorBar(x, ymean, ysem, {'Color', 'r'});
            ylim([-1 1]);
        end
        line([premotorWind premotorWind], ylim);
        lt_plot_zeroline;
        title(PlotTypeList{pp});
        
    end
    if ~isempty(plotname)
    lt_subtitle(['expt' num2str(ee) ' (' plotname{1} ')']);
    else
        lt_subtitle(['expt' num2str(ee)]);
    end
    linkaxes(hsplots, 'xy');
    
end


%% PLOTS



%% ===================================== PLOT [AVERAGE ACROSS EXPTS]
lt_figure;
hsplots = [];
NumPairTypes = max(AllPairtype);
for pp=1:NumPairTypes
    hsplot = lt_subplot(3,2,pp); hold on;
    hsplots = [hsplots hsplot];
    inds = AllPairtype==pp;
    
    
    % --- runs
        if plotanova==1
        NeurSimChangeCell = AllAnova_change(inds);
        XbinCell = AllAnova_xbin(inds);
        else
        NeurSimChangeCell = AllDprime_change(inds);
                XbinCell = AllXbin_dprime(inds);
        end
    
    func = @(X) length(X);
    tmp = min(cellfun(func, NeurSimChangeCell));
    Ymat = nan(length(NeurSimChangeCell), tmp); % trial x bins
    for i=1:length(NeurSimChangeCell)
        Ymat(i,:) = NeurSimChangeCell{i}(1:tmp);
    end
    
    % -- plot
    if ~isempty(Ymat)
    x = XbinCell{1}(1:tmp);
        % raw
        plot(x, Ymat', '-', 'Color', [0.6 0.6 0.6]);
    end
    if size(Ymat,1)>1
        % mean
        ymean = mean(Ymat,1);
        ysem = lt_sem(Ymat);
        shadedErrorBar(x, ymean, ysem, {'Color', 'r'});
        ylim([-1 1]);
    end
    line([premotorWind premotorWind], ylim);
    lt_plot_zeroline;
            title(PlotTypeList{pp});

end
lt_subtitle(['all expt']);
linkaxes(hsplots, 'xy');






