function lt_neural_v2_ANALY_VocModel_plot(VOCALSTRUCTall, birdnum, neuronnum, timebin, plotSmoothed)
% for a single bird, neuron, timebin, plots across all syls and presyls distributions of FR.
% plotSmoothed = 1; % then plots smoothed tiemocurse (will still plot binned tcourse_)

%%
binedges = VOCALSTRUCTall.global.binedges;
Pretime = VOCALSTRUCTall.global.pretime;
Posttime = VOCALSTRUCTall.global.posttime;
mindattoplot = VOCALSTRUCTall.global.mindattoplot;
Binsize = binedges(2)-binedges(1);
%%     
z=birdnum; % pick out one to plot
    zz=neuronnum;
%     timebin = 27;
        
        VOCALSTRUCT = VOCALSTRUCTall.birds(z).neurons(zz).dat;
        
        
        lt_figure; hold on;
        plotMean = 0;
        
        % ==== 1) for single time bin, single neuron, plot data
        lt_subplot(2,1,1); hold on;
        tmp1 = binedges(timebin);
        tmp2 = binedges(timebin+1);
        title(['ignoring context; timebin: ' num2str(tmp1) ' to ' num2str(tmp2) ' (onset=' num2str(Pretime) ') sec']);
        ylabel('spks/bin');
        
        % x axis for syl, vary color for sequence (1-back)
        syltypes = unique({VOCALSTRUCT.data_vocalrend.syl});
        for i=1:length(syltypes)
            syl = syltypes{i};
            inds = strcmp({VOCALSTRUCT.data_vocalrend.syl}, syl);
            
            if sum(inds)<mindattoplot
                continue % not enough samples
            end
            
            tmp = {VOCALSTRUCT.data_vocalrend(inds).Spk_bincounts};
            SpkCntMat = cell2mat(tmp');
            
            % -- extract each trial dat
            Y = SpkCntMat(:, timebin);
            X = i;
            
            if plotMean ==1
                lt_plot(X, mean(Y), {'Errors', lt_sem(Y), 'Color', 'k'});
            else
                distributionPlot(Y, 'xValues', X, 'addSpread', 0, 'showMM', 0);
            end
            %     plot(X, Y, 'o');
        end
        lt_plot_zeroline;
        
        
        
        % ==== 2) for single time bin, single neuron, plot data [including sequence
        % info]
        lt_subplot(2,1,2); hold on;
        title('colors = context');
        xlabel('syl');
        ylabel('spks/bin');
        grid on
        % x axis for syl, vary color for sequence (1-back)
        syltypes = unique({VOCALSTRUCT.data_vocalrend.syl});
        plotcols = lt_make_plot_colors(length(syltypes), 0, 0);
        
        for i=1:length(syltypes)
            syl = syltypes{i};
            %     inds = strcmp({VOCALSTRUCT.data_vocalrend.syl}, syl);
            
            % --- break out into preceding syl
            %         presyltypes = unique({VOCALSTRUCT.data_vocalrend.presyl});
            for ii=1:length(syltypes)
                sylpre = syltypes{ii};
                
                inds = strcmp({VOCALSTRUCT.data_vocalrend.syl}, syl) & ...
                    strcmp({VOCALSTRUCT.data_vocalrend.presyl}, sylpre);
                

                if sum(inds)<mindattoplot
%                     disp('NO')
                    continue % not enough samples
                else
%                     disp('YES')
                end
                
                tmp = {VOCALSTRUCT.data_vocalrend(inds).Spk_bincounts};
                SpkCntMat = cell2mat(tmp');
                
                % -- extract each trial dat
                Y = SpkCntMat(:, timebin);
                X = i + 0.3-0.6*rand;
                
                if plotMean==1
                    lt_plot(X, mean(Y), {'Errors', lt_sem(Y), 'Color', plotcols{ii}});
                    lt_plot_text(X+0.01, mean(Y), sylpre, plotcols{ii});
                else
                    distributionPlot(Y, 'xValues', X, 'addSpread', 0, 'color', plotcols{ii}, ...
                        'distWidth', 0.3, 'showMM', 0);
                    %     distributionPlot(Y, 'xValues', X, 'addSpread', 0, 'color', plotcols{ii});
                end
                %     plot(X, Y, 'o');
            end
        end
        lt_plot_zeroline;
        set(gca, 'XTick', 1:length(syltypes), 'XTickLabel', syltypes);
        
        %% DIAGNOSTIC (WAS COMPARING IT TO FORMER CODE, M AKING SURE WAS SIMILAR
        % =============== PLOT RUNNING TIMECOURSE [no smoothuing]
        figcount=1;
        subplotrows=8;
        subplotcols=3;
        fignums_alreadyused=[];
        hfigs=[];
        
        % x axis for syl, vary color for sequence (1-back)
        syltypes = unique({VOCALSTRUCT.data_vocalrend.syl});
        plotcols = lt_make_plot_colors(length(syltypes), 0, 0);
        
        for i=1:length(syltypes)
            syl = syltypes{i};
            %     inds = strcmp({VOCALSTRUCT.data_vocalrend.syl}, syl);
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title(['to ' syl]);
            
            % --- break out into preceding syl
            presyltypes = unique({VOCALSTRUCT.data_vocalrend.presyl});
            for ii=1:length(syltypes)
                sylpre = syltypes{ii};
                
                inds = strcmp({VOCALSTRUCT.data_vocalrend.syl}, syl) & ...
                    strcmp({VOCALSTRUCT.data_vocalrend.presyl}, sylpre);
                
                if sum(inds)<mindattoplot
                    continue % not enough samples
                end
                
                if (0)
                    SpkCntMat = reshape([VOCALSTRUCT.data_vocalrend(inds).Spk_bincounts], sum(inds), ...
                        length(VOCALSTRUCTall.global.binedges)); % collect all into matrix of spk counts (N x bins)
                else
                    tmp = {VOCALSTRUCT.data_vocalrend(inds).Spk_bincounts};
                    SpkCntMat = cell2mat(tmp');
                end
                
                % -- extract each trial dat
                Y = mean(SpkCntMat, 1);
                X = VOCALSTRUCTall.global.binedges(1:end)+Binsize/2;
                
                plot(X, Y./Binsize, 'Color', plotcols{ii});
                
                line([VOCALSTRUCTall.global.pretime VOCALSTRUCTall.global.pretime], ylim ,'Color', 'k')
            end
        end
        lt_plot_zeroline;
        set(gca, 'XTick', 1:length(syltypes), 'XTickLabel', syltypes);
        
        
        % =============== PLOT RUNNING TIMECOURSE [smoothing]
        if plotSmoothed==1;
            figcount=1;
            subplotrows=8;
            subplotcols=3;
            fignums_alreadyused=[];
            hfigs=[];
            
            % x axis for syl, vary color for sequence (1-back)
            syltypes = unique({VOCALSTRUCT.data_vocalrend.syl});
            plotcols = lt_make_plot_colors(length(syltypes), 0, 0);
            
            for i=1:length(syltypes)
                syl = syltypes{i};
                %     inds = strcmp({VOCALSTRUCT.data_vocalrend.syl}, syl);
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                title(['to ' syl]);
                
                % --- break out into preceding syl
                %             presyltypes = unique({VOCALSTRUCT.data_vocalrend.presyl});
                for ii=1:length(syltypes)
                    sylpre = syltypes{ii};
                    
                    inds = strcmp({VOCALSTRUCT.data_vocalrend.syl}, syl) & ...
                        strcmp({VOCALSTRUCT.data_vocalrend.presyl}, sylpre);
                    
                    if sum(inds)<mindattoplot
                        continue % not enough samples
                    end
                    
                    Yspks = {VOCALSTRUCT.data_vocalrend(inds).Spktimes_relonset};
                    
                    window_fr = 0.04;
                    windshift_fr = 0.002;
                    [xbin, ~, ~, ymean_hz, ysem_hz, ~, ystd_hz] = ...
                        lt_neural_plotRastMean(Yspks, window_fr, windshift_fr, 0, '');
                    
                    shadedErrorBar(xbin, ymean_hz, ysem_hz, {'Color', plotcols{ii}}, 1);
                    
                    line([VOCALSTRUCTall.global.pretime VOCALSTRUCTall.global.pretime], ylim ,'Color', 'k')
                end
            end
            lt_plot_zeroline;
            set(gca, 'XTick', 1:length(syltypes), 'XTickLabel', syltypes);
    end