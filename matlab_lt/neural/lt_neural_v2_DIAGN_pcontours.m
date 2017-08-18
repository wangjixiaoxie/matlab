function lt_neural_v2_DIAGN_pcontours(SummaryStruct, useDiffColors)



%%  LT 8/9/17 - plots all pitch contours and extracted window that was used
% for learning experiments, also plots PC over time.



numbirds = length(SummaryStruct.birds);

for i=1:numbirds
    numneurons = length(SummaryStruct.birds(i).neurons);
    
    figcount=1;
    subplotrows=7;
    subplotcols=4;
    fignums_alreadyused=[];
    hfigs=[];
    
    birdname = SummaryStruct.birds(i).birdname;
    
    for ii=1:numneurons
        
        % -- go to data dir
        cd(SummaryStruct.birds(i).neurons(ii).dirname);
        
        
        
        % -- load all PC data
        load('extractFF.mat'); % FFvals
        load('extractFF_params.mat'); % Params
        load('extractPC.mat'); % PCvals
        
        cd ..
        batchf = SummaryStruct.birds(i).neurons(ii).batchfilename;
        chan = SummaryStruct.birds(i).neurons(ii).channel;
        [SongDat] = lt_neural_ExtractDat(batchf, chan, 0);
        
        
        tvals = Params.tfirstbin:Params.tbin:0.3; % make this long, then shorten
        
        % --- for every syllable (ignore context) overlay all PCs + lines
        % for windows
        sylList = unique(SongDat.AllLabels);
        for j=1:length(sylList)
            syl = sylList(j);
            
            if strcmp(syl, '-')
                continue
            end
            
            if ~any(strcmp(Params.cell_of_FFtimebins, syl))
                continue
                % since this syl FF not gotten
            end
            
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title([birdname '-n' num2str(ii) '-' syl]);
            
            
            % - confirm that all labels and FFvals are matched
            assert(length(FFvals) == length(SongDat.AllLabels), 'asdfads')
            
            % -- get inds for this syl
            inds = strfind(SongDat.AllLabels, syl);
            
            % ----- 1) Plot PCs
            pcvals_toplot = PCvals(inds);
            assert(~(any(cellfun('isempty', pcvals_toplot) ==1) & any(cellfun('isempty', pcvals_toplot) ==0)), ...
                'problem: should be either all empty or all full');
            
            
            for jj=1:length(pcvals_toplot)
                
                if useDiffColors==1
                    plotcol = [rand rand rand];
                    plotstyle = '-';
                else
                    plotcol = [0.5 0.5 0.5];
                    plotsyle = ':';
                end
                
                plot(tvals(1:length(pcvals_toplot{jj})), pcvals_toplot{jj}, plotstyle, 'Color', plotcol);
                
            end
            
            axis tight;
            
            % ------- 2) Plot lines for windows
            timewindow = Params.cell_of_FFtimebins{find(strcmp(Params.cell_of_FFtimebins, syl))+1};
            line([timewindow(1) timewindow(1)], ylim, 'Color', 'k');
            line([timewindow(2) timewindow(2)], ylim, 'Color', 'k');
            
            ffwindow = Params.cell_of_freqwinds{find(strcmp(Params.cell_of_freqwinds, syl))+1};
            line(xlim, [ffwindow(1) ffwindow(1)], 'Color', 'k');
            line(xlim, [ffwindow(2) ffwindow(2)], 'Color', 'k');
            
            line([timewindow(1) timewindow(2)], [mean(FFvals(inds)) mean(FFvals(inds))], 'Color','k');
            
            % --- zoom into region of interest
            flankdur = 0.018;
            xlim([timewindow(1)-flankdur, timewindow(2)+flankdur]);
            
        end
        
        
        
        % --- automatically get distribution of WN onsets and make sure that
        % is outside of FF calculation window
        
        
        
        % --- plot FF distribution, look for clear outliers, which indicate
        % errors when getting FF
        
        
        
        
        SummaryStruct.birds(i).neurons(ii).POSTINFO;
        
    end
end
