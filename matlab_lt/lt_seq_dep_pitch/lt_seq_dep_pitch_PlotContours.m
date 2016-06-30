function lt_seq_dep_pitch_PlotContours(AllDays_RawDatStruct, Params, day, checkMUSC);
% checkMUSC = 1; then checks musc data, not PBS - NOTE NOT ACTUAL TIME
% WINDOWS USED FOR MUSC (OR PBS)


%%

if ~exist('checkMUSC','var');
    checkMUSC = 0;
end


%% LT 9/7/15 - quick, to plot all syls for a given day, overlaid with time windows used to calculate pitch


if checkMUSC==1;
    datafield='data_MUSC';
else
    datafield='data_WithOutlier';
end

count=1;
SubplotsPerFig=9;
subplotrows=3;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];

% ====
syls_list=fieldnames(AllDays_RawDatStruct{day}.(datafield));

for i=1:length(syls_list);
    syl=syls_list{i};
    
    if isfield(AllDays_RawDatStruct{day}.(datafield), syl);
        
        [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
        title(syl);
        
        PCmat=cell2mat(AllDays_RawDatStruct{day}.(datafield).(syl)(:,2));
        
        plot(PCmat');
        
        N=size(PCmat,1);
        
        Ylim=ylim;
        
        
        axis tight
        
        % ============================ put lines for window
        if checkMUSC==0;
            ind=i;
        else
            % becuase order of stuff in list of time windows is only
            % correct for pbs
            ind=find(strcmp(fieldnames(AllDays_RawDatStruct{day}.data_WithOutlier), syl));
        end
        
        if isfield(Params, 'RecalculateFF'); % then means changed window, e.g. for WN songs - all days will use same
            twindow=Params.RecalculateFF.pc_time_window_list(:, ind);
            
            line([twindow(1) twindow(1)], ylim, 'Color','r', 'LineWidth',2);
            line([twindow(2) twindow(2)], ylim, 'Color','r', 'LineWidth', 2);
            
        end
        
        twindow=Params.SeqFilter.pc_time_window_list{day}(:,ind);
        
        
        line([twindow(1) twindow(1)], ylim, 'Color','k', 'LineWidth', 2);
        line([twindow(2) twindow(2)], ylim, 'Color','k', 'LineWidth', 2);
        
        % == annotate number of outliers removed (out of total)
        if checkMUSC==0;
            if isfield(Params, 'RecalculateFF'); % then means changed window, e.g. for WN songs - all days will use same
                n_outliers   =         numel(Params.RecalculateFF.OutlierInds{day}.(syl));
            else
                n_outliers=numel(Params.SeqFilter.OutlierInds{day}.(syl));
                
            end
            
            
            Ylim=ylim;
            lt_plot_text(5, Ylim(2)-(Ylim(2)-Ylim(1))/6, ['N = ' num2str(N) ' (' num2str(n_outliers) ' o.l.)'], 'k')
        end;
    end
end