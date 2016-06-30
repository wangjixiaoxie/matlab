function lt_seq_dep_pitch_ListOutliers(Params);
%% LT 11/5/15 - List outliers

AllSyls=Params.PlotLearning.SylFieldsAll;
for i=1:length(AllSyls);
    
    syl=AllSyls{i};
    Ndays=size(Params.SeqFilter.OutlierInds,2);
    % to list outliers detected
    for ii=1:Ndays
        try
            if isfield(Params, 'RecalculateFF');
                % recalcluated windows, use this instead
                OutInds=Params.RecalculateFF.OutlierInds{ii}.(syl);
            else
                
                OutInds=Params.SeqFilter.OutlierInds{ii}.(syl);
                
            end
            
            
            if ~isempty(OutInds);
                disp([num2str(length(OutInds)) ' outliers for ' syl ' on day ' num2str(ii)]);
            end
        catch err
        end
    end
    
end
