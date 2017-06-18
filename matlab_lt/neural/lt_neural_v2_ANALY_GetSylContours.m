%% lt 6/16/17 - given SegmentsExtract, output vector (0s and 1s) for syl/gaps
function SylContours = lt_neural_v2_ANALY_GetSylContours(segextract, numsyltrialstoplot, maxdur)

% numsyltrialstoplot = random subset of trials
% maxdur = [] , in sec, clips data.

% OTUPUT
% SylContours = []; % trial x timebin

% NOTE: !!! can optimize to make faster
% CURRENTLY ONLY OUTPUTS MATRIX (so forced to have same num time bins)
% time bins is 1ms

%%
TrialInds = randi(length(segextract), 1, numsyltrialstoplot); % which random trials to get

SylContours = [];
for kk=1:length(TrialInds)
    trialindtmp = TrialInds(kk);
    trialdur = max(segextract(trialindtmp).FRsmooth_xbin);
    tmpons = ones(1, ceil((trialdur)*1000)); % one bin per ms, for onsets
    
    % fill in all gaps (after first offset)
    for k=1:length(segextract(trialindtmp).sylOffTimes_RelDataOnset)
        tmpind = find(segextract(trialindtmp).sylOnTimes_RelDataOnset>...
            segextract(trialindtmp).sylOffTimes_RelDataOnset(k), 1, 'first');
        
        if ~isempty(tmpind)
            % fil in this gap
            ind1 = ceil(1000*segextract(trialindtmp).sylOffTimes_RelDataOnset(k));
            ind2 = ceil(1000*segextract(trialindtmp).sylOnTimes_RelDataOnset(tmpind));
            tmpons(ind1:ind2) = 0;
        else
            % end of song, 
            tmpons(ceil(1000*segextract(trialindtmp).sylOffTimes_RelDataOnset(k)):end) ...
                = 0;
        end
    end
    
    % fill in if trial started on a gap
    if segextract(trialindtmp).sylOnTimes_RelDataOnset(1)<segextract(trialindtmp).sylOffTimes_RelDataOnset(1);
        % then starts during off
        tmpons(1:ceil(segextract(trialindtmp).sylOnTimes_RelDataOnset(1)*1000))=0;
    end
    
    % -- cut to max dur to fit into matrix
    tmpons = tmpons(1:ceil(1000*(maxdur)));
    
    % =========== OUTPUT
    SylContours = [SylContours; tmpons];
 %%   
    if (0) % troubleshooting (overlaying onset offset plotted as balls with outputted contours
        lt_figure; hold on;
        plot(segextract(trialindtmp).sylOnTimes_RelDataOnset, 1, 'og');
        plot(segextract(trialindtmp).sylOffTimes_RelDataOnset, 1, 'or');
        x = (1:length(tmpons))/1000;
        plot(x, tmpons, 'Color', plotcols{nn});
    end
end
