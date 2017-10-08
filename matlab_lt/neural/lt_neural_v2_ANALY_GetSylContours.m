%% lt 6/16/17 - given SegmentsExtract, output vector (0s and 1s) for syl/gaps
function [SylContours, x] = lt_neural_v2_ANALY_GetSylContours(segextract, numsyltrialstoplot, maxdur)

% numsyltrialstoplot = random subset of trials
% maxdur = [] , in sec, clips data.

% OTUPUT
% SylContours = []; % trial x timebin

% NOTE: !!! can optimize to make faster
% CURRENTLY ONLY OUTPUTS MATRIX (so forced to have same num time bins)
% time bins is 1ms

%% shorten desired trial number if not enough
numsyltrialstoplot = min([length(segextract) numsyltrialstoplot]);

%%
TrialInds = randi(length(segextract), 1, numsyltrialstoplot); % which random trials to get

SylContours = nan(numsyltrialstoplot, ceil(1000*(maxdur))); % fill with nans, if actual contyours too short, then edges will have nans


% --- figure out min trial dur
trialdur = min([segextract.global_offtime_motifInclFlank] - ...
    [segextract.global_ontime_motifInclFlank]) - 0.002;
% - old, too slow
if (0)
    if ~isfield(segextract, 'FRsmooth_xbin')
        segextract = lt_neural_SmoothFR(segextract, []);
    end
    trialdur = max(segextract(1).FRsmooth_xbin_CommonTrialDur);
end

maxdur = min([maxdur trialdur]);

for kk=1:length(TrialInds)
    trialindtmp = TrialInds(kk);
    
    % --- figure out trial dur (use all clust, since just want trial dur)
    
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
    %     SylContours = [SylContours; tmpons];
    SylContours(kk, 1:length(tmpons)) = tmpons;
    x = (1:length(tmpons))/1000;
    %%
    if (0) % troubleshooting (overlaying onset offset plotted as balls with outputted contours
        lt_figure; hold on;
        plot(segextract(trialindtmp).sylOnTimes_RelDataOnset, 1, 'og');
        plot(segextract(trialindtmp).sylOffTimes_RelDataOnset, 1, 'or');
        x = (1:length(tmpons))/1000;
        plot(x, tmpons, 'Color', plotcols{nn});
    end
end
