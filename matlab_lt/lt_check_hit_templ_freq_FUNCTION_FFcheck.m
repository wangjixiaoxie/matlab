threshold = [10 30 50 70 90]; %Percentiles for new threshold

%% Calculates FF for all labeled syllables
check_stuff.(sofinterest).freq.vals = evtaf_freq2_LT(parameters.batchfile, parameters.freq_range, parameters.sofinterest,...
    parameters.sofinterest_pre,parameters.sofinterest_post, 128, 'obs0', 1,evtaf_ver);


%% Calculates and plot various statistics
if size(check_stuff.(sofinterest).freq.vals,1)>0; % i.e. make sure got data
    % == Calculates running FF for all labeled syllables
    try
        running_avg_window = 10;
        check_stuff.(sofinterest).freq.running_avg = db_runningaverage(check_stuff.(sofinterest).freq.vals(:,2),running_avg_window);
        lt_figure; hold on;
        grid on
        title(['Running avg of FF (window: ' num2str(running_avg_window) ')'])
        ylabel('FF')
        xlabel('Syllable number')

        plot(check_stuff.(sofinterest).freq.running_avg,'-','LineWidth',2,'MarkerSize',4);
        lt_plot(check_stuff.(sofinterest).freq.running_avg);
        
        saveas(figure(gcf),['check_template_' parameters.today_date '/Running_avg_' parameters.sofinterest{i} '_' parameters.today_date], 'fig')
    catch err
    end
    
    % == calculates cumulative FF for all labeled syllables
%     for jj = 1:length(check_stuff.(sofinterest).freq.vals(:,2))
%         check_stuff.(sofinterest).freq.cumulative_avg(jj) = mean(check_stuff.(sofinterest).freq.vals(1:jj,2));
%     end
%     figure
%     grid on
%     title('Cumulative avg of FF')
%     ylabel('FF')
%     xlabel('Syllable number')
%     hold on
%     plot(check_stuff.(sofinterest).freq.cumulative_avg,'-','LineWidth',2,'MarkerSize',4);
%     saveas(figure(gcf),[SaveDir3 '/Cumulative_avg_' parameters.sofinterest '_' parameters.today_date], 'fig')
    
    %== Summary statistics of FF for labeled syllables
    check_stuff.(sofinterest).freq.mean = mean(check_stuff.(sofinterest).freq.vals(:,2));
    check_stuff.(sofinterest).freq.sd = std(check_stuff.(sofinterest).freq.vals(:,2));
    check_stuff.(sofinterest).freq.median = median(check_stuff.(sofinterest).freq.vals(:,2));
    check_stuff.(sofinterest).freq.iqr = iqr(check_stuff.(sofinterest).freq.vals(:,2));
    
    check_stuff.(sofinterest).freq.prctiles = num2str(threshold);
    check_stuff.(sofinterest).freq.threshold = prctile(check_stuff.(sofinterest).freq.vals(:,2), threshold);
    
    
    %To see if FF  is normally distributed - LT commented out, not needed
    %     try
    %         if kstest(syllables.(check_pitch.parameters.sofinterest{i}).freq.vals(:,2)) == 1;
    %             syllables.(check_pitch.parameters.sofinterest{i}).freq.vals_normal = 'no';
    %         elseif kstest(syllables.(check_pitch.parameters.sofinterest{i}).freq.vals(:,2)) == 0;
    %             syllables.(check_pitch.parameters.sofinterest{i}).freq.vals_normal = 'yes';
    %         end
    %     catch err
    %         continue
    %     end
    
    % == Histogram of FF distrobution
    lt_figure
    [f,x] = hist(check_stuff.(sofinterest).freq.vals(:,2),30);
    bar(x,f/trapz(x,f))
    title(['FF Probability density function for ' parameters.sofinterest])
    xlabel('Frequency (Hz)')
    ylabel('Density')
    saveas(figure(gcf), [SaveDir3 '/FF_PDF_' parameters.sofinterest '_' parameters.today_date], 'fig')
    
    % == Plots FF of labeled syllables
    lt_figure(), hold on; grid on
    lt_plot(24.*(check_stuff.(sofinterest).freq.vals(:,1)-floor(check_stuff.(sofinterest).freq.vals(:,1))),...
        check_stuff.(sofinterest).freq.vals(:,2));
    title(['FF of ' parameters.sofinterest])
    xlabel('Time (hours)')
    ylabel('Frequency (Hz)')
    saveas(figure(gcf),[SaveDir3 '/check_pitch_shift_' parameters.sofinterest '_' parameters.today_date],'fig')
    
    % Plot running percentiles (30th and 70th) to determine new
    % evtaf threshold
    total_rends=length(check_stuff.(sofinterest).freq.vals(:,2));
    for kk=1:total_rends; % will give percentiles for all possible size rendition windows going back in time
        rends_to_take=total_rends-kk+1:total_rends;
        percentiles_AllWindSizes(kk,:)=prctile(check_stuff.(sofinterest).freq.vals(rends_to_take,2),[5 30 70 95]);
    end
    lt_figure; hold on;
    plot(flipud(percentiles_AllWindSizes(:,2)),'og');
    plot(flipud(percentiles_AllWindSizes(:,3)),'ob');
    xlabel('how far back the beginning of the sample window goes (e.g. most right is closest to the most recent data)');
    ylabel('Pitch (hz)');
    title('30th and 70th percentiles of pitch: Each point is a recalculation, using different windows');
    saveas(figure(gcf),[SaveDir3 '/moving_threshold_percentiles_' parameters.sofinterest '_' parameters.today_date],'fig')
    
    
    if get_WN_hits==1 && get_FF==1;
        % PLOT - FF distribution (offline), distinguishing hits and
        % escapes (actual).
        % Note: assumes that offline match and online hits use the
        % same template and cntrng parameters, and thus the
        % timeline of hits and escapes (as calculated based on
        % labels and
        HitOrEscape=[];
        MatchLogicals=[];
        for kk=1:size(check_stuff.(sofinterest).hit.trigger,2);
            HitOrEscape=[HitOrEscape; check_stuff.(sofinterest).hit.trigger(kk).Labeled_HitsEscapes];
            MatchLogicals=[MatchLogicals; check_stuff.(sofinterest).match.trigger(kk).Labeled_HitsEscapes];
        end
        % Matchlogical contains all labeled things that are matched offline.
        %                 HitOrEscape counts all labeled things that were hit online. However, evtaf_freq calcualtes
        %                 frequency only for those things that are labeled, and matched offline. So, to ask what frequency
        %                 values correspond to hit vs. miss trials, I need to take the HitorEscape logicals that are
        %                 for the matched trials only:
        HitOrEscape_MatchTrials=HitOrEscape(logical(MatchLogicals));
        
        try
            HitFreqVals=check_stuff.(sofinterest).freq.vals(logical(HitOrEscape_MatchTrials),:);
            
            lt_figure; hold on;
            [f,x] = hist(check_stuff.(sofinterest).freq.vals(:,2),30);
            bar(x,f/trapz(x,f),'g');
            f = hist(HitFreqVals(:,2),x);
            bar(x,f/trapz(x,f),'r');
            title('PDF of FF (all labeled and matched offline): Actual hits (red) and escapes (blue)')
            xlabel('Frequency (Hz)')
            ylabel('Density')
            saveas(figure(gcf), [SaveDir3 '/Hit_Escape_FF_PDF_' parameters.sofinterest '_' parameters.today_date], 'fig')
            
            check_stuff.HitFreqVals=HitFreqVals;
            try
                check_stuff.EscapeFreqVals=check_stuff.(sofinterest).freq.vals(logical(1-HitOrEscape_MatchTrials),:);
            catch err
                disp('Error, number of matches offline and online are likely not same.  likely becuase online and offlien templates are different. To solve, can perform offline and online analyses separately');
            end
            
            % Plot scatter plot as well
            lt_figure; hold on;
            lt_plot(24.*(check_stuff.(sofinterest).freq.vals(:,1)-floor(check_stuff.(sofinterest).freq.vals(:,1))),...
                check_stuff.(sofinterest).freq.vals(:,2),{'Color','g'});
            lt_plot(24.*(HitFreqVals(:,1)-floor(HitFreqVals(:,1))),...
                HitFreqVals(:,2), {'Color','r'});
            
            title(['FF of ' parameters.sofinterest '; actual hits (red) and escapes (blue)'])
            xlabel('Time (hours)')
            ylabel('Frequency (Hz)')
            saveas(figure(gcf),[SaveDir3 '/FF_scatter_Hits_escapes_' parameters.sofinterest '_' parameters.today_date],'fig')
        catch err
            disp('NOTE: Error likely due to following: trigger is slightly outside of boundaries you created based on your segmentation thresholds.  TRIGLABEL_2_LT decides that that hit is hit the closest syllable. i.e. takes the average, i.e. the syl that is closest to trigger')
            disp('Sometimes that is actually the syllable it really hit (see line 209, in triglabel2_LT_evtafv4_v2, where it takes that average). That leads to error. reason follows:');
            disp('However, evtaf_freq considers that a miss. Therefore it does not calculate the frequency for that syl')
            disp('The current function reaches an error because it thinks that that trigger SHOULD (based on triglabel) have a freq calculated for it, when it fact it does not. The array made by evtaf freq is too small');
            disp(' A solution would be to not consider this a hit. but that is undercounting hits. best way is to change code to calculate frequency for all hits, even if they are not at label, then just sort out labels using waht is defined a label by triglabel');
        end
        
    end
end
