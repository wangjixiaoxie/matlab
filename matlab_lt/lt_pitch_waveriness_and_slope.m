%% LT 11/5/13 - looking at pitch waveriness and slope.

% need to first get pitch contour data using
% lt_db_contour_and_FF_analysis_over_time_v4.,
% which uses jc_pitchcontourFV to get the P.C.

% PC should be saved in a structure of:
% multiple_pitch.pc_all.[syllable letter]{day #}(data points, renditions)


%% jsut visualizing the PCs (not clean, this is borrowed from rd3gr35_analysis script
% go to all_days... folder and load [birdname].mat.
% looking at day 3, which is postlesion 2

figure; plot(multiple_pitch.pc_all.t{3});
figure; plot(multiple_pitch.pc_all.t{2});

% do the same, but for each contour subtract the mean (of that single
% contour)

contour_start=220;
contour_end=650;

for day=1:10;
    for i=1:size(multiple_pitch.pc_all.t{day},2);
        mean_pitch_of_contour{day}(i)=mean(multiple_pitch.pc_all.t{day}(contour_start:contour_end,i),1);
        multiple_pitch.pc_all_TrialMeanSubtracted.t{day}(:,i)=multiple_pitch.pc_all.t{day}(:,i)-mean_pitch_of_contour{day}(i);
    end
    figure(day); plot(multiple_pitch.pc_all_TrialMeanSubtracted.t{day}(:,1:50)); title(['P.C. for day' num2str(day)]);
end

% finding STD of deviation from mean (of PC) at each timepoint (across
% renditions)

for day=1:10;
    for ii=1:size(multiple_pitch.pc_all.t{day},1);
        STDAcrossRenditions_OfDeviationFromMean_ofPC{day}(ii)=std(multiple_pitch.pc_all_TrialMeanSubtracted.t{day}(ii,:));
    end
     figure(day+3); plot(STDAcrossRenditions_OfDeviationFromMean_ofPC{day}(180:600)); title(['STD of P.C. across renditions for day' num2str(day)]);   
end


%% quantification

%% FIRST, for each rendition find the how different it is from a "flat" pitch


clear pitch_waveriness

for i = 1:length(multiple_pitch.parameters.syllables);
    pitch_waveriness.deviation_from_flat.(multiple_pitch.parameters.syllables(i)).parameters.window_range=[150 220];
    for j = multiple_pitch.parameters.days{3}:multiple_pitch.parameters.days{4};
        jj=j-multiple_pitch.parameters.days{3}+1;
        try
            pitch_waveriness.deviation_from_flat.(multiple_pitch.parameters.syllables(i)).windowed_PC{jj}=...
                multiple_pitch.pc_all.(multiple_pitch.parameters.syllables(i)){jj}...
                (pitch_waveriness.deviation_from_flat.(multiple_pitch.parameters.syllables(i)).parameters.window_range(1):...
                pitch_waveriness.deviation_from_flat.(multiple_pitch.parameters.syllables(i)).parameters.window_range(2),:); % getting the pc window first
            pitch_waveriness.deviation_from_flat.(multiple_pitch.parameters.syllables(i)).mean_pitch{jj}=...
                mean(pitch_waveriness.deviation_from_flat.(multiple_pitch.parameters.syllables(i)).windowed_PC{jj}); % getting the mean pitch acros window
            pitch_waveriness.deviation_from_flat.(multiple_pitch.parameters.syllables(i)).windowed_PC_minus_mean{jj}=...
                pitch_waveriness.deviation_from_flat.(multiple_pitch.parameters.syllables(i)).windowed_PC{jj}...
                -repmat(pitch_waveriness.deviation_from_flat.(multiple_pitch.parameters.syllables(i)).mean_pitch{jj},...
                size(pitch_waveriness.deviation_from_flat.(multiple_pitch.parameters.syllables(i)).windowed_PC{jj},1),1); % subtract mean for each rendition
            pitch_waveriness.deviation_from_flat.(multiple_pitch.parameters.syllables(i)).area_under_mean_subtracted_PC{jj}=...
                sum(abs(pitch_waveriness.deviation_from_flat.(multiple_pitch.parameters.syllables(i)).windowed_PC_minus_mean{jj}),1)./...
                size(pitch_waveriness.deviation_from_flat.(multiple_pitch.parameters.syllables(i)).windowed_PC{jj},1); % find area under curve  (sum/time range)
            pitch_waveriness.deviation_from_flat.(multiple_pitch.parameters.syllables(i)).day_mean_of_area_under_curve{jj}=...
                mean(pitch_waveriness.deviation_from_flat.(multiple_pitch.parameters.syllables(i)).area_under_mean_subtracted_PC{jj});
        catch err
            continue
        end
    end % getting the pc window first
end

% Plot area under curve of all renditions
%   Calculates baseline
%     fprintf('\nCalculating baseline...')
    for i = 1:length(multiple_pitch.parameters.syllables)
        try
            for j = 1:multiple_pitch.parameters.lesion_begin;
                temp_base.mean{j} = pitch_waveriness.deviation_from_flat.(multiple_pitch.parameters.syllables(i)).area_under_mean_subtracted_PC{j};
            end
            temp_base.mean = temp_base.mean(~cellfun('isempty',temp_base.mean));
            temp_base.mean = cell2mat(temp_base.mean);
            pitch_waveriness.deviation_from_flat.(multiple_pitch.parameters.syllables(i)).baseline_mean = mean(temp_base.mean);
        catch err
            continue
        end
        
        clear('temp_base')
        
    end
%     fprintf('done!\n')
    
    % Makes a figure for each day 
%     fprintf('\nCreating figures of FF over time...')
    for i = 1:length(multiple_pitch.parameters.syllables)
        figure(), hold on
        title([multiple_pitch.parameters.nameofbird '   duration: ' num2str(multiple_pitch.parameters.days{4}-multiple_pitch.parameters.days{3}+1)...
            ' days   syllable: ' multiple_pitch.parameters.syllables(i)])
        ylabel('Frequency (Hz)')
        xlabel('Time (days)')
        if isempty(multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).time_and_FF{end}) == 0
            xlim([-1-multiple_pitch.parameters.lesion_begin...
                ceil(multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).time_and_FF{end}(end,1))+1])
        else
            xlim([-1-multiple_pitch.parameters.lesion_begin ...
                -1-multiple_pitch.parameters.lesion_begin + length(multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).time_and_FF)])
        end
        for j = multiple_pitch.parameters.days{3}:multiple_pitch.parameters.days{4}
            try
                
                %plots FF for time range
                plot(multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).time_and_FF{j-multiple_pitch.parameters.days{3}+1}(:,1),...
                    pitch_waveriness.deviation_from_flat.(multiple_pitch.parameters.syllables(i)).area_under_mean_subtracted_PC{j-multiple_pitch.parameters.days{3}+1},...
                    '.')                
                plot(multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).time_and_FF{j-multiple_pitch.parameters.days{3}+1}(:,1),...
                    pitch_waveriness.deviation_from_flat.(multiple_pitch.parameters.syllables(i)).area_under_mean_subtracted_PC{j-multiple_pitch.parameters.days{3}+1},...
                    '.')

%                 %plots FF mean per day
%                 plot(median(multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).time_and_FF{j-multiple_pitch.parameters.days{3}+1}(:,1)),...
%                     multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).mean_FF{j-multiple_pitch.parameters.days{3}+1},...
%                     'ko',...
%                     'MarkerFaceColor','k',...
%                     'LineWidth',2)
%                 
%                 %plots FF sd per day
%                 errorbar(median(multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).time_and_FF{j-multiple_pitch.parameters.days{3}+1}(:,1)),...
%                     multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).mean_FF{j-multiple_pitch.parameters.days{3}+1},...
%                     multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).sd_FF{j-multiple_pitch.parameters.days{3}+1},...
%                     'ko',...
%                     'LineWidth',2)
%                 plot day mean
                
            catch err
                continue
            end
        end
        temp=cell2mat(pitch_waveriness.deviation_from_flat.(multiple_pitch.parameters.syllables(i)).day_mean_of_area_under_curve);
        figure; plot(1:length(temp),temp,'o'); ylim([0 60])
        
        % line for baseline pitch
        try
            line(get(gca, 'XLim'),...
                [pitch_waveriness.deviation_from_flat.(multiple_pitch.parameters.syllables(i)).baseline_mean...
                 pitch_waveriness.deviation_from_flat.(multiple_pitch.parameters.syllables(i)).baseline_mean],...
                'Color', [.3 .3 .3],...
                'LineStyle', '--')
        catch err
            continue
        end
        %line for end of baseline period
        try
            line([0-.25 0-.25],...
                get(gca, 'YLim'),...
                'Color', [0 0 0],...
                'LineStyle', '--')
        catch err
            continue
        end
        
%         % LINES demarcating lesion days (all are on end of day because
%         % generlaly have pre-lesion recordings for that same day, but not
%         % post)
%         for l=1:multiple_pitch.parameters.lesion_amount-1;
%             line_index(l)=multiple_pitch.parameters.lesion_dates_datenum{l+1}-multiple_pitch.parameters.lesion_dates_datenum{1};
%         end
%         for i=1:length(line_index);
%             try
%                 line([line_index(i)+0.25...
%                     line_index(i)]+0.25,...
%                     get(gca, 'YLim'),...
%                     'Color', [0 0 0],...
%                     'LineStyle', '--')
%             catch err
%                 continue
%             end
%         end
%         %SAVE, each fig for each syllable
%         saveas(figure(gcf), ['all_days_pitch_' multiple_pitch.parameters.phrase...
%             '/FF_over_time_' multiple_pitch.parameters.nameofbird '_' 'duration_'...
%             num2str(multiple_pitch.parameters.days{4}-multiple_pitch.parameters.days{3}+1)...
%             'days' '_' multiple_pitch.parameters.syllables(i)], 'fig')
        
    end
    

%% CHECK: plot random 
figure;
i=1;
j=5;
k=10;
plot(pitch_waveriness.deviation_from_flat.(multiple_pitch.parameters.syllables(i)).windowed_PC{jj}(:,k),'r');
pitch_waveriness.deviation_from_flat.(multiple_pitch.parameters.syllables(i)).mean_pitch{jj}(k)
figure;
plot(pitch_waveriness.deviation_from_flat.(multiple_pitch.parameters.syllables(i)).windowed_PC_minus_mean{jj}(:,k),'g');
%%
    
        
        
    

for day=1:10;
    for i=1:size(multiple_pitch.pc_all.t{day},2);
        mean_pitch_of_contour{day}(i)=mean(multiple_pitch.pc_all.t{day}(contour_start:contour_end,i),1);
        multiple_pitch.pc_all_TrialMeanSubtracted.t{day}(:,i)=multiple_pitch.pc_all.t{day}(:,i)-mean_pitch_of_contour{day}(i);
    end
    figure(day); plot(multiple_pitch.pc_all_TrialMeanSubtracted.t{day}(:,1:50)); title(['P.C. for day' num2str(day)]);
end

% finding STD of deviation from mean (of PC) at each timepoint (across
% renditions)

for day=1:10;
    for ii=1:size(multiple_pitch.pc_all.t{day},1);
        STDAcrossRenditions_OfDeviationFromMean_ofPC{day}(ii)=std(multiple_pitch.pc_all_TrialMeanSubtracted.t{day}(ii,:));
    end
     figure(day+3); plot(STDAcrossRenditions_OfDeviationFromMean_ofPC{day}(180:600)); title(['STD of P.C. across renditions for day' num2str(day)]);   
end