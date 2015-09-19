%% LT 1/21/14 - changed to not use findwnote2tw data but lt_db_get_label.  Latter better as does not repeat single song data.

%% LT 1/7/14 - getting gap duration for different sequences (use lt_all...various_v2 to get labels/onsets/offsets,
% and then this (better to do this in the context of lt_all_days_all_analysis

function [days, gap_durations_over_day]=lt_all_days_get_gap_durations_v2(transitions, date_and_time_of_various_struct, bird_folder);

%input:
% transitions; 
% date_and_time_of_various_struct % date and time that is in the file name of the all_days_all_various structure containing required data e.g. '06Jan2014_1737'} 
% bird_folder % home folder for the bird (path directory)

%output:
% gap_durations_over_day % per day, all durations in one array


%% LOAD data structure
load([bird_folder '/all_days_various_calculations/all_days_various_' date_and_time_of_various_struct])


%% CALCULATE GAP DURATIONS
num_days=size(all_days_various.lt_db_get_labels.labels,2);

for jj=1:num_days; % days
    if size(all_days_various.lt_db_get_labels.labels{jj},2)>0;
        for ii=1:size(all_days_various.lt_db_get_labels.labels{jj},2) % number of songs
            % GET DATA PER SONG            
            labels{jj}{ii}=all_days_various.lt_db_get_labels.labels{jj}{ii};
            onsets{jj}{ii}=all_days_various.lt_db_get_labels.onsets{jj}{ii};
            offsets{jj}{ii}=all_days_various.lt_db_get_labels.offsets{jj}{ii};
        end
    end
end
    

% FOR EACH SONG find the indices for occurances of the transitions and then find the
% gap duration.
for tr=1:size(transitions,2);
    for jj=1:num_days;
        if size(all_days_various.lt_db_get_labels.labels{jj},2)>0; % continue if day contains songs
            for ii=1:size(all_days_various.lt_db_get_labels.labels{jj},2); % iterate over songs
                transition_index.(transitions{tr}){jj}{ii}=findstr(labels{jj}{ii},transitions{tr}); % gives you the index of the first letter in the transition
                % for each occurence of transition, get onset of 2nd syl minus offset of 1st syl
                if transition_index.(transitions{tr}){jj}{ii}>0; % continue if there exist desired transitions in this song
                    for kk=1:size(transition_index.(transitions{tr}){jj}{ii},2);
                        gap_durations.(transitions{tr}){jj}{ii}(kk)=onsets{jj}{ii}((transition_index.(transitions{tr}){jj}{ii}(kk))+1)...
                            -offsets{jj}{ii}(transition_index.(transitions{tr}){jj}{ii}(kk));
                    end
                else
                    gap_durations.(transitions{tr}){jj}{ii}=[];
                end
            end
            
            % FOR EACH DAY, concatenate across songs
            dummy=[]; % new dummy array for each day
            for ii=1:size(gap_durations.(transitions{tr}){jj},2); % again, index for song #
                dummy=[dummy gap_durations.(transitions{tr}){jj}{ii}];
            end
            gap_durations_over_day{jj}.(transitions{tr})=dummy;
        end
    end
end

%%
days=all_days_various.parameters.days;


% find stats for each day
% 
% for tr=1:size(transitions);
%     
%     for jj=1:num_days;
%         if size(all_days_various.findwnote2tw.fvalsstr{jj},2)>0;
%             gap_durations_mean.(transitions{tr})(jj)=mean(gap_durations_over_day{jj}.(transitions{tr}));
%             gap_durations_STD.(transitions{tr})(jj)=std(gap_durations_over_day{jj}.(transitions{tr}));
% %         else
% %             gap_durations_mean.(transitions{tr})(7)=7;
% %             gap_durations_STD.(transitions{tr})(jj)=[];
% %             
%         end
%     end
% end
% % 
% %     % plot
% %     figure; hold on; title(['transition: ' transitions{tr} ', gap duration mean for each day.']);
%     scatter(gap_durations_mean.(transitions{tr}))
%     
% end
% 

end