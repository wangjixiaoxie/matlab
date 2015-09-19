%% LT 5/19/14 - use all_days_all_analysis structure. Run this in the bird folder.
% gives plots of 1) how effective was your template over days? 2) how well
% did the bird learn to avoid WN (per song and per syl sang).
% NEEDS: lt_get_all_trans... and triglabel data.

function lt_all_days_PLOT_WN_stats(ADAAind)

all_days_all_analysis=lt_load_all_days_all_analysis(ADAAind);
syls=fieldnames(all_days_all_analysis.data{1}.WN_hit_rate.per_song);
syl=syls{1};

for i=1:size(all_days_all_analysis.data,2);
    num_songs=size(all_days_all_analysis.data{i}.WN_hit_rate.per_song.(syl).hits_labels_allhits,1);
    if num_songs>0;
        hit_rate_Mean_PerSong(i)=all_days_all_analysis.data{i}.WN_hit_rate.all_songs.(syl).hits_labels_allhits(1)/num_songs;
        rend_rate_Mean_PerSong(i)=all_days_all_analysis.data{i}.WN_hit_rate.all_songs.(syl).hits_labels_allhits(2)/num_songs;
        song_length_Mean(i)=length(all_days_all_analysis.data{i}.all_songs.all_trans.syl_order_compiled_across_songs)/num_songs;
        
        hit_rate_Mean_PerSyl(i)=hit_rate_Mean_PerSong(i)/song_length_Mean(i);
        rend_rate_Mean_PerSyl(i)=rend_rate_Mean_PerSong(i)/song_length_Mean(i);
        
    elseif num_songs==0;
        hit_rate_Mean_PerSong(i)=nan;
        rend_rate_Mean_PerSong(i)=nan;
        song_length_Mean(i)=nan;
        
        hit_rate_Mean_PerSyl(i)=nan;
        rend_rate_Mean_PerSyl(i)=nan;
    end
end

figure; hold on;
plot(hit_rate_Mean_PerSong,'sr'); 
plot(rend_rate_Mean_PerSong,'ob');
title('Mean per song: RED: hits/song; BLUE: renditions of target/song; (all catch)');

figure; hold on;
plot(hit_rate_Mean_PerSyl,'sr'); 
plot(rend_rate_Mean_PerSyl,'ob');
title('Mean per all syllables sang: RED: hits/syl; BLUE: renditions of target/syl; (all catch)');

figure; hold on;
plot(song_length_Mean,'or');
title('mean song length across days'); 
ylim([0 max(song_length_Mean)+5])




% plotting WN hit rate over all things, see bk51bk59 analysis script for details.

syllables=all_days_all_analysis.parameters.WN_trig.trig_syl;
for i=1:length(syllables);
    figure; hold on; title(['Template accuracy: ' syllables(i) ' green= WN hits/Labels; black= FalsePos/TotalWN'])
    for j=1:length(all_days_all_analysis.data); % num of days
        try
            scatter(j,all_days_all_analysis.data{j}.WN_hit_rate.all_songs.(syllables(i)).hits_divide_labels,'g','MarkerFaceColor','g');
            scatter(j,all_days_all_analysis.data{j}.WN_hit_rate.all_songs.(syllables(i)).FalsePos_divide_TotalWN,'k','MarkerFaceColor','k')
            %         scatter(j,all_days_all_analysis.data{j}.WN_hit_rate.all_songs.a.FalsePos/all_days_all_analysis.data{j}.WN_hit_rate.all_songs.a.hits_labels_allhits(1),'b','MarkerFaceColor','b')
        catch err
        end
    end
end

