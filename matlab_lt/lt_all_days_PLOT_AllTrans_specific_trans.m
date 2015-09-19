%% LT 5/15/14 modified from code written within lt_all_days_all_analysis_PLOT
% use this to plot any specific transitions (e.g. all div from a) from all_days_all_analysis
% structure, from lt_get_all_transitions data.
% instructions: run lt_all_days_all_various, choose lt_get_all_transitions,
% then run lt_all_days_all_analysis, and use that structure here as input.
% RUN THIS IN THE BIRD FOLDER
% INPUTS:
% con_div= 'convergence' or 'divergence'
% plot_all= 1(if plot all div or con) or 0
% syl_list= {'a','b'} (syl you want to plot contingent on (e.g. div from) IF plot_all is 0;
% ADAAindex = '28Apr2014_1419' (index for structure).

function lt_all_days_PLOT_AllTrans_specific_trans(con_div,plot_all,syl_list,ADAAindex)
% Load structure
cd all_days_all_analysis/
fname=dir(['all_days_all_analysis_' ADAAindex '*']);
load(fname.name)
cd ..


% INITIALIZE syllables
% extract all used syls
try
    syl_list_all=all_days_all_analysis.data{1}.all_songs.all_trans.syllable_list;
catch err; % error if there is no data for day 1
    syl_list_all=all_days_all_analysis.data{2}.all_songs.all_trans.syllable_list;
end

if plot_all==1;
    
    syl_listT=syl_list_all;
else
    syl_listT=syl_list;
end


% Make plots
% DIVERGENCE
if strcmp(con_div,'divergence')==1;
    syl1_list=syl_listT;
    for kk=1:length(syl1_list);
        syl1=syl1_list{kk};
        % syl2='b';
        syl2_list=syl_list_all; % plots div to all syl
        num_syl=length(syl2_list);
        % make syl colors
        syl_color=lt_make_plot_colors(length(syl2_list));
        figure;  hold on; title([con_div ' transition probability from ' syl1 ' to (see legend): ']);
        for k=1:length(syl2_list);
            syl2=syl2_list{k};
            dummy=[];
            for i=1:size(all_days_all_analysis.data,2);
                try
                    dummy=[dummy; [i all_days_all_analysis.data{i}.all_songs.all_trans.(syl1).transition_to_.(syl2).relative_fraction_of_trans]];
                catch err
                    continue
                end
            end
            %     plot(dummy,'Marker','o','Color',syl_color{k},'MarkerEdgeColor','black','MarkerFaceColor',syl_color{k},'MarkerSize',8);
            hplot(k)=plot(dummy(:,1), dummy(:,2),'Marker','o','Color',syl_color{k},'MarkerEdgeColor','black','MarkerFaceColor',syl_color{k},'MarkerSize',8);
        end
        legend(hplot,[syl2_list]);
        ylim([0 1]); xlabel(['days (day 1 = ' all_days_all_analysis.parameters.first_day]); ylabel('probability of transition')
    end
    
elseif strcmp(con_div,'convergence')==1;
    % CONVERGENCE
    syl2_list=syl_listT;
    for kk=1:length(syl2_list)
        syl2=syl2_list{kk};
        syl1_list=syl_list_all;
        num_syl=length(syl1_list);
        % make syl colors
        syl_color=lt_make_plot_colors(length(syl1_list));
        figure;  hold on; title(['Convergence transition probability to ' syl2 ' from (see legend): ']);
        for k=1:length(syl1_list);
            syl1=syl1_list{k};
            dummy = [];
            for i=1:size(all_days_all_analysis.data,2);
                try
                    dummy = [dummy; [i all_days_all_analysis.data{i}.all_songs.all_trans.convergent_to.(syl2).from.(syl1).fraction_of_trans]];
                catch err
                    continue
                end
            end
            hplot2(k)= plot(dummy(:,1), dummy(:,2),'Marker','o','Color',syl_color{k},'MarkerEdgeColor','black','MarkerFaceColor',syl_color{k},'MarkerSize',8);
        end
        legend(hplot2,[syl1_list]);
        ylim([0 1]); xlabel(['days (day 1 = ' all_days_all_analysis.parameters.first_day]); ylabel('probability of transition')
    end
end
end