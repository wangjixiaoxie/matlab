%% LT 7/23/14 - Runs lt_all_days_various... for both normal syllables (e.g. a, or ab) and motifs (e.g. abc[c] or dgg[g]g). 
% THIS basically is a compact way to run both
% lt_all_days_various_calculations_FUNCTION and
% lt_all_days_various_calculations_FUNCTION_MOTIF together

% OUTPUT: One combined structure
% Enter the variables by hand below.

%% GENERAL INPUTS
% phrase = 'SyntaxDepPitchShift';
% first_day= '05Jul2014';
% last_day= '23Jul2014';
% save_results=1;
% 
% % functions to run (SAME FOR ALL MOTIFS)
% FcnAll={'check_HTF','calcdaypitch'};


%% INPUTS specifically for motif switches:
% i.e. each entry is for a separate motif (if want to get both ab(b) and
% ab(c), then syl_targ_list would be {'b','c'}). function will run
% separately for each motif.

% motif_pre_list={'dcc','dccb','ab','bcc','bccb'};
% syl_targ_list={'b','b','b','b','b'};
% motif_post_list={'','','','',''};
% syl_replace_list={'M','N','Z','P','Q'}; % important, labels will temporarily be changed to this. has to be unused (e.g. all caps).

%% INPUTS FOR analysis (lt_all_days_various)
% enter as normally would, for each syl target (each cell entry is for separate motif switch)

% lt_calc_pitch...:
% freq_rangeDP_list={{[3050 3750]},{[3050 3750]},{[3050 3750]},{[3050 3750]},{[3050 3750]}}; % i.e. first cell is for first replaced syl
% pc_time_windowDP_list={{[60 150]},{[60 150]},{[60 150]},{[55 150]},{[60 150]}};
% plot_resultDP=0;
% pc_windowDP=2000;
% 
% % lt_check_HTF...:
% freq_rangeHTF={{[3050 3750]},{[3050 3750]},{[3050 3750]},{[3050 3750]},{[3050 3750]}};
% evtaf_verHTF='v4';
% get_WN_hitsHTF=1;
% get_offline_matchHTF=1;
% get_FFHTF=1;
% if get_offline_matchHTF==1; % template params if also doing offline check
%     template_name='pu11wh87ab_SyntaxDepPitchShift_v1_6';
%     cntrng_values{1}={[5 14 1] 'or' 'n' 'y'};
%     cntrng_values{2}={[5 14 0.1] 'or' 'n' 'y'};
%     cntrng_values{3}={[5 14 1.1] 'and' 'n' 'y'};
%     cntrng_values{4}={[1 4 1.5] 'or' 'n' 'n'};
%     cntrng_values{5}={[1 4 1.5] 'or' 'n' 'n'};
%     cntrng_values{6}={[1 4 1.5] 'or' 'n' 'n'};
%     col_logic='(a+b+c)*(d+e+f)';
% end

%% INPUT: WITHIN FUNCTION VARIABLES SETTING (DIFFERENT FOR ALL MOTIFS)
% DO not change the values within each function, but add and subtract V's if you need.

% num_motifs=size(syl_targ_list,2);
% for i=1:num_motifs;
%     % lt_check_HTF...:
%     if get_offline_matchHTF==1;
%         V1_mot{i}={'sylHTF',{syl_replace_list{i}},'syl_preHTF',{''},'syl_postHTF',{''},'freq_rangeHTF',freq_rangeHTF{i},'evtaf_verHTF',evtaf_verHTF,'get_WN_hitsHTF',get_WN_hitsHTF,'get_offline_matchHTF',...
%             get_offline_matchHTF,'get_FFHTF',get_FFHTF,'template_name',template_name,'cntrng_values',cntrng_values,'col_logic',col_logic};
%     else
%         V1_mot{i}={'sylHTF',{syl_replace_list{i}},'syl_preHTF',{''},'syl_postHTF',{''},'freq_rangeHTF',freq_rangeHTF{i},'evtaf_verHTF',evtaf_verHTF,'get_WN_hitsHTF',get_WN_hitsHTF,'get_offline_matchHTF',...
%             get_offline_matchHTF,'get_FFHTF',get_FFHTF};
%     end
%     
%     % lt_calc_pitch...:
%     V2_mot{i}={'syl_targetDP',syl_replace_list{i}, 'syl_preDP',{''},'phraseDP',phrase,'freq_rangeDP',freq_rangeDP_list{i},'pc_time_windowDP',pc_time_windowDP_list{i},'plot_resultDP',plot_resultDP,'pc_windowDP',pc_windowDP};
%     
%     % combine
%     Vall_mot{i}=[V1_mot{i},V2_mot{i}];
% end
% 
% 


%% INPUTS FOR REGULAR (NON_MOTIF) SYLLABLES:

% lt_check_pitch:
% syl_preHTF={'a'};
% sylHTF={'b'};
% syl_postHTF={''};
% freq_rangeHTF={[3050 4200]};
% 
% % lt_calc_pitch
% syl_preDP={'a'};
% syl_targetDP='b';
% freq_rangeDP={[3050 4200]};
% pc_time_windowDP={[50 150]};
% 
% 
% % Within function VARIABLES SETTING:
% % lt_check...:
% if get_offline_matchHTF==1;
%     V1={'sylHTF',sylHTF,'syl_preHTF',syl_preHTF,'syl_postHTF',syl_postHTF,'freq_rangeHTF',freq_rangeHTF,'evtaf_verHTF',evtaf_verHTF,'get_WN_hitsHTF',get_WN_hitsHTF,'get_offline_matchHTF',get_offline_matchHTF,'get_FFHTF',get_FFHTF,'template_name',template_name,...
%         'cntrng_values',cntrng_values,'col_logic',col_logic};
% else
%     V1={'sylHTF',sylHTF,'syl_preHTF',syl_preHTF,'syl_postHTF',syl_postHTF,'freq_rangeHTF',freq_rangeHTF,'evtaf_verHTF',evtaf_verHTF,'get_WN_hitsHTF',get_WN_hitsHTF,'get_offline_matchHTF',get_offline_matchHTF,'get_FFHTF',get_FFHTF};
% end
% % lt_calc...:
% V2={'syl_targetDP',syl_targetDP, 'syl_preDP',syl_preDP,'phraseDP',phrase,'freq_rangeDP',freq_rangeDP,'pc_time_windowDP',pc_time_windowDP,'plot_resultDP',plot_resultDP,'pc_windowDP',pc_windowDP};
% Vall=[V1,V2];


%% RUN ANALYSIS - part 1: MOTIFS separately
% First ask about motifs, so mistakes won't lead to permanent changes to
% batch files;
motif_counter=0;
for i=1:num_motifs;
    disp([motif_pre_list{i} syl_targ_list{i} motif_post_list{i} ' will be temporarily changed to '...
        motif_pre_list{i} syl_replace_list{i} motif_post_list{i} ' in the following analysis.']);
    if strcmp(input('type y to continue, anything else to quit: ','s'),'y')==0;
        forced2quit; % force quit
    end
end

% Second, run function once for each motif.
for i=1:num_motifs;
    MOTIF_params={syl_targ_list{i},motif_pre_list{i},motif_post_list{i},syl_replace_list{i}};
    
    % RUN
    [filename_save, all_days_various(motif_counter+1)]=lt_all_days_various_calculations_FUNCTION_MOTIF(phrase,first_day,last_day,FcnAll, Vall_mot{i},save_results,MOTIF_params);
    motif_counter=motif_counter+1;
    
    % OUTPUT
    motif_table{motif_counter,1}=[motif_pre_list{i} syl_targ_list{i} motif_post_list{i}];
    motif_table{motif_counter,2}=[motif_pre_list{i} syl_replace_list{i} motif_post_list{i}];
end


%% RUN ANALYSIS- part 2, regular syllables, will only run once.

% RUN
% [filename_save, all_days_various(motif_counter+1)]=lt_all_days_various_calculations_FUNCTION(phrase, first_day, last_day, FcnAll, Vall, save_results);
% motif_counter=motif_counter+1;
% 
% for i=1:length(syl_targetDP);
%     motif_table{motif_counter+i-1,1}=[syl_preDP{i} syl_targetDP(i)];
% end


%% COMBINE DATA STRUCTURES (i.e. the motif ADV (multiple, one for each motif), and regular ADV (one containing all syls));
bird_dir=['/bluejay3/lucas/birds/' all_days_various(1).parameters.nameofbird '/'];
ADV_dir=[bird_dir 'all_days_various_calculations/'];
all_days_various_COMPILED.all_days_various=all_days_various;

% GET list of motifs their structure index  - for lt_calc_day_pitch only
c=1;
try
for i=1:size(all_days_various_COMPILED.all_days_various,2);
    try % in case no data in first day.
    fields=fieldnames(all_days_various_COMPILED.all_days_various(i).lt_calc_day_pitch{1}.FF);
    catch err
    fields=fieldnames(all_days_various_COMPILED.all_days_various(i).lt_calc_day_pitch{2}.FF);
    end        
        num_fields=length(fields);
    for ii=1:num_fields;
        syllables_COMPILED{c,1}=i;
        syllables_COMPILED{c,2}=fields{ii};
        c=c+1;
    end
end
parameters.syllables_COMPILED=syllables_COMPILED;

catch err
end

% OUTPUT
% compiled data structure
parameters.motif_table=motif_table;
timestamp=lt_get_timestamp(0);
parameters.first_day=first_day;
parameters.last_day=last_day;
parameters.time_of_analysis=timestamp;
savename=[ADV_dir 'all_days_various_MOTIF_COMPILED_' timestamp];
things_analyzed=all_days_various_COMPILED.all_days_various(1).parameters.things_analyzed;
parameters.things_analyzed=things_analyzed;
parameters.savename=savename;
parameters.birdname=all_days_various(1).parameters.nameofbird;
all_days_various_COMPILED.parameters_COMPILATION=parameters;

save(savename,'all_days_various_COMPILED');

% add to log.txt (assumes that all motifs underwent saem analyses)
cd(ADV_dir);
log_note=['MOTIF_COMPILED_' timestamp '__|' things_analyzed '|__' strjoin(parameters.motif_table(:,2)')  '|__' first_day '_to_' last_day];
fid=fopen('log.txt','a');
fprintf(fid,'%s \n \n',log_note);
fclose(fid);
