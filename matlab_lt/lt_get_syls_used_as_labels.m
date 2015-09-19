%% LT 7/7 - Given batch file, spits out all alphabateical symbols used as label in a single day
% Run from within day
% e.g. syl_list = lt_get_syls_used_as_labels(batch_file);
% Note: variable names a little weird because this was extracted from
% lt_get_all_transitions.

function syl_list_final = lt_get_syls_used_as_labels(batch_file);

%% RUN

[syl_order_temp_WithDashes syl_times_temp song_filenames syl_onsets_temp syl_offsets_temp]=lt_db_get_labels(batch_file);

for j=1:length(syl_order_temp_WithDashes);
        syl_order_temp{j}=syl_order_temp_WithDashes{j}(syl_order_temp_WithDashes{j}~='-'); % removes nonlabeled syllables (i.e. "-' or spacebar);
        syl_order_temp{j}=syl_order_temp{j}(syl_order_temp{j}~=' '); 
        
        syl_times{j}=syl_times_temp{j}(syl_order_temp_WithDashes{j}~='-');
        syl_times{j}=syl_times_temp{j}(syl_order_temp{j}~=' ');

        syl_onsets{j}=syl_onsets_temp{j}(syl_order_temp_WithDashes{j}~='-');
        syl_onsets{j}=syl_onsets_temp{j}(syl_order_temp{j}~=' ');
        
        syl_offsets{j}=syl_offsets_temp{j}(syl_order_temp_WithDashes{j}~='-');
        syl_offsets{j}=syl_offsets_temp{j}(syl_order_temp{j}~=' ');
end

clock=1;
for j=1:length(syl_order_temp);
    if length(syl_order_temp{j})>0;  % removes day if there are no labeled syl
        all_trans.syl_order{clock}=syl_order_temp{j};    
        all_trans.syl_order_WithDashes{clock}=syl_order_temp_WithDashes{j};
        % for this day, adds data on syl absolute times, onsets, and
        % offsets, and song filename
        all_trans.syl_times{clock}=syl_times{j};
        all_trans.syl_onsets{clock}=syl_onsets{j};
        all_trans.syl_offsets{clock}=syl_offsets{j};
        all_trans.song_filenames{clock}=song_filenames{j};
        clock=clock+1;
    end
end

% puts all syl sang that day into one long vector and gets a list of all
% letters used as syll.
    all_trans.syl_order_compiled_across_songs=cell2mat(all_trans.syl_order);
%     disp('getting all labels used.  this only works if you only used lower and upper case alphabet letter.');
    alphabet={'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z',...
        'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'};

clock=1;
for i=1:length(alphabet);
    instances_temp=strfind(all_trans.syl_order_compiled_across_songs,alphabet{i});
    if sum(instances_temp)>0;
        all_trans.syllable_list{clock}=alphabet{i};
        clock=clock+1;
    elseif sum(instances_temp==0);
    end
    clear instances_temp
end

syl_list_final=all_trans.syllable_list;


end
