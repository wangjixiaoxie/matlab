function [] = lt_get_all_transition_probabilities_FUNCTION(batch,bird_name,phrase_folder,phrase_filename,date_expt);
%% LT 5/6 updated to include onsets, offsets, etc.


% LT 3/19 - converted script to function.
% LT 10/2013 - takes labeled songs (and a batch file) and pulls out 1)
% all letters that have been used as labels and 2) calculates fraction of
% all possible transitions in those songs (output is transition matrix).  

% batch_question=input('do you want to use batch.catch.keep? (y or n) ', 's');
% if batch_question=='y';
%     all_trans.parameters.batchfile='batch.catch.keep';
% elseif batch_question=='n';
%     all_trans.parameters.batchfile=input('what is the name of the batch file? ','s');
% end
% clear batch_question;
all_trans.parameters.batchfile=batch;

% computer=input('which bluejay computer? ');
% all_trans.parameters.computer=['/bluejay' num2str(computer) '/lucas/birds'];

curr_day_dir=pwd;
cd ../../
all_trans.parameters.computer=pwd;
cd(curr_day_dir);


% all_trans.parameters.bird_name=input('What is the name of your bird? ', 's');
% all_trans.parameters.phrase=input('Phrase (for save folder)?  ','s');
% all_trans.parameters.filename_marker=input('Label to add to filename (e.g. to identify different batches within same day)? ','s');
% all_trans.parameters.date_of_exp=input('what is the date of your experiment?\n(ex: 03Aug2012)  ', 's');
all_trans.parameters.bird_name=bird_name;
all_trans.parameters.phrase=phrase_folder;
all_trans.parameters.filename_marker=phrase_filename;
all_trans.parameters.date_of_exp=date_expt;


%% MAKING a list of all letters used for syllables

% syl_order_temp_WithDashes=getlabels(all_trans.parameters.batchfile);
[syl_order_temp_WithDashes syl_times_temp song_filenames syl_onsets_temp syl_offsets_temp]=lt_db_get_labels(all_trans.parameters.batchfile);

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

% NOTE: This is obsolete. I used to get transitions using this, but now I
% do for each individual song and then combine across songs.  It is still
% useful for getting list of alphabets.

% puts all syl sang that day into one long vector and gets a list of all
% letters used as syll.
    all_trans.syl_order_compiled_across_songs=cell2mat(all_trans.syl_order);
    disp('getting all labels used.  this only works if you only used lower and upper case alphabet letter. type return if fine')
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

%% Getting all syllables that occur in conv and div transitions.
% for each syllable, get a list of all the syllables that follow it.
for i=1:length(all_trans.syllable_list);
    for j=1:length(all_trans.syl_order); % j = song #
        try
            all_trans.(all_trans.syllable_list{i}).index_of_syl{j}=strfind(all_trans.syl_order{j},all_trans.syllable_list{i}); % this gets the indices of all instances of the syl.
            all_trans.(all_trans.syllable_list{i}).all_subsequent_syl{j}=all_trans.syl_order{j}(all_trans.(all_trans.syllable_list{i}).index_of_syl{j}+1); %this gets the lable of all syllables that follow the syl of interest.
        catch % have this try-catch here to override error when the syl is the last syl of song (and thus no syl follow it)
            all_trans.(all_trans.syllable_list{i}).all_subsequent_syl{j}=all_trans.syl_order{j}(all_trans.(all_trans.syllable_list{i}).index_of_syl{j}(1:end-1)+1);
            continue
        end
    end
    for j=1:length(all_trans.syl_order); % j = song #
        try
            all_trans.(all_trans.syllable_list{i}).all_preceding_syl{j}=all_trans.syl_order{j}(all_trans.(all_trans.syllable_list{i}).index_of_syl{j}-1); %this gets the lable of all syllables that precede the syl of interest.
        catch % have this try-catch here to override error when the syl is the first syl of song (and thus no syl precede it)
            all_trans.(all_trans.syllable_list{i}).all_preceding_syl{j}=all_trans.syl_order{j}(all_trans.(all_trans.syllable_list{i}).index_of_syl{j}(2:end)-1);
            continue
        end
    end
    all_trans.(all_trans.syllable_list{i}).all_subseq_syl_compiled_over_songs=cell2mat(all_trans.(all_trans.syllable_list{i}).all_subsequent_syl);
    all_trans.(all_trans.syllable_list{i}).all_preceding_syl_compiled_over_songs=cell2mat(all_trans.(all_trans.syllable_list{i}).all_preceding_syl);
end

%% Getting probabilities.
% DIVERGERNT: For each syllable, gets the proportion of transitions that go to each one
% of the other syllables.   
for j=1:length(all_trans.syl_order); % number of songs.
    for i=1:length(all_trans.syllable_list); %first syllable
        all_trans.(all_trans.syllable_list{i}).total_number_of_transitions_persong{j}=length(all_trans.(all_trans.syllable_list{i}).all_subsequent_syl{j});%total number of div transitions recorded for each syl
        for ii=1:length(all_trans.syllable_list); % second syllable
            all_trans.(all_trans.syllable_list{i}).transition_to_.(all_trans.syllable_list{ii}).indices_persong{j}=strfind(all_trans.(all_trans.syllable_list{i}).all_subsequent_syl{j}, all_trans.syllable_list{ii}); %find the indices in the compiled list of the 2nd desired syl
            all_trans.(all_trans.syllable_list{i}).transition_to_.(all_trans.syllable_list{ii}).number_of_instances_persong{j}=length(all_trans.(all_trans.syllable_list{i}).transition_to_.(all_trans.syllable_list{ii}).indices_persong{j}); % counts the number of instances of this trans.
            all_trans.(all_trans.syllable_list{i}).transition_to_.(all_trans.syllable_list{ii}).relative_fraction_of_trans{j}=all_trans.(all_trans.syllable_list{i}).transition_to_.(all_trans.syllable_list{ii}).number_of_instances_persong{j}/all_trans.(all_trans.syllable_list{i}).total_number_of_transitions_persong{j}; % proportion of total trans.
            
        end
    end
end

% CONVERGENT: For each syllable, gets the proportion of transitions that
% come from each otehr syllable.
for j=1:length(all_trans.syl_order); % number of songs.
    for i=1:length(all_trans.syllable_list); %possible syllables
        all_trans.convergent_to.(all_trans.syllable_list{i}).total_number_of_transitions_persong{j}=length(all_trans.(all_trans.syllable_list{i}).all_preceding_syl{j});%total number of div transitions recorded for each syl
        for ii=1:length(all_trans.syllable_list); % preceding syllables
            all_trans.convergent_to.(all_trans.syllable_list{i}).from.(all_trans.syllable_list{ii}).indices_persong{j}=strfind(all_trans.(all_trans.syllable_list{i}).all_preceding_syl{j}, all_trans.syllable_list{ii}); %find the indices in the compiled list of the preceding syl
            all_trans.convergent_to.(all_trans.syllable_list{i}).from.(all_trans.syllable_list{ii}).number_of_instances_persong{j}=length(all_trans.convergent_to.(all_trans.syllable_list{i}).from.(all_trans.syllable_list{ii}).indices_persong{j}); % counts the number of instances of this trans.
            all_trans.convergent_to.(all_trans.syllable_list{i}).from.(all_trans.syllable_list{ii}).fraction_of_trans{j}=all_trans.convergent_to.(all_trans.syllable_list{i}).from.(all_trans.syllable_list{ii}).number_of_instances_persong{j}/all_trans.convergent_to.(all_trans.syllable_list{i}).total_number_of_transitions_persong{j}; % proportion of total trans.
            
        end
    end
end


%% compile data across songs to get day measure of transition probabilities.
% DIVERGENT
for i=1:length(all_trans.syllable_list); %first syllable
    all_trans.(all_trans.syllable_list{i}).total_number_of_transitions=sum(cell2mat(all_trans.(all_trans.syllable_list{i}).total_number_of_transitions_persong));%total number of div transitions recorded for each syl
    for ii=1:length(all_trans.syllable_list); % second syllable
%         all_trans.(all_trans.syllable_list{i}).transition_to_.(all_trans.syllable_list{ii}).indices_for_compiled_list=strfind(all_trans.(all_trans.syllable_list{i}).all_subseq_syl_compiled_over_songs, all_trans.syllable_list{ii}); %find the indices in the compiled list of the 2nd desired syl
        all_trans.(all_trans.syllable_list{i}).transition_to_.(all_trans.syllable_list{ii}).number_of_instances=sum(cell2mat(all_trans.(all_trans.syllable_list{i}).transition_to_.(all_trans.syllable_list{ii}).number_of_instances_persong)); % counts the number of instances of this trans.
        all_trans.(all_trans.syllable_list{i}).transition_to_.(all_trans.syllable_list{ii}).relative_fraction_of_trans=all_trans.(all_trans.syllable_list{i}).transition_to_.(all_trans.syllable_list{ii}).number_of_instances/all_trans.(all_trans.syllable_list{i}).total_number_of_transitions; % proportion of total trans.
    end
end

% CONVERGENT.
for i=1:length(all_trans.syllable_list); %first syllable
    all_trans.convergent_to.(all_trans.syllable_list{i}).total_number_of_transitions=sum(cell2mat(all_trans.convergent_to.(all_trans.syllable_list{i}).total_number_of_transitions_persong));%total number of div transitions recorded for each syl
    for ii=1:length(all_trans.syllable_list); % second syllable
%         all_trans.(all_trans.syllable_list{i}).transition_to_.(all_trans.syllable_list{ii}).indices_for_compiled_list=strfind(all_trans.(all_trans.syllable_list{i}).all_subseq_syl_compiled_over_songs, all_trans.syllable_list{ii}); %find the indices in the compiled list of the 2nd desired syl
        all_trans.convergent_to.(all_trans.syllable_list{i}).from.(all_trans.syllable_list{ii}).number_of_instances=sum(cell2mat(all_trans.convergent_to.(all_trans.syllable_list{i}).from.(all_trans.syllable_list{ii}).number_of_instances_persong)); % counts the number of instances of this trans.
        all_trans.convergent_to.(all_trans.syllable_list{i}).from.(all_trans.syllable_list{ii}).fraction_of_trans=all_trans.convergent_to.(all_trans.syllable_list{i}).from.(all_trans.syllable_list{ii}).number_of_instances/all_trans.convergent_to.(all_trans.syllable_list{i}).total_number_of_transitions; % proportion of total trans.
    end
end

% BELOW PROBLEM: did with compiled syl across songs. should do each song
% individually. 
% for i=1:length(all_trans.syllable_list); %first syllable
%     all_trans.(all_trans.syllable_list{i}).total_number_of_transitions=length(all_trans.(all_trans.syllable_list{i}).all_subseq_syl_compiled_over_songs);%total number of div transitions recorded for each syl
%     for ii=1:length(all_trans.syllable_list); % second syllable
%         all_trans.(all_trans.syllable_list{i}).transition_to_.(all_trans.syllable_list{ii}).indices_for_compiled_list=strfind(all_trans.(all_trans.syllable_list{i}).all_subseq_syl_compiled_over_songs, all_trans.syllable_list{ii}); %find the indices in the compiled list of the 2nd desired syl
%         all_trans.(all_trans.syllable_list{i}).transition_to_.(all_trans.syllable_list{ii}).number_of_instances=length(all_trans.(all_trans.syllable_list{i}).transition_to_.(all_trans.syllable_list{ii}).indices_for_compiled_list); % counts the number of instances of this trans.
%         all_trans.(all_trans.syllable_list{i}).transition_to_.(all_trans.syllable_list{ii}).relative_fraction_of_trans=all_trans.(all_trans.syllable_list{i}).transition_to_.(all_trans.syllable_list{ii}).number_of_instances/all_trans.(all_trans.syllable_list{i}).total_number_of_transitions; % proportion of total trans.
%     end
% end

% make a "summary" field with an array
% convert that array to a csv file.
all_trans.summary.syl_labels_in_order=all_trans.syllable_list;
all_trans.summary.legend='1st dim=1st syl; 2nd dim=2nd syl; ';

for i=1:length(all_trans.syllable_list); %first syllable
    for ii=1:length(all_trans.syllable_list); % second syllable
        all_trans.summary.divergence.matrix_of_fractions(i,ii)=all_trans.(all_trans.syllable_list{i}).transition_to_.(all_trans.syllable_list{ii}).relative_fraction_of_trans;
        all_trans.summary.divergence.matrix_of_amounts(i,ii)=all_trans.(all_trans.syllable_list{i}).transition_to_.(all_trans.syllable_list{ii}).number_of_instances;
        
        all_trans.summary.convergence.matrix_of_fractions(i,ii)=all_trans.convergent_to.(all_trans.syllable_list{ii}).from.(all_trans.syllable_list{i}).fraction_of_trans;
        all_trans.summary.convergence.matrix_of_amounts(i,ii)=all_trans.convergent_to.(all_trans.syllable_list{ii}).from.(all_trans.syllable_list{i}).number_of_instances;
        
    end

end
all_trans.summary.matrix_of_fraction_of_all_trans=all_trans.summary.divergence.matrix_of_amounts./(sum(squeeze(sum(all_trans.summary.divergence.matrix_of_amounts))));
% all_trans.summary.matrix_of_fraction_of_all_trans=all_trans.summary.convergence.matrix_of_amounts./(sum(squeeze(sum(all_trans.summary.convergence.matrix_of_amounts))));

%         do this also per song
%         for all songs
%             repeat as above

%% plotting heat map of transition matrix

% COMMENTED OUT - old version, not good
% figure; hold on
% 
% subplot(1,3,1); imagesc(all_trans.summary.divergence.matrix_of_fractions)
% labels_temp=cell2mat(all_trans.summary.syl_labels_in_order);
% title({['divergence transition matrix on ' all_trans.parameters.date_of_exp '. V-axis=Syl1; H-axis=Syl2.'],['labels: ' labels_temp]})
% 
% subplot(1,3,2); imagesc(all_trans.summary.convergence.matrix_of_fractions)
% labels_temp=cell2mat(all_trans.summary.syl_labels_in_order);
% % title({['convergence transition matrix on ' all_trans.parameters.date_of_exp '. V-axis=Syl1; H-axis=Syl2.'],['labels: ' labels_temp]})
% title('convergence transition matrix')
% 
% subplot(1,3,3); imagesc(all_trans.summary.matrix_of_fraction_of_all_trans)
% labels_temp=cell2mat(all_trans.summary.syl_labels_in_order);
% % title({['total transition matrix (same for con and div) on ' all_trans.parameters.date_of_exp '. V-axis=Syl1; H-axis=Syl2.'],['labels: ' labels_temp]})
% title('total transition matrix (same for con and div)')

figure; hold on;

%divergent fractions
subplot(1,3,1); 
syl_mat = 100*all_trans.summary.divergence.matrix_of_fractions;
syl_num=size(all_trans.summary.divergence.matrix_of_fractions,2);
syl_labels=all_trans.summary.syl_labels_in_order;

imagesc(syl_mat);
colormap(flipud(gray))
textStrings = num2str(syl_mat(:),'%0.0f');
textStrings = strtrim(cellstr(textStrings));
[x,y]=meshgrid(1:syl_num);
hStrings = text(x(:),y(:),textStrings(:), 'HorizontalAlignment','center');
midValue = mean(get(gca,'CLim'));
textColors = repmat(syl_mat(:) > midValue,1,3);
set(hStrings,{'Color'},num2cell(textColors,2));
set(gca,'XTick',1:syl_num,'XTickLabel',syl_labels, 'YTick',1:syl_num,...
    'YTickLabel',syl_labels, 'TickLength',[0 0]);
ylabel('transition from:')
xlabel('transition to:')
title('divergent fractions');


%convergent fractions
subplot(1,3,2); 
syl_mat_con = 100*all_trans.summary.convergence.matrix_of_fractions;

imagesc(syl_mat_con);
colormap(flipud(gray))
textStrings = num2str(syl_mat_con(:),'%0.0f');
textStrings = strtrim(cellstr(textStrings));
[x,y]=meshgrid(1:syl_num);
hStrings = text(x(:),y(:),textStrings(:), 'HorizontalAlignment','center');
midValue = mean(get(gca,'CLim'));
textColors = repmat(syl_mat_con(:) > midValue,1,3);
set(hStrings,{'Color'},num2cell(textColors,2));
set(gca,'XTick',1:syl_num,'XTickLabel',syl_labels, 'YTick',1:syl_num,...
    'YTickLabel',syl_labels, 'TickLength',[0 0]);
ylabel('transition from:')
xlabel('transition to:')
title('convergent fractions');


%absolute number of transitions
subplot(1,3,3); 
syl_mat_amounts = all_trans.summary.convergence.matrix_of_amounts;

imagesc(syl_mat_amounts);
colormap(flipud(gray))
textStrings = num2str(syl_mat_amounts(:),'%0.0f');
textStrings = strtrim(cellstr(textStrings));
[x,y]=meshgrid(1:syl_num);
hStrings = text(x(:),y(:),textStrings(:), 'HorizontalAlignment','center');
midValue = mean(get(gca,'CLim'));
textColors = repmat(syl_mat_amounts(:) > midValue,1,3);
set(hStrings,{'Color'},num2cell(textColors,2));
set(gca,'XTick',1:syl_num,'XTickLabel',syl_labels, 'YTick',1:syl_num,...
    'YTickLabel',syl_labels, 'TickLength',[0 0]);
ylabel('transition from:')
xlabel('transition to:')
title('number of transitions');


%% Saving

if exist([all_trans.parameters.computer '/' all_trans.parameters.bird_name '/all_days_transition_matrix_' all_trans.parameters.phrase], 'dir') ~= 7
    mkdir([all_trans.parameters.computer '/' all_trans.parameters.bird_name '/all_days_transition_matrix_' all_trans.parameters.phrase]);
else
end

save([all_trans.parameters.computer '/' all_trans.parameters.bird_name '/all_days_transition_matrix_' all_trans.parameters.phrase '/'...
    all_trans.parameters.bird_name '_' all_trans.parameters.date_of_exp '_' all_trans.parameters.filename_marker], 'all_trans');

saveas(figure(gcf),[all_trans.parameters.computer '/' all_trans.parameters.bird_name '/all_days_transition_matrix_' all_trans.parameters.phrase '/three_transition_matrices_for' all_trans.parameters.date_of_exp '_NEW'],'fig')

% export transition matrix as csv for excel
% xlswrite([all_trans.parameters.computer '/' all_trans.parameters.bird_name '/all_days_transition_matrix_' all_trans.parameters.phrase '/trans_matrix_' all_trans.parameters.bird_name '_' all_trans.parameters.date_of_exp],all_trans.summary.transition_matrix.matrix_of_fraction_of_all_trans)
% %%
% all_trans.summary.syl_labels_in_order
% all_trans.summary.divergence.matrix_of_fractions
% all_trans.summary.convergence.matrix_of_fractions

end