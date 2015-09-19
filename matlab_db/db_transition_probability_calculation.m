function [trans_per_song, con_or_div, boot_results, entropy_results, syl_pool, time_per_song] =...
    db_transition_probability_calculation( all_syllables, motifs, boot_yes_or_no, number_bootstraps, time_syl )
%db_transition_probability_calculation Given a cell of syllables and a set 
%of motifs (needs to be a cell, ex: {'ab' 'ac'}), it will give the following:
% 1. a list of the transition syllables per song (trans_per_song)
% 2. whether you input a convergent or divergent sequence
% 3. if 'boot_yes_or_no is not specified or set to 'y', a bootstrap
% structure of each motif's transition probability. The number of
% iterations of the bootstrap can be set by number_bootstraps (default is
% 10,000). 
% 4. Transition entropy(-Sigma p log2 p)
% 5. Pooled syllables for the day.
% 6. If time of syllables is input, it will give you time of syllables of
% interest

%sets number_bootstraps if not specified and will run bootstrap by default
if nargin < 3
    boot_yes_or_no = 'y';
    number_bootstraps = 10000;
    time_syl = '';
elseif nargin < 4 && strcmpi(boot_yes_or_no, 'y') == 1
    number_bootstraps = 10000;
    time_syl = '';
elseif nargin < 4 && strcmpi(boot_yes_or_no, 'n') == 1
    time_syl = '';
end

for i = 1:length(motifs)
    %specifies the length of the motif excluding the first syllable
    length_motif{i} = length(motifs{i})-1;
end

%checks if motifs form a convergent or divergent sequence
[con_or_div modified_motifs] = db_con_or_div(motifs);

%loop goes through each song to get the transition syllables of
%interest
for k = 1:length(all_syllables)
    %creates a variable that will be used to eliminate all syllables
    %that we do not care about
    interest_syllables = [];
    
    %goes through each type of motif
    for ii = 1:length(motifs)
        if strcmpi(con_or_div, 'div') == 1
            %for divergent sequences, finds when the divergent sequence
            %occured in each song
            interest_syllables = [interest_syllables strfind(all_syllables{k},motifs{ii})+length_motif{i}];
        elseif strcmpi(con_or_div, 'con') == 1
            %finds when the convergent sequence occured in each song
            interest_syllables = [interest_syllables strfind(all_syllables{k},motifs{ii})];
        end
    end
    %sorts the interest_syllables so that the order is correct
    interest_syllables = sort(interest_syllables);
    
    %creates a string per song of the transitions of interest
    trans_per_song{k} = all_syllables{k}(interest_syllables);
    
    %says the time (in serial date time) when the transition occured
    if ~isempty(time_syl) == 1
        time_per_song{k} = time_syl{k}(interest_syllables);
    end
    
    clear interest_syllables
end

%gets rid of any song that did not contain the transitions of interest
trans_per_song = trans_per_song(~cellfun('isempty',trans_per_song));
if ~isempty(time_syl) == 1
    time_per_song = time_per_song(~cellfun('isempty',trans_per_song));
end

%creates a structure for bootstrap data if option is selected
if strcmpi(boot_yes_or_no,'y') == 1
    
    %creates a pool of all transistions for the day
    syl_pool = [];
    for kk = 1:length(trans_per_song)
        syl_pool = [syl_pool trans_per_song{kk}];
    end
    
    %bootstrap procedure for calculating transition probabilities
    for iii = 1:number_bootstraps
        %creates a vector that indicates which syllables will be
        %sampled during this iteration
        sampling = randi(length(syl_pool),[1 length(syl_pool)]);
        
        %creates a list of syllables by sampling from syl_pool with
        %sampling variable above
        boot_pool = syl_pool(sampling);
        
        for jjj = 1:length(motifs)
            %goes through the different motifs and calculates how many
            %times it occured for this bootstrap iteration
            boot_results.(motifs{jjj})(iii) = length(strfind(boot_pool,modified_motifs{jjj}))./length(boot_pool);
        end
        
        %calculates transition entropy for this bootstrap iteration
        %(entropy = Sigma -p * log2 p)
        entropy_iteration = [];
        for nnn = 1:length(motifs)
            entropy_iteration = [entropy_iteration boot_results.(motifs{nnn})(iii).*log2(boot_results.(motifs{nnn})(iii))];
        end
        
        entropy_results(iii) = sum(entropy_iteration)*-1;
        clear entropy_iteration
        
    end
    
    %replaces all NaN in entropy_results (when p = 0) to 0 (definition)
    entropy_results(isnan(entropy_results)) = 0;
    
    for mmm = 1:length(motifs)
        %transposes bootstrap results so it is easier to read
        boot_results.(motifs{mmm}) = boot_results.(motifs{mmm})';
    end
    
    entropy_results = entropy_results';
    
    %now you have a structure of bootstrap results. the median is the
    %transition probability for the day, and you can add 95% or 99% CI
end


        
    
    

end

