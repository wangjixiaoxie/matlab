function [ history_dependence ] = db_history_dependence_calculation( syllables_per_song, motifs, number_bootstraps)
%db_history_dependence_calculation Given a cell of syllables (each cell is a song), a set of motifs,
%and the number of bootstraps, this will calculate the history dependence and a 95% CI. 
%Scrambles the syllables_per_song to calculate the shuffled history
%dependence. Data is presented as a structure history_dependence.motif{1}
%history_dependence.motif{2}.
%  HD = abs( p(ab|ab) - p(ab|ac) )

%% sets bootstrap number to 10,000 if not specified
if nargin < 3
    number_bootstraps = 10000;
end

%% See how many history dependent transitions there are
first_motifs = motifs;

%% See if it is a convergent or divergent sequence
con_or_div = db_con_or_div(motifs);

%% Counts the number of a certain hd transition per song

%goes through each song
for i = 1:length(syllables_per_song)
    %goes through the first motif
    for j = 1:length(first_motifs)
        %goes through the second set of motifs
        for k = 1:length(motifs)
            %if the sequence is divergent, it looks for the number of a
            %certain transition ex: 'bbbc' (ab_to_ac = 1)
            if strcmpi(con_or_div,'div')
                raw_count.([first_motifs{j} '_to_' motifs{k}]){i} = length(strfind(syllables_per_song{i},[first_motifs{j}(end) motifs{k}(end)]));
            elseif strcmpi(con_or_div, 'con')
                raw_count.([first_motifs{j} '_to_' motifs{k}]){i} = length(strfind(syllables_per_song{i},[first_motifs{j}(1) motifs{k}(1)]));
            end
        end
    end
end

%% Calculates the relative frequency of a certain hd tranistion per song

%lists the name of hd transitions to make comparisons
hd_motifs = fieldnames(raw_count);

%finds the underscores in hd_motifs
for i = 1:length(hd_motifs)
    location_underscores{i} = strfind(hd_motifs{i},'_');
end

%finds which cells to compare (last motif in hd_motif is the same)
for i = 1:length(hd_motifs)
    %determines which comparisons for history dependence (last motifs are the
    %same) ((ab|ab) - (ab|ac))
    which_comparisons_hd{i} = find(~cellfun('isempty',strfind(hd_motifs,hd_motifs{i}(location_underscores{i}(2):end))));
    
    %determines which comparisons for relative frequency (first motifs are
    %the same) ((ab -> ac vs. ab -> ab)
    which_comparisons_rf{i} = find(~cellfun('isempty',strfind(hd_motifs,hd_motifs{i}(1:location_underscores{i}(2)))));
end

%calculates the relative frequency of a certain hd_transition
% (# ab|ab)./ (# ab|ab + # ab|ac)
for j = 1:length(syllables_per_song)
    for i = 1:length(hd_motifs)
        total_per_song = [];
        for k = 1:length(which_comparisons_rf{i})
            total_per_song = [total_per_song raw_count.(hd_motifs{which_comparisons_rf{i}(k)}){j}];
        end
        total_per_song = sum(total_per_song);
        
        rel_freq.(hd_motifs{i}){j} = raw_count.(hd_motifs{i}){j}./total_per_song;
        
    end
end
        
%gets rid of NaN songs and convert to an array       
for i = 1:length(hd_motifs)
    rel_freq.(hd_motifs{i}) = rel_freq.(hd_motifs{i})(~cellfun(@(V) any(isnan(V(:))),rel_freq.(hd_motifs{i})));
    rel_freq.(hd_motifs{i}) = cell2mat(rel_freq.(hd_motifs{i}));
end

%% Calculating History dependence using a bootstrap procedure
for i = 1:length(motifs)
    for j = 1:number_bootstraps
        for n = 1:length(hd_motifs)
            sampling = randi(length(rel_freq.(hd_motifs{n})),[1 length(rel_freq.(hd_motifs{n}))]);
            boot_run(j,n) = mean(rel_freq.(hd_motifs{n})(sampling));
        end
    end

    history_dependence.(motifs{i}) = abs(diff(boot_run(:,which_comparisons_hd{i})'));
        
end
end

