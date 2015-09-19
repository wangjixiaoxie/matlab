%% To find out which trials had misses between hit/match and labeled 
%
%Written by DM Brady 07/2012

seq_or_pitch = input('Seq or pitch?  ','s');


if strcmpi(seq_or_pitch, 'seq')
    for i = 1:length(tpofinterest)
        hit_misses.(tpofinterest{i}) = find(syllables.(tpofinterest{i}).hit.values(:,1)...
            - syllables.(tpofinterest{i}).hit.values(:,2) <= -1);
        match_misses.(tpofinterest{i}) = find(syllables.(tpofinterest{i}).match.values(:,1)...
            - syllables.(tpofinterest{i}).match.values(:,2) <= -1);
        
        display(' ')
        display([tpofinterest{i} ' hit misses (song number): '])
        display(hit_misses.(tpofinterest{i}))
        
        display(' ')
        display([tpofinterest{i} ' match misses (song number): '])
        display(match_misses.(tpofinterest{i}))
    end
    
elseif strcmpi(seq_or_pitch, 'pitch')
        for i = 1:length(check_pitch.parameters.sofinterest)
        hit_misses.(check_pitch.parameters.sofinterest{i}) = find(syllables.(check_pitch.parameters.sofinterest{i}).hit.values(:,1)...
            - syllables.(check_pitch.parameters.sofinterest{i}).hit.values(:,2) <= -1);
        match_misses.(check_pitch.parameters.sofinterest{i}) = find(syllables.(check_pitch.parameters.sofinterest{i}).match.values(:,1)...
            - syllables.(check_pitch.parameters.sofinterest{i}).match.values(:,2) <= -1);
       
        
        display(' ')
        display([check_pitch.parameters.sofinterest{i} ' hit misses (song number): '])
        display(hit_misses.(check_pitch.parameters.sofinterest{i}))
        
        display(' ')
        display([check_pitch.parameters.sofinterest{i} ' match misses (song number): '])
        display(match_misses.(check_pitch.parameters.sofinterest{i}))
        end
end
