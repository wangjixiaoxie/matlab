function [trans_per_song, time_per_song, con_or_div, boot_results, entropy_results, syl_pool] =...
    db_transition_probability_for_sequence( batchfile, motifs, boot_yes_or_no, number_bootstraps )
%db_transition_probability_per_day Given a batch file or single song (cbin or
%cbin.not.mat) and a set of motifs (needs to be a cell, ex: {'ab' 'ac'}), 
%it will give the following:
% 1. a list of the transition syllables per song (trans_per_song)
% 2. a list of when the transition syllables were sung per song (time_per_song)
% 3. whether you input a convergent or divergent sequence
% 4. if 'boot_yes_or_no is not specified or set to 'y', a bootstrap
% structure of each motif's transition probability. The number of
% iterations of the bootstrap can be set by number_bootstraps (default is
% 10,000). 
% 5. Transition entropy(-Sigma p log2 p)
% 6. Pooled syllables for the day.

%sets number_bootstraps if not specified and will run bootstrap by default
if nargin < 3
    boot_yes_or_no = 'y';
    number_bootstraps = 10000;
elseif nargin < 4 && strcmpi(boot_yes_or_no, 'y') == 1
    number_bootstraps = 10000;
end

%gets all the labels for the day in batchfile
try
    [all_syllables, time_syl] = db_get_labels(batchfile);
catch err
    try
        %if there was an error running db_get_labels as if it were a batch
        %file, it runs it as if it were a single song
        [all_syllables, time_syl] = db_get_labels(batchfile,'y');
    catch err
        display('Something is wrong with you batchfile')
        return
    end
end

%checks to see if motifs are good, and if running a convergent or divergent
%transition probability calculation
if length(motifs) < 2
    display('Not enough motifs')
    return
else
    [trans_per_song, con_or_div, boot_results, entropy_results, syl_pool, time_per_song] = ...
        db_transition_probability_calculation( all_syllables, motifs, boot_yes_or_no, number_bootstraps, time_syl );
end
        
        
        
    
    

end

