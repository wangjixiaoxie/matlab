function [  ] = db_seq_func_day_save_v2(batchfile, computer, birdname, date, data_type, motif, phrase, varargin)
%db_seq_func_day_save Will perform the function db_motif_function and save
%the data to the appropriate folder
%   batchfile, computer, birdname, date, data_type, motif, phrase

if ischar(motif)
    motif = {motif};
end

if data_type == 1;
    data_type = {'dur' 'amp' 'frac' 'ent'};
end

if ischar(data_type)
    data_type = {data_type};
end


%if computer, batchfile, and birdname are empty it will fill in the apporpriate names
%(assuming that you are in a song folder for that day)
if isempty(batchfile) == 1
    batchfile = 'batch.catch.keep';
end

if isempty(computer) == 1
    temp = pwd;
    cd ..
    computer = pwd;
    cd(temp)
end

if isempty(birdname) == 1
    temp = pwd;
    slashes = strfind(temp,'/');
    birdname = temp(slashes(end-1)+1:slashes(end)-1);
end

for k = 1:length(data_type)
    for j = 1:max(size(motif))
        
        %% For motif timing/duration
        if ~isempty(strfind(data_type{k}, 'dur')) == 1;
            %checks to see if there is a folder to put timing data
            if exist(['/' computer '/all_days_gap_timing_' phrase],'dir') == 7
            else
                mkdir(['/' computer '/all_days_gap_timing_' phrase])
            end
            
            %runs db_motif_timing
            [data] = db_motif_timing(batchfile,motif{j});
            
            %saves data and time in to timing folder
            save(['/' computer '/all_days_gap_timing_' phrase '/' birdname '_' date '_gapduration_' motif{j} '.mat'], 'data', 'time')
            
            %% For amplitude and amplitude ratio
        elseif ~isempty(strfind(data_type{k},'amp')) == 1;
            %checks to see if there is a folder for amplitude data
            if exist(['/' computer '/all_days_amplitude_' phrase],'dir') == 7
            else
                mkdir(['/' computer '/all_days_amplitude_' phrase])
            end
            
            %if varargin is not empty, it will be the base syllable, otherwise
            %it is the first syllable of the motif.
            if ~isempty(varargin) == 1
                base_syl = varargin{1};
            elseif length(motif{j}) == 1
                base_syl = motif{j};
            else
                base_syl = motif{j}(1);
            end
            
            %runs through and measures the amplitude for each syllable in the
            %motif.
            for i = 1:length(motif{j})
                [amps, times, ratios] = db_motif_amplitude(batchfile,motif{j},base_syl);
                time = times(:,i);
                data = amps(:,i);
                
                %saves the amp and time for the current syllable
                save(['/' computer '/all_days_amplitude_' phrase '/'...
                    birdname '_' date '_amplitude_' motif{j}(i) '_from_' motif{j} '.mat'], 'data', 'time')
                %             display(['/' computer '/dbrady/' birdname '/all_days_amplitude/'...
                %                 birdname '_' date '_amplitude_' motif{j}(i) '_from_' motif{j} '.mat'])
                
                %if the motif is longer than one syllable, a ratio is
                %calculated using current syllable and base syllable
                if motif{j}(i) == base_syl
                else
                    data = ratios(:,i);
                    save(['/' computer '/all_days_amplitude_' phrase '/'...
                        birdname '_' date '_ratio_amp_' motif{j}(i) '_over_' base_syl '_from_' motif{j} '.mat'], 'data', 'time')
                    %                 display(['/' computer '/dbrady/' birdname '/all_days_amplitude/'...
                    %                     birdname '_' date '_ratio_amp_' motif{j}(i) '_over_' base_syl '_from_' motif{j} '.mat'])
                end
            end
            
            %% For calculating fraction of motif sung (# reps motif/ total number of syllables)
        elseif ~isempty(strfind(data_type{k},'frac')) == 1;
            
            %checks to see if a folder exists for the fraction data
            if exist(['/' computer '/all_days_fraction_' phrase],'dir') == 7
            else
                mkdir(['/' computer '/all_days_fraction_' phrase])
            end
            
            %calculates the fraction of syllable sung by bootstrap method
            [percentiles bootstrap frac_per_song] = db_fraction_syl(batchfile,motif{j});
            
            %saves fraction data
            save(['/' computer '/all_days_fraction_' phrase '/' birdname '_' date '_fraction_' motif{j} '.mat'],...
                'percentiles', 'bootstrap', 'frac_per_song')
            
            %% For calculating entropy of a motif (from onset of motif to offset)
        elseif ~isempty(strfind(data_type{k},'ent')) == 1;
            
            %checks to see if a folder exists for entropy data
            if exist(['/' computer '/all_days_entropy_' phrase],'dir') == 7
            else
                mkdir(['/' computer '/all_days_entropy_' phrase])
            end
            
            %%calculates entropy of the motif
            [data, time] = db_motif_entropy(batchfile,motif{j});
            
            %saves entropy data
            save(['/' computer '/all_days_entropy_' phrase '/' birdname '_' date '_entropy_' motif{j} '.mat'],...
                'data', 'time')
            
        end
    end
end


end

