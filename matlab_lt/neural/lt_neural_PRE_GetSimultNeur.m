function SummaryStruct = lt_neural_PRE_GetSimultNeur(SummaryStruct)
%% lt 12/2017 - for each "experiment" get all sets of simultaneous neurons
% gets those with some temporal overlap. for each set record the song files with
% overlap.

%% run
numbirds = length(SummaryStruct.birds);
for i=1:numbirds
    
    ListOfExpt = unique({SummaryStruct.birds(i).neurons.exptID});
    
    exptcounter = 1;
    for exptthis = ListOfExpt
        
        % ==== which neurons are in this expt?
        NeurInThisExpt = find(strcmp({SummaryStruct.birds(i).neurons.exptID}, exptthis));
        
        % ==== FIND "SETS" OF NEURONS.
        
        % ----------- METHOD 1: same batch name
        
        
        % ---------------- METHOD 2: overlapping song files
        SongFilenamesByNeur = {SummaryStruct.birds(i).neurons(NeurInThisExpt).Filedatestr_unsorted};
        
        % ------ go thru each song file. for that song file get set of neurons.
        % finally get all unique sets of neurons and their corresponding
        % song files
        % keep the sets that have at least 2 neurons
        
        % ------------------ 1) get list of unique filenames
        SongFilenamesByNeur_ALL = {};
        for j=1:length(SongFilenamesByNeur)
            sfiles = SongFilenamesByNeur{j};
            
            SongFilenamesByNeur_ALL = [SongFilenamesByNeur_ALL sfiles];
        end
        SongFilenamesUnique = unique(SongFilenamesByNeur_ALL);
        
        % ----------------- 2)  go thru each song, get neurons that active
        NeuronsInEachSong = {};
        NeuronsInEachSongStr = {};
        for sfile = SongFilenamesUnique
            
            func = @(X) any(strcmp(X, sfile)); % for each neur, asks whether it ahs this file.
            goodneur = cellfun(func, SongFilenamesByNeur);
            
            NeuronsInEachSong = [NeuronsInEachSong goodneur];
            NeuronsInEachSongStr = [NeuronsInEachSongStr num2str(goodneur)];
        end
        
        % ------------------- 3)  get all sets
        [~, inds1, inds2] = unique(NeuronsInEachSongStr); % unique id
        
        % -- neurons for each set
        Sets_neurons = NeuronsInEachSong(inds1);
        % -- convert back to original neuron ID
        for setnum=1:length(Sets_neurons)
            Sets_neurons{setnum} = NeurInThisExpt(Sets_neurons{setnum});
        end
        
        
        % -- songfiles for each set
        Sets_songfiles = {};
        for setnum = 1:length(Sets_neurons)
            % -- for each set, find out corrresponding inds in original
            % file
            
            songsinset = SongFilenamesUnique(inds2==setnum);
            Sets_songfiles =[Sets_songfiles {songsinset}];
        end
        
        % -- REMOVE SETS WITH <2 NEURONS
        setToRemove = [];
        for setnum = 1:length(Sets_neurons)
            if length(Sets_neurons{setnum})<2
                
                setToRemove = [setToRemove setnum];
            end
        end
        Sets_neurons(setToRemove) = [];
        Sets_songfiles(setToRemove) = [];
        
        
        % ============================= SAVE INFORMATION
        SummaryStruct.birds(i).exptnum_pop(exptcounter).Sets_neurons = Sets_neurons;
        SummaryStruct.birds(i).exptnum_pop(exptcounter).Sets_songfiles = Sets_songfiles;
        SummaryStruct.birds(i).exptnum_pop(exptcounter).exptname = exptthis{1};
        exptcounter = exptcounter+1;
        
        
        % =========================== sanity check
        if (0)
            lt_figure; hold on;
            % -- for each neuron plot dots for times of songs
            for neurnum = NeurInThisExpt
                x = SummaryStruct.birds(i).neurons(neurnum).Filedatenum_unsorted;
                y = neurnum;
                
                plot(x, y, 'ok');
                lt_plot_text(x(end), max(y), ['neur' num2str(y)]);
            end
            
            
            for ss = 1:length(Sets_neurons)
                disp(['=== Set num ' num2str(ss)]);
                disp(['Neurons: ' num2str(Sets_neurons{ss})]);
                disp(['First file: ' Sets_songfiles{ss}{1}]);
                disp(['Last file: ' Sets_songfiles{ss}{end}]);
                line([datenum(Sets_songfiles{ss}{1}, 'yymmdd_HHMMSS') ...
                    datenum(Sets_songfiles{ss}{1}, 'yymmdd_HHMMSS')], ylim);
                line([datenum(Sets_songfiles{ss}{end}, 'yymmdd_HHMMSS') ...
                    datenum(Sets_songfiles{ss}{end}, 'yymmdd_HHMMSS')], ylim);
                lt_plot_text(datenum(Sets_songfiles{ss}{1}, 'yymmdd_HHMMSS'), ...
                    Sets_neurons{ss}(end), ['set ' num2str(ss)], 'b')
                lt_plot_text(datenum(Sets_songfiles{ss}{end}, 'yymmdd_HHMMSS'), ...
                    Sets_neurons{ss}(end), ['set ' num2str(ss)], 'b')
            end
            pause;
            close all
        end
    end
end
