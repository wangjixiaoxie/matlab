function SummaryStruct = lt_neural_DISP_AllPopUnits(SummaryStruct)
%% summary plot of all population units

%% run
numbirds = length(SummaryStruct.birds);
for i=1:numbirds
    birdname = SummaryStruct.birds(i).birdname;
    numexpt = length(SummaryStruct.birds(i).exptnum_pop);
    
    for ii = 1:numexpt
        exptname = SummaryStruct.birds(i).exptnum_pop(ii).exptname;
        
        lt_figure; hold on;
        title([birdname '-' exptname]);
        
        NeurInThisExpt = find(strcmp({SummaryStruct.birds(i).neurons.exptID}, exptname));
        
                % -- for each neuron plot dots for times of songs
                for neurnum = NeurInThisExpt
                    x = SummaryStruct.birds(i).neurons(neurnum).Filedatenum_unsorted;
                    y = neurnum;
        
                    plot(x, y, 'ok');
                    lt_plot_text(x(end), max(y), ['neur' num2str(y)]);
                end
        
        % ================== figure out sets
        Sets_neurons = SummaryStruct.birds(i).exptnum_pop(ii).Sets_neurons;
        Sets_songfiles = SummaryStruct.birds(i).exptnum_pop(ii).Sets_songfiles;
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
%             lt_plot_text(datenum(Sets_songfiles{ss}{end}, 'yymmdd_HHMMSS'), ...
%                 Sets_neurons{ss}(end), ['set ' num2str(ss)], 'b')
        end
    end
end
