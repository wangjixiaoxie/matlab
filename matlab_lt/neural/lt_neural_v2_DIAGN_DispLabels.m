function lt_neural_v2_DIAGN_DispLabels(SummaryStruct, stoponbird)


%% lt 9/25/17 - displays all songs and labels

numbirds = length(SummaryStruct.birds);

for i=1:numbirds
    
    
    % ########################## display things about this bird
    disp(' ');
    disp([' ############################################## ' ...
        SummaryStruct.birds(i).birdname ' ################# ']);
        
   % === go thru all neurons
    numneurons = length(SummaryStruct.birds(i).neurons);
    for ii=1:numneurons
   
        cd(SummaryStruct.birds(i).neurons(ii).dirname);
        batchfile = SummaryStruct.birds(i).neurons(ii).batchfilename;
        eval(['!cp ' batchfile ' ..']); % copy batch file one dir up. [as this is most accurate batch file]
        cd .. % labeled files are up one.
        
        disp(' ');
        disp(' +++++++++++++++++++ ')
        
        location = SummaryStruct.birds(i).neurons(ii).NOTE_Location;
        disp(['neuron ' num2str(ii) ' out of ' num2str(numneurons) ' [' ...
            location ']']);
        disp(SummaryStruct.birds(i).neurons(ii).dirname)
        
        lt_batch_disp_all_labels(batchfile);
        
        
    end
    
    if stoponbird==1
       
        input(' !!!!!!!!!! PAUSED, press enter !!');
        close all
    end
end
