function lt_neural_BatchSmth_Clean(DATAllSwitches)
%% ============================== NOTE DOWN WHICH ARE NOISY TRIALS
% ---- will plot raw data for each channel/motif, all trials. type the
% trials that are noise. will save a file noting down noise trials dor each
% channel/motif.
% ---- NOTE: the saved file will be over the original song files, so that this
% code does not have to be rerun for every dat struct

%%
MotifsAll = {DATAllSwitches.switch(1).motif.motif};
ChansAll = find(~cellfun('isempty', DATAllSwitches.switch(1).motif(1).batchinorder(1).DatAll));


%%
for cc = ChansAll
    for mm = 1:length(MotifsAll)
        
        
        numswitches = length(DATAllSwitches.switch);
        for i=1:numswitches
            numbatches = length(DATAllSwitches.switch(i).motif(mm).batchinorder);
            
            for bb=1:numbatches
                
                % ======= already have? if so then skip
                if isfield(DATAllSwitches.switch(i).motif(mm).batchinorder(bb), ...
                        'NoiseTrialsToIgnore')
                    if ~isempty(DATAllSwitches.switch(i).motif(mm).batchinorder(bb).NoiseTrialsToIgnore)
                        % if done, then will either have numbers or have
                        % "all trials good" flag.
                        disp('SKIPPING!! - already dione');
                    end
                end
                
                
                cond = DATAllSwitches.switch(i).motif(mm).batchinorder(bb).condition;
                datraw = DATAllSwitches.switch(i).motif(mm).batchinorder(bb).DatAllRaw{cc};
                t = DATAllSwitches.switch(i).motif(mm).batchinorder(bb).t;
                % ================ plot each trial
                numtrials = size(datraw,1);
                
                % ------ plot N songs at a time
                Nsongs = 10;
                
                nextNsongs = 1:Nsongs;
                trialsToRemove = [];
                
                while nextNsongs(end)<=numtrials
                    % keep going until reach the end of list of songs
                    
                    
                    % ------------------- plot the next N songs
                    figcount=1;
                    subplotrows=Nsongs;
                    subplotcols=1;
                    fignums_alreadyused=[];
                    hfigs=[];
                    
                    for j=nextNsongs
                        
                        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                        title(['[ch' num2str(cc) ']-' MotifsAll{mm} ', trial num ' num2str(j)]);
                        plot(t, datraw(j,:));
                        axis tight;
                        ylim([-200 200]);
                    end
                    
                    tt = input('Which trials to remove (e.g. 1 2 3 or 5:10)? ', 's');
                    tt = str2num(tt);
                    disp(['(to remove) added trials: ' num2str(tt)]);
                    trialsToRemove = [trialsToRemove tt];
                    
                    
                    % ------ get list of next N songs
                    nextNsongs = nextNsongs+10;
                    close all;
                    
                end
                disp(['---- TRIALS TO REMOVE, for ch' num2str(cc) '-' MotifsAll{mm} ...
                    '-sw' num2str(i) '-' cond ' = ' num2str(trialsToRemove)]);
                
                if isempty(trialsToRemove)
                    trialsToRemove = 'all trials good';
                end
                
                %% === output
                % ########################## 1) update compiled structure 
                % ====== then update
                DATAllSwitches.switch(i).motif(mm).batchinorder(bb).NoiseTrialsToIgnore = trialsToRemove;
                
                % ===== save DATAllSwitches again
                disp(['SAVED to ' DATAllSwitches.savename]);
                save(DATAllSwitches.savename, 'DATAllSwitches')
                
                % ########################### 2) save an array linked to
                % original song file
               
                
            end
            
        end
        
    end
end

