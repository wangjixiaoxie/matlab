function [SummaryStruct NeurToRemove] = lt_neural_v2_DIAGN_RemoveOlap(SummaryStruct)
%% TO DO

% 1) instead of removing entire neuron, just remove data that is overlaping
% (e.g from SongDat)

%%
% NeurToRemove cel array, each cell one bird, and within cell array of neurons to remove


%%


numbirds = length(SummaryStruct.birds);
NeurToRemove = {};

for i=1:numbirds
   
    birdname = SummaryStruct.birds(i).birdname;
    
    % === same units if have same channel, clustnum, and song filenames
    numneurons = length(SummaryStruct.birds(i).neurons);
    
    neurtoremove = [];
    
    disp(birdname);
    
    for ii=1:numneurons
        
        
        for iii=ii+1:numneurons
            
            chan1 = SummaryStruct.birds(i).neurons(ii).channel;
            clust1 = SummaryStruct.birds(i).neurons(ii).clustnum;
            elecdepth1 = SummaryStruct.birds(i).neurons(ii).electrode_depth;
            
            chan2 = SummaryStruct.birds(i).neurons(iii).channel;
            clust2 = SummaryStruct.birds(i).neurons(iii).clustnum;
            elecdepth2 = SummaryStruct.birds(i).neurons(iii).electrode_depth;
            
            if chan1~=chan2 | clust1~=clust2 | elecdepth1~=elecdepth2
                continue
            end
            
            % ============== CHECK IF OVERLAPPING SONG FILENAMES
            filedates1 = SummaryStruct.birds(i).neurons(ii).Filedatestr_unsorted;
            filedates2 = SummaryStruct.birds(i).neurons(iii).Filedatestr_unsorted;
            
            cd(SummaryStruct.birds(i).neurons(ii).dirname);
            tmp1 = load('times_data.mat');
            cd(SummaryStruct.birds(i).neurons(iii).dirname);
            tmp2 = load('times_data.mat');
            
            if ~strcmp(tmp1.par.detection, tmp2.par.detection)
                % neg/pos spike detection, if different then keep both.
                continue
            end
            
            
            if ~isempty(intersect(filedates1, filedates2))
                % ====== THEN PROBLEM --- REMOVE ONE OF THE NEURONS
                % REMOVE THE ONE WITH LESS DATA (I.E. NUM SONGS)
                
                disp('pair: ');
                disp(SummaryStruct.birds(i).neurons(ii).dirname);
                disp(SummaryStruct.birds(i).neurons(iii).dirname);
                
                numsongs1 = length(filedates1);
                numsongs2 = length(filedates2);
                
                if numsongs1>=numsongs2
                    neurtoremove = [neurtoremove iii];
                else
                    neurtoremove = [neurtoremove ii];
                end
                
            end

            
            
%             % ================== OLD VERSION
% %             tic
% %             cd(SummaryStruct.birds(i).neurons(ii).dirname);
% %             mtdat1 = load('MetaDat.mat');
% %             toc
%             tic
%             [~, NeurDat1, ~] = lt_neural_ExtractDat2(SummaryStruct, i, ii);
%             toc
%             [~, NeurDat2, ~] = lt_neural_ExtractDat2(SummaryStruct, i, iii);
%             
%             if ~isempty(intersect({NeurDat1.metaDat.filename}, {NeurDat2.metaDat.filename}))
%                 % ====== THEN PROBLEM --- REMOVE ONE OF THE NEURONS
%                 % REMOVE THE ONE WITH LESS DATA (I.E. NUM SONGS)
%                 
%                 numsongs1 = length(NeurDat1.metaDat);
%                 numsongs2 = length(NeurDat2.metaDat);
%                 
%                 if numsongs1>=numsongs2
%                     neurtoremove = [neurtoremove iii];
%                 else
%                     neurtoremove = [neurtoremove ii];
%                 end
%                 
%             end
        end
    end
    
    if ~isempty(neurtoremove)
       disp(['=== REMOVING NEURONS [duplicate of another] ' num2str(neurtoremove) ]); 
       SummaryStruct.birds(i).neurons(neurtoremove) = [];
    end
    
    NeurToRemove{i} = neurtoremove;
    
end