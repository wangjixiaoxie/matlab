function lv_givelabels(batchname,newclusterlabels,newclusnum)

fid = fopen(batchname,'r');
[~,numlines] = fscanf(fid,'%s');
frewind(fid);
for i=1:numlines
    fn{i}=fgetl(fid);
end

% newclusnumwithwn = nan(size(oldbasismtx,1),1);
% newclusnumwithwn(~wnidx) = newclusnum;
% newclusnumwithwn(wnidx) = max(newclusnum)+1;
% newclusterlabels = [newclusterlabels 'w'];

% labelsout=[newclusterlabels{newclusnumwithwn}]; %remove here if not
% removed white noise
labelsout=[newclusterlabels{newclusnum}];

     
     %not doing this anymore because I'm shortening batch.keep.lvrand.keep
     %to 3000 syllables. still useful to do if labelling new file;
% all2labels = [];
% for i = 1:numlines
%     load([fn{i} '.not.mat']);
%     all2labels = [all2labels labels];
% end
%      
% assert(length(all2labels)==length(labelsout))

for i=1:numlines
    load([fn{i} '.not.mat']);
    try
    labels=labelsout(1:length(labels));
    
    labelsout(1:length(labels))=[];
    catch
        i
        break
        
    end
    save([fn{i} '.not.mat'],'min_int','Fs','min_dur','sm_win','onsets','offsets','threshold','labels');
end

fprintf('done labelling \n')