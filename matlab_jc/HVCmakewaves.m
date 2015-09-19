function HVCmakewaves(folders,indfolders)

for i=1:length(indfolders) %each folder
    index=indfolders(i);
    clearvars -except folders indfolders index
    currfolder=folders{index};
    cmd=['cd ' currfolder]
    eval(cmd);
    dirf('*.song.mat','batchsongs')
    ff=load_batchf('batchsongs');
    for ifn=1:length(ff) % for each song
        fnn=ff(ifn).name;
        load(fnn)
        % This takes the raw song
        dat=Song;
        SongN=Song-mean(Song);
        SongN=SongN/max(SongN*1.0001);
        wavwrite(SongN,32000,[fnn(1:end-9) '.wav'])
    end
end
