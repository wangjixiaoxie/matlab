function lv_resegment_fixed(batchname,filetype,threshold)

fid = fopen(batchname,'r');
[~,numlines] = fscanf(fid,'%s');
frewind(fid);

for i=1:numlines
    fn{i}=fgetl(fid);
end

SylMin = 25;
SylMax = 250;
gapmin = 3;
gapmax = 150;
FileType = filetype;

for i=1:length(fn)
    
    fprintf('working on %d of %d \n',i,length(fn))
    
    min_int=gapmin;
    min_dur=SylMin;
    sm_win=2.5;
    
    [song, Fs]=ReadAllSongFiles('./',fn{i},FileType);
    %     song = bandpass(song,Fs,500,15000); %quatsch, wird in evsmooth
    %     gemacht 500-10k
    sm=evsmooth(song,Fs,0.01,512,0.8, sm_win); %
    
    
    [onsets,offsets]=SegmentNotes(sm,Fs,min_int,min_dur,threshold);
    onsets=onsets*1000; % in S for EVSONGANALY
    offsets=offsets*1000; % in S for EVSONGANALY
    syllength=(offsets)-(onsets); % in ms.
    
    shortkillidx=find(syllength < SylMin);
    longkillidx=find(syllength > SylMax);
    killidx=[shortkillidx; longkillidx];
    onsets(killidx)=[];
    offsets(killidx)=[];
    
    gaplength = [500; onsets(2:end)-offsets(1:end-1); 500];
    gapkillidx = find(gaplength(1:end-1)>gapmax & gaplength(2:end)>gapmax);
    onsets(gapkillidx)=[];
    offsets(gapkillidx)=[];
    
    if ~isempty(offsets)
    if offsets(end)>=length(sm)
        offsets(end) = [];
        onsets(end) = [];
    end
    end
    
    
    labels=[];
    
    for j=1:length(onsets);
        labels=[labels '-'];
    end;
    
    save([fn{i} '.not.mat'],'min_int','Fs','min_dur','sm_win','onsets','offsets','threshold','labels');
end