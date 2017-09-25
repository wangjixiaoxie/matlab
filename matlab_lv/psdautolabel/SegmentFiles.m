function thresheffective = SegmentFiles(fn,SylMin,SylMax,gapmin,gapmax,FileType)
%Lena modified this (from Hamish?) to have Otsu method
%and to have lower threshold 2SD over noise

for i=1:length(fn)
    
    fprintf('working on %d of %d \n',i,length(fn))
    
    min_int=gapmin;
    min_dur=SylMin;
    sm_win=2.5;
    
    [song, Fs]=ReadAllSongFiles('./',fn{i},FileType);
%     song = bandpass(song,Fs,500,15000); %quatsch, wird in evsmooth
%     gemacht 500-10k
    sm=evsmooth(song,Fs,0.01,512,0.8, sm_win); %

    
    %Otsu method / Dave
    imagesm = log(sm);
    minim = min(imagesm);
    maxim = max(imagesm);
    imla = (imagesm-minim)./(maxim-minim);
    [th, thresheffective(i)] = graythresh(imla);
    th = th.*(maxim-minim) + minim;
    th = exp(th);
    threshold = th;
    
    %Feelab method: second lower threshold
    lognoise = imagesm(imagesm<log(threshold));
    lowerthresh = mean(lognoise) + 2*std(lognoise);
    lowerthresh = exp(lowerthresh);
    
    if thresheffective(i)<0.7
        onsets=[];
        offsets=[];
    else
        [onsets,offsets]=SegmentNotes_LV(sm,Fs,min_int,min_dur,threshold,lowerthresh);
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
        threshold = lowerthresh;
    end
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
end;
