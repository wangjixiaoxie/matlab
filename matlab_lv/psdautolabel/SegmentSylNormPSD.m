function [PSDSylOut, SylOut, alllabels, Fs, gapdur]=SegmentSylNormPSD(filenames,filetype,nparts)
%modified from Hamish to have multiple psds, lower freq resolution, gap duration....

PSDSylOut=[];
syllengthout=[];
alllabels=[];
SylOut={};
zeropad=zeros(1,1050);
gapdur = [];


for i=1:length(filenames)
    
    fprintf('Working on %d of %d',i,length(filenames))
    notdata=load([filenames{i} '.not.mat']);
    notdata.onsets=notdata.onsets/1000;
    notdata.offsets=notdata.offsets/1000;
    dt=(1/notdata.Fs); %*1000; % in ms;
    
    alllabels=[alllabels notdata.labels];
    
    [song, Fs]=ReadAllSongFiles('./',filenames{i},filetype);
    song = bandpass(song,Fs,300,15000);
    if size(song,2)==1
        song=song';
    end;
    clear songdata;
    
    syllength=notdata.offsets-notdata.onsets
    gapdur_pre = [nan; notdata.onsets(2:end) - notdata.offsets(1:end-1)];
    gapdur_post = [ notdata.onsets(2:end) - notdata.offsets(1:end-1); nan];
    
    for j=1:length(notdata.onsets)
        nopadsyl = song((notdata.onsets(j)*notdata.Fs):(notdata.offsets(j)*notdata.Fs)-1);
        
        syl=[zeropad nopadsyl zeropad];
        syl=syl-mean(syl);
        nopadsyl = nopadsyl - mean(nopadsyl);
        
        SylOut{end+1}=syl;
        try
            
            
            partdur = floor(length(nopadsyl)/nparts);
            psdsyl = [];
            for k = 1:nparts
                thispart = [zeropad nopadsyl((k-1)*partdur+1:k*partdur) zeropad];
                psdsyl = [psdsyl pwelch(thispart,2048,2044,(200:8:8000),Fs)];  %hier war 2,8000 teste weniger fine
                
                
            end
        catch
            keyboard
        end
        
        
        PSDSylOut(end+1,:)=(psdsyl-mean(psdsyl))./std(psdsyl);
        
        
    end
end