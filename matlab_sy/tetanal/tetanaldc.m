function [data,fs,spkind,spkamp]=tetanal(fname,TH,song_chan,tet_chans);
% [data,fs,spkind,spkamp]=tetanal(fname,TH,song_chan,tet_chans);
% 10/05 changed by TW to take peaktopeak amplitude

%SOME VALS
refrac=0.25e-3;
PKVAL=1;
SpkSrch=[-0.25e-3,1e-3];
SameSpkWin = 0.1e-3;

if (~exist('song_chan'))
	song_chan=1;
end

if (~exist('tet_chans'))
	tet_chans=[2:5];
end
[data,fs]=ReadCbinFile(fname);


for kk=1:4
    ind=tet_chans(kk);
    test=mean(data,1);
    data(:,ind)=data(:,ind)-test(ind);
end

%TH=-2000;
NChan=length(tet_chans);
pk_win_sz=ceil(refrac*fs);
spkind=FindSpikeTimes(data,fs,tet_chans,TH,PKVAL,SpkSrch,SameSpkWin);

spkamp=zeros([size(spkind,1),NChan]);
for iChan = 1:NChan
    for ii=1:length(spkind)
        if ((spkind(ii)>2)&(spkind(ii)<length(data)-2))
            if (PKVAL==1)
                spkamp(ii,:) = -min(data(spkind(ii)+[-2:2],tet_chans));
            else
                spkamp(ii,:) = max(data(spkind(ii)+[-2:2],tet_chans));
            end
        end
    end
end
return;
