function [data,fs,spkind,spkamp,waves,tmstmps]=tetanaltw2(fname,TH,song_chan,tet_chans,wavesflag);
% [data,fs,spkind,spkamp]=tetanal(fname,TH,song_chan,tet_chans);
%

%SOME VALS
refrac=0.25e-3;
PKVAL=0;
SpkSrch=[-0.25e-3,1e-3];
SameSpkWin = 0.1e-3;
spklng=40;
edgsz=(spklng)/2;
sfct=(1e3/2^15)
waves=[];
tmstmps=[];
if (~exist('song_chan'))
	song_chan=1;
end

if (~exist('tet_chans'))
	tet_chans=[2:5];
end
[data,fs]=ReadCbinFile(fname);

for ii = 1:length(tet_chans)
    mn=mean(data(:,tet_chans(ii)));
    data(:,tet_chans(ii))=data(:,tet_chans(ii))-mn;
end

%TH=-2000;
NChan=length(tet_chans);
pk_win_sz=ceil(refrac*fs);
spkind=FindSpikeTimes(-sfct*data,fs,tet_chans,TH,PKVAL,SpkSrch,SameSpkWin);

spkind=spkind(find((spkind>(edgsz))&(spkind<length(data)-edgsz)));

if(wavesflag)
    for ii=1:length(spkind)
    
        %Waves(numpoints, number of electrode channels, size of spkind=number
        %of spikes);
        waves(ii,:,:)=-sfct*data(spkind(ii)+[(-edgsz+1):edgsz], tet_chans);
    end    
    tmstmps=uint32(spkind);
end
spkamp=zeros([size(spkind,1),NChan]);
for iChan = 1:NChan
	for ii=1:length(spkind)
		if (PKVAL==1)
			spkamp(ii,:) = -sfct*min(data(spkind(ii)+[-2:2],tet_chans));
		else
			spkamp(ii,:) = -sfct*max(data(spkind(ii)+[-2:2],tet_chans));
		end
	end
end
return;
