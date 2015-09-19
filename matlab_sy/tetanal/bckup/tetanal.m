function [data,fs,spkind,spkamp]=tetanal(fname,TH,song_chan,tet_chans);
%

%SOME VALS
refrac=0.25e-3;
PKVAL=1;
SpkSrch=[-0.25e-3,1e-3];
SameSpkWin = 0.1e-3;

if (~exist('song_chan'))
	song_chan=1;
end

if (~exist('tet_chans'))
	tet_chan=[2:5];
end
[data,fs]=ReadCbinFile(fname);

%TH=-2000;
NChan=length(tet_chans);
pk_win_sz=ceil(refrac*fs);
spkind=FindSpikeTimes(data,fs,tet_chans,TH,PKVAL,SpkSrch,SameSpkWin);

spkamp=zeros([size(spkind,1),NChan]);
for iChan = 1:NChan
	for ii=1:length(spkind)
		if (PKVAL==1)
			spkamp(ii,:) = min(data(spkind(ii)+[-2:2],tet_chans));
		else
			spkamp(ii,:) = max(data(spkind(ii)+[-2:2],tet_chans));
		end
	end
end
return;
