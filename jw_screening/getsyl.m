function syl1=getsyl(i,soundfile)
%get syl1 (a 10ms smoothed syl vector) using syl1=getsyl(i,soundfile)
%where i is the index into the labels, onsets, or offsets vectors in a .not.mat
%file that matches the soundfile (soundfile is wav only)

[song, Fs]=evsoundin('',soundfile,'w');
if ~exist([soundfile,'.not.mat'])
	disp([soundfile,'.not.mat not found.'])
	return
end

load([soundfile,'.not.mat']);

syl=song(floor(onsets(i)*Fs/1000):ceil(offsets(i)*Fs/1000));
syl1=smoother(syl,Fs,10);
