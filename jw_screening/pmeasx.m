function [peakinds]=pmeasx(filenm)

%830 to 937 w83pk50_010411_065602.wav

[sounddata, Fs] = evsoundin('',filenm,'w');
sounddata = highpass(sounddata,200,Fs);
fftbin=512;
[s,f,t,p]=spectrogram(sounddata,fftbin,300,[],Fs); %sounddata,512,300,[],44100 means 212samples per segment.
peakinds=zeros(size(p,2),floor(.33*size(p,1)));
pitchs=zeros(size(p,2),1);
for loc=830:937
spec=log(smooth(p(1:floor(.33*size(p,1)),loc),3));
highs=spec>(max(spec)-min(spec))*.6+min(spec);
highinds=find(highs);
slopes=diff(spec);
ups=slopes>0;
minmax=conv([1 -1],single(ups));
count=0;
for i=1:length(highinds)
	if minmax(highinds(i))==-1
		count=count+1;
		%peakinds(loc-905,:)
		peakinds(loc,count)=highinds(i);
	end
end
peakinds(loc,1:count)
diff(peakinds(loc,1:count))

pitchs(loc)=mean(diff(peakinds(loc,1:count)));
end
%figure;plot(spec);hold on;plot(minmax)
pitchs=pitchs.*((Fs/2)/(fftbin/2));  %get freqs in Hz instead of bin#
%figure;plot(pitchs)

return

