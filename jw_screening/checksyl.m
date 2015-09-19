function checksyl(soundfile)
% ( [labels,scores,lengths,peakpos,pkheight,htvar]= )   checksyl(soundfile)

[song, Fs]=evsoundin('',soundfile,'w');
load([soundfile,'.not.mat']);
for i=1:length(labels)
	index=i;
	syl1=song(floor(onsets(index)*Fs/1000):ceil(offsets(index)*Fs/1000));

	smoothed=smoother(syl1,Fs,18);
	%%%plot(smoothed);hold;
	%testing(smoothed)
	%return
	[syl1score, syl1mins, sylmaxs]= mpeaks(smoothed);
	syl1mins=transpose(syl1mins);
	%%%plot(syl1mins, 1000000,'ro')
	
	%output variables calc here-- first is scores
	scores(i)=syl1score;
	
	%then lengths
	lengths(i)=offsets(i)-onsets(i);
	
	%then peakpos
	peakvals=[];
	for m=1:length(sylmaxs)
		peakvals(m)=smoothed(sylmaxs(m));
	end
	peakinds=(find(peakvals==max(peakvals))); %thats indices into sylmaxs that represent the syl's highest vals (peakinds often a scalar)
	inds=0;
	for m=1:length(peakinds)
		inds=inds+sylmaxs(peakinds(m));
	end
	peakpos(i)=(inds/length(peakinds))/length(smoothed); % mean of indices of abs maxes, divided by syl length.. a scale 0 to 1.
	
	%now peakht and ht variability.
	%peakht is NOT the height of the peak, it is "how much higher is the peak compared to avg ht of syl."
	peakheight(i)=max(peakvals)/mean(smoothed);
	htvar(i)=std(smoothed);
	
	disp([labels(index),',',num2str(index),',',num2str(scores(i)),',',num2str(lengths(index)),',',num2str(peakpos(index)),',',num2str(peakheight(index)),',',num2str(htvar(index))])
	%disp(['The ', labels(index), '-', num2str(index), ' had ', num2str(syl1score), ' peaks.'])
end

