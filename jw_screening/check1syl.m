function [labs,scores,lengths,peakpos,pkheight,htvar]=check1syl(songf,syl,s1,s2)
% ( [labs,scores,lengths,peakpos,pkheight,htvar]=check1syl(songf,syl,s1,s2)
%s1 and s2 are the 2 smoothing windows in ms. If == then only one smoother() is run.
%syl is an integer-- an index into the onsets vector in songf.not.mat (songf is a wav file).

index=1;	%holdover from checksyl2, a multiple syl version.

[song, Fs]=evsoundin('',songf,'w');	
load([songf '.not.mat']);
	if ischar(labels(syl))
	labs(index)=labels(syl);
	else
	labs(index)=0;
	end
	syl1=song(floor(onsets(syl)*Fs/1000):ceil(offsets(syl)*Fs/1000));
	length(syl1)

	smoothed=smoother(syl1,Fs,s1);	%this "hypersmoothing" is for the mpeaks analysis, and for the pkheight. Pkheight because
					%identical syls so often have their biggest peak of varying heights.
					%now peakht and ht variability.
	

	[syl1score, syl1mins, sylmaxs]= mpeaks_del(smoothed);
	plot(smoothed);hold;
	%syl1mins=transpose(syl1mins);
	
	%output variables calc here-- first is scores
	scores(index)=syl1score;
	%peakht is NOT the height of the peak, it is "how much higher is the peak compared to avg ht of syl."
	%for m=1:length(sylmaxs)
	%	peakvals(m)=smoothed(sylmaxs(m));
	%end
	%tried to do pkheight in a multi-peak way
	%maxpk=max(peakvals); %maxpk for pkheight --- different from maxpk for peakpos which will be the actual max
	%maxpknum=length(find(peakvals==maxpk));
	
	%pkheight(index)=maxpk/mean(smoothed);
	%htvar(index)=std(smoothed);

	
	%then lengths
	lengths(index)=offsets(syl)-onsets(syl);
	
	%%%for peakpos, pkht, ht var, will use a new smoothing job-- less smoothing.
	
	%%%BRING THESE 2 LINES BACK FOR 2 SMOOTHING STEPS
	if s1~=s2
	smoothed=smoother(syl1,Fs,s2);
	[syl1score, syl1mins, sylmaxs]= mpeaks(smoothed);
	semilogy(smoothed,'r');
	end
	%now peakpos. will avg the indices of the peaks in the top 15% of max syl height.
	peakvals=[];
	for m=1:length(sylmaxs)
		peakvals(m)=smoothed(sylmaxs(m));	%pkvals: the y-values of smoothed... index pos in sylmaxs is meaningless.
		logpkvals=log10(peakvals);
	end
	%maxpk=max(peakvals);
	%peakinds=(find(peakvals>(.7*maxpk))); %thats indices into sylmaxs that represent the syl's highest vals (peakinds often a scalar)
	maxpk=max(logpkvals);			%log version
	peakinds=(find(logpkvals>(.94*maxpk)));	%log version...... .93 is good... .94 better..
	sylmaxs(peakinds)
	inds=0;
	for m=1:length(peakinds)
		inds=inds+sylmaxs(peakinds(m));	%gets out the indices into smoothed. sums for averaging.
	end
	peakpos(index)=(inds/length(peakinds))/length(smoothed); % mean of indices of abs maxes, divided by syl length.. a scale 0 to 1.
	
	%pkheight here
	%maxpk=max(peakvals);
	%new: averaging the top 2 peaks instead of just maxpk.. unless second max is too low.
	secmax=max(logpkvals(find(logpkvals<maxpk)));
	if isempty(secmax)|secmax<.9*maxpk
	secmax=maxpk;
	end
	%maxpk
	%secmax
	%secmax=max(peakvals(find(peakvals<max(peakvals))))
	pkheight(index)=((maxpk+secmax)/2)/mean(log10(smoothed));
	htvar(index)=std(log10(smoothed));
	disp(['lab ',labels(syl),',','score ',num2str(scores(index)),',','len ',num2str(lengths(index)),',','pkpos ',num2str(peakpos(index)),',','pkht ',num2str(pkheight(index)),',','htvar ',num2str(htvar(index))])
	
	%no good: disp(['The ', labels(index), '-', num2str(index), ' had ', num2str(syl1score), ' peaks.'])
	
