function [labs,scores,lengths,peakpos,pkheight,htvar]=checksyl2(batchfile,s1,s2)
% ( [labs,scores,lengths,peakpos,pkheight,htvar]=checksyl2(batchfile,s1,s2)
%s1 and s2 are the 2 smoothing windows in ms. If == then only one smoother() is run.
index=0;
fid=fopen(batchfile);
while 1
soundfile=fgetl(fid);
if ~ischar(soundfile), break, end
	
if soundfile(length(soundfile)-7:length(soundfile))=='.not.mat'
	if ~exist([soundfile])
		disp([soundfile,' not found.'])
		fclose(fid);
		return
	end
	load([soundfile]);
	disp(['computing ' soundfile])
	soundfile=soundfile(1:length(soundfile)-8);
	else
	if ~exist([soundfile,'.not.mat'])
		disp([soundfile,'.not.mat not found.'])
		fclose(fid);
		return
	end
	load([soundfile,'.not.mat']);disp(['computing ' soundfile '.not.mat'])
end
[song, Fs]=evsoundin('',soundfile,'w');	


for i=1:length(onsets)
	index=index+1;
	if ischar(labels(i))
	labs(index)=labels(i);
	else
	labs(index)=0;
	end
	syl1=song(floor(onsets(i)*Fs/1000):ceil(offsets(i)*Fs/1000));
%disp(['now on ',soundfile,' ',labels(i),'-',num2str(i)])
	smoothed=smoother(syl1,Fs,s1);	%this "hypersmoothing" is for the mpeaks analysis, and for the pkheight. Pkheight because
					%identical syls so often have their biggest peak of varying heights.
					%now peakht and ht variability.
	

	[syl1score, syl1mins, sylmaxs]= mpeaks_del(smoothed);
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
	lengths(index)=offsets(i)-onsets(i);
	
	%%%for peakpos, pkht, ht var, will use a new smoothing job-- less smoothing.
	
	%%%BRING THESE 2 LINES BACK FOR 2 SMOOTHING STEPS
	if s1~=s2
	smoothed=smoother(syl1,Fs,s2);
	[syl1score, syl1mins, sylmaxs]= mpeaks(smoothed);
	end
	%now peakpos. will avg the indices of the peaks in the top 6% of max syl log-height.
	peakvals=[];
	for m=1:length(sylmaxs)
		peakvals(m)=smoothed(sylmaxs(m));	%pkvals: the y-values of smoothed... index pos in sylmaxs is meaningless.
		logpkvals=log10(peakvals);
	end
	%maxpk=max(peakvals);
	%peakinds=(find(peakvals>(.7*maxpk))); %thats indices into sylmaxs that represent the syl's highest vals (peakinds often a scalar)
	maxpk=max(logpkvals);			%log version
	peakinds=(find(logpkvals>(.94*maxpk)));	%log version
	inds=0;
	for m=1:length(peakinds)
		inds=inds+sylmaxs(peakinds(m));%gets out the indices into smoothed. sums for averaging.
	end
	peakpos(index)=(inds/length(peakinds))/length(smoothed); % mean of indices of abs maxes, divided by syl length.. a scale 0 to 1.
	
	%pkheight here
	%maxpk=max(peakvals);
	%new: averaging the top 2 peaks instead of just maxpk.. unless second max is too low.
	secmax=max(logpkvals(find(logpkvals<maxpk)));
	if isempty(secmax)|secmax<.9*maxpk
	secmax=maxpk;
	%disp(['secmax set to maxpk for',soundfile,' ',labels(i),'-',num2str(i)])
	end
	if length(((maxpk+secmax)/2)/mean(log10(smoothed)))~=1
	disp(['WARNING pkheight computation in ',soundfile,' ',labels(i),'-',num2str(i)])
	maxpk
	secmax
	end
	pkheight(index)=((maxpk+secmax)/2)/mean(log10(smoothed));
	htvar(index)=std(log10(smoothed));
%can del if r38g19 c's aren't bifricated anymore in the scores dimension..
%	if labels(i)=='c'&scores(index)<.8
%		disp(['Found little c in ',soundfile,' ',labels(i),'-',num2str(i),' =',num2str(scores(index))])
%		if labels(i-1)=='c'
%		disp('Not the first in train.')
%		else
%		disp('first in train.')
%		end
%	end
	
	%disp([labels(i),',',num2str(i),',',num2str(scores(index)),',',num2str(lengths(index)),',',num2str(peakpos(index)),',',num2str(pkheight(index)),',',num2str(htvar(index))])
	
	%no good: disp(['The ', labels(index), '-', num2str(index), ' had ', num2str(syl1score), ' peaks.'])
end
end %for while loop - (batchfile fgetl)
fclose(fid);
