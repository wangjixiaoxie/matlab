function [scores]=sylcomp(syl1,batchfile)
%[scores]=sylcomp(syl1,batchfile)
%syl1 is a vector--smoothed (10ms) syllable to compare syls with
%get syl1 using syl1=getsyl(i,soundfile)
%
%[scores] has 4 columns: song index, sylindex, label (use char(scores(1,3))).. is 0 if no label, comparisonresultscore.

%make truncated syl1-- will be used for the score calc.
[row1,col1]=size(syl1);
	if row1>1
		if col1>1
			disp('comparison std syllable is not a vector'), return, end
	syl1=transpose(syl1);
	end
trunc1=floor(.08*length(syl1));
syl1tr=syl1(trunc1:length(syl1)-trunc1);
index=0; %the row-index of the output matrix
filenumb=0;
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

%%%%%%%%%%below here original	
[song, Fs]=evsoundin('',soundfile,'w');

if ~exist([soundfile,'.not.mat'])
	disp([soundfile,'.not.mat not found.'])
	break
end
load([soundfile,'.not.mat']);
filenumb=filenumb+1;
for i=1:length(onsets)
%if i~=1, break,break,end
	index=index+1;
	sylx=song(floor(onsets(i)*Fs/1000):ceil(offsets(i)*Fs/1000));
	sylx=smoother(sylx,Fs,10);
	%plot(sylx);hold;plot(syl1);hold
	
	%find optimal offset to use to compare these 2 vectors.
	[corr,lag]=xcorr(syl1,sylx,floor(.5*length(sylx)));
	maxcorrind=find(corr==max(corr));
	if length(maxcorrind)>1      % in case of tie... take xcorr peak of lowest abs(lag)
		maxcorrdists=(abs(maxcorrind-(.5*(length(corr)+1))));
		maxcorrind=maxcorrind(find(maxcorrdists==min(maxcorrdists)));
	end
	sylxoffs=lag(maxcorrind);
	
	%build new vectors to compute the sylcomp score with. Is deleting elements the right thing to do?
	[row1,col1]=size(syl1);
	[rowx,colx]=size(sylx);
	if row1>1
		if col1>1
			disp('comparison std syllable is not a vector'), break, end
		syl1=transpose(syl1);
	end
	if rowx>1
		if colx>1
			disp('song syllable is not a vector'), break, end
		sylx=transpose(sylx);
	end

	%first and last 8% of each syl won't be used.
	%syl1tr from beginning is 8% truncated.
	truncx=floor(.08*length(sylx));
	sylxtr=sylx(truncx:length(sylx)-truncx);
	%trunc1-truncx is the offset (element to start from) needed for sylx to be at 0 offset.. should be at sylxoffs offset.
	%since leftward shift should be neg number, using truncx-trunc1 for offset.
	%to adjust length-mismatch penalty-- change p. p=1 is low penalty, drop it towards 0 for more.
	p=.03;
	sylxoffs=sylxoffs+(truncx-trunc1);
	if sylxoffs<0
		sylxtr=sylxtr(abs(sylxoffs):length(sylxtr));
		sylxoffs=0;
	elseif sylxoffs>0
		
		sylxtr=[p*(mean(sylxtr)*ones(1,sylxoffs)), sylxtr];
	end
		%modify syllables for score computation.
	offset2=length(sylxtr)-length(syl1tr); %if offset2 is neg, syl1 is longer.
	if offset2<0
		%sylxtr=[sylxtr, (mean(sylxtr)*ones(1,abs(offset2)))];
		%above is less extreme penalties
		sylxtr=[sylxtr, (p*mean(sylxtr)*ones(1,abs(offset2)))];
		syl1trc=syl1tr;
	else	
		%less extreme: syl1trc=[syl1tr, (mean(syl1tr)*ones(1,abs(offset2)))];
		syl1trc=[syl1tr, (p*mean(syl1tr)*ones(1,abs(offset2)))];
	end 	
	
	%now compute score: variance of observed ratios in the 2 smoothed syls.
	scores(index,1)=filenumb;
	scores(index,2)=i;
	if ischar(labels(i))
	scores(index,3)=labels(i);
	else
	scores(index,3)=0;
	end
	%offset2
	%plot(syl1tr);hold;plot(sylxtr)
	%figure;plot(syl1tr(sylxoffs+1:length(syl1tr)+offset2));hold;plot(sylxtr(sylxoffs+1:length(sylxtr)))
	if offset2<0
		%mat1=syl1tr(sylxoffs+1:length(syl1tr)+offset2);
		%mat2=sylxtr(sylxoffs+1:length(sylxtr));
		%old version rats=syl1trc(sylxoffs+1:length(syl1trc)+offset2)./sylxtr(sylxoffs+1:length(sylxtr));
		rats=log2(syl1trc)./log2(sylxtr);
	else
		%old verion rats=syl1trc(sylxoffs+1:length(syl1trc))./sylxtr(sylxoffs+1:length(sylxtr)-offset2);
		rats=log2(syl1trc)./log2(sylxtr); % should get rid of this archaic if loop!
	end
	%here correction should randomly vary for the duration of (offset2+sylxoffs).. what range
	%should it vary across? min(rats) to max(rats)? that is lower penalty than unbounded.

	%%%%%%%%%%%%
	%deleted to allow offsets to be part of sylx
	%randoms=(rand(1,abs(offset2)+sylxoffs));
	%
	%flr=0;
	%%uncomment for <min(rats),max(rats)> version
	%%flr=((1-randoms)*min(rats)); %flr+randoms raises the floor
	%correction=randoms.*(max(rats))^2; %that's random <0,max(rats)^2>
	%correction=flr+correction; %that's random <min(rats),max(rats)>
	
	%rats=[rats, correction];
	%%%%%%%%%%%%%%
	
	scores(index,4)=var(rats);
	%length(correction)
	
end
end %for while loop - (batchfile fgetl)
fclose(fid);
