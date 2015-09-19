function score=sylcomp_s(syl1,syl2)
%[scores]=sylcomp(syl1,syl2)
%syl1 is a vector--smoothed (10ms) syllable to compare syls with
%get syl1 using syl1=getsyl(i,soundfile)
%
%[score] has 1 column, comparisonresultscore.

%make truncated syl1-- will be used for the score calc.

%make var to change amount of trunc:

[row1,col1]=size(syl1);
	if row1>1
		if col1>1
			disp('comparison std syllable is not a vector'), return, end
	syl1=transpose(syl1);
	end
trunc1=floor(.08*length(syl1));
syl1tr=syl1(trunc1:length(syl1)-trunc1);
index=0; %the row-index of the output matrix
	
	i=1;
	filenumb=1;
	labels='x';
	
	index=index+1;
	sylx=syl2;
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
			disp('comparison std syllable is not a vector'), return, end
		syl1=transpose(syl1);
	end
	if rowx>1
		if colx>1
			disp('song syllable is not a vector'), return, end
		sylx=transpose(sylx);
	end

	%first and last 8% of each syl won't be used.
	%syl1tr from beginning is 10% truncated.
	truncx=floor(.08*length(sylx));
	sylxtr=sylx(truncx:length(sylx)-truncx);
	%trunc1-truncx is the offset (element to start from) needed for sylx to be at 0 offset.. should be at sylxoffs offset.
	%since leftward shift should be neg number, using truncx-trunc1 for offset.
	sylxoffs=sylxoffs+(truncx-trunc1);
	if sylxoffs<0
		sylxtr=sylxtr(abs(sylxoffs):length(sylxtr));
		sylxoffs=0;
	elseif sylxoffs>0
		%will keep track of how many elements of new syltr are "fake", but will ignore them first,
		%then get score for the portion of syl after the offset, then change score to reflect the beginning of sylxtr.
		sylxtr=[zeros(1,sylxoffs), sylxtr];
	end
		%modify syllables for score computation.
	offset2=length(sylxtr)-length(syl1tr); %if offset2 is neg, syl1 is longer. 	
	
	%now compute score: variance of observed ratios in the 2 smoothed syls, from a semilog point of view.
	%scores(index,1)=filenumb;
	%scores(index,2)=i;
	%if ischar(labels(i))
	%scores(index,3)=labels(i);
	%else
	%scores(index,3)=0;
	%end
	
	%offset2
	%plot(syl1tr);hold;plot(sylxtr)
	%figure;plot(syl1tr(sylxoffs+1:length(syl1tr)+offset2));hold;plot(sylxtr(sylxoffs+1:length(sylxtr)))
	if offset2<0
		%mat1=syl1tr(sylxoffs+1:length(syl1tr)+offset2);
		%mat2=sylxtr(sylxoffs+1:length(sylxtr));
		rats=(syl1tr(sylxoffs+1:length(syl1tr)+offset2))./(sylxtr(sylxoffs+1:length(sylxtr)));
	else
		rats=(syl1tr(sylxoffs+1:length(syl1tr)))./(sylxtr(sylxoffs+1:length(sylxtr)-offset2));
	end
	%here correction should randomly vary for the duration of (offset2+sylxoffs).. what range
	%should it vary across? min(rats) to max(rats)? that is lower penalty than unbounded.

	
	randoms=(rand(1,abs(offset2)+sylxoffs));
	flr=0;
	%uncomment for <min(rats),max(rats)> version
	%flr=((1-randoms)*min(rats)); %flr+randoms raises the floor
	%correction=randoms.*(max(rats)); %that's random <0,max(rats)>
	correction=randoms.*(max(rats))^2;
	correction=flr+correction; %that's random <min(rats),max(rats)>
	
	rats=[rats, correction];
	score=var(rats);
	%length(correction)
	
