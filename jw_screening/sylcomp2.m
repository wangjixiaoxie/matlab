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

%%if ~exist([soundfile,'.not.mat'])
%%	disp([soundfile,'.not.mat not found.'])
%%	break
%%end
%%load([soundfile,'.not.mat']);
filenumb=filenumb+1;
for i=1:length(onsets)
%if i~=1, break,break,end
	index=index+1;
	sylx=song(floor(onsets(i)*Fs/1000):ceil(offsets(i)*Fs/1000));
	sylx=smoother(sylx,Fs,10);
	%plot(sylx);hold;plot(syl1);hold
	
	scores(index,4)=sylcomp_s(syl1,sylx);
	scores(index,1)=filenumb;
	scores(index,2)=i;
	if ischar(labels(i))
	scores(index,3)=labels(i);
	else
	scores(index,3)=0;
	end
	
end
end %for while loop - (batchfile fgetl)
fclose(fid);
