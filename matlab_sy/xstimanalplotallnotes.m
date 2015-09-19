%Xstimanalscript
%analysis of acute effects of X stimulation
%first thing to do is just to compare the average note to catch note for
%frequency and amplitude...(bracket amplitude)

%also want to compare total number of songs [bracket]
%proportion of as in all songs
%distribution of number of repeats

%to do this for each batch
clear strcmd;
clear bt;
dirchange=0;
dir{1}='doyale4/twarren/bu92bu1/stim2'
%dir{2}='doyale4/twarren/bu92bu1/stim1'
%dir{3}='doyale4/twarren/bu92bu1'

trigcheck=1;





%compare the frequency of the notes
tbinshft=0.01;
NFFT=1024;%number of data points to FFTstrcmd=strcat('!cd ' dir{i})
fbins=[3000,4500; 6000,9000];
save BINS_B NFFT fbins tbinshft
% frequency analysis just for 'b'
load BINS_B

vals=getvals(fv,2,'TRIG');
vals=[vals]

indtrig=find(vals(:,3)==1);
indnotrig=find(vals(:,3)==0);
figure
plot(vals(indtrig,1)-fix(vals(1,1)),vals(indtrig,2)/2,'r.')
hold on;
plot(vals(indnotrig,1)-fix(vals(1,1)),vals(indnotrig,2)/2,'k.');

