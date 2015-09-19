%look at the two different b's and their frequency dependence in the
%screened data so i will use the batch file, batch45.keep.rand

%first make the avn spectrogram of 'b'- first one
[avn,t,f]=get_avn('batch45.keep.rand','b',0.2,0.2,'','b','obs0');
imagesc(t,f,log(avn));syn;ylim([0,1e4]);

%make avn spectrogram of 'c-second one
[avn2,t,f]=get_avn('batch45.keep.rand','b',0.2,0.2,'b','','obs0');
imagesc(t,f,log(avn2));syn;ylim([0,1e4]);


tbinshft=0.01;
NFFT=1024;%number of data points to FFT
fbins=[1000,2000; 2500,3750; 4000,5500; 5700,7250];
save BINS_B NFFT fbins tbinshft

% frequency analysis just for 'b'

load BINS_B
bt='batch45.keep.rand';
NT='b';PRENT='';PSTNT='b';
fv=findwnote4(bt,NT,PRENT,PSTNT,tbinshft,fbins,NFFT,1,'obs0');
%plot the histogram of all 'b's
plt_ff2(fv,tbinshft,fbins,'b',4,30);

%now for 'c'
NT='b';PRENT='b';PSTNT='';
fv2=findwnote4(bt,NT,PRENT,PSTNT,tbinshft,fbins,NFFT,1,'obs0');
%plot the histogram of all 'c's
plt_ff2(fv2,tbinshft,fbins,'b',4,30);


%plot freq with time for 'b'
%to look at values for 4th band
vals=getvals(fv,1,'TRIG');
NN=1
size(vals)

%there were notes 382'b' and then three things out put, time index, the frequency, and whether it triggered or not 
%plots the frequency of the 1st harmonic with the time passed since the
%first song sung
plot(vals(:,1)-fix(vals(1,1)),vals(:,2)/2,'ro')
% now for 'c'
vals2=getvals(fv2,1,'TRIG');
plot(vals2(:,1)-fix(vals2(1,1)),vals2(:,2)/2,'ro');

%get mean fund frequency for 'b'
mean(vals(:,2))/2

ans =

  799.0855
  
%get mean fund frequency for 'c'
mean(vals2(:,2))/2
ans =

  799.0552
  
%so those are pretty close

