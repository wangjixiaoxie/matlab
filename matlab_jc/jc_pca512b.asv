function [PC,SCORE,LATENT,a]=jc_pca512b(pitch_data)

%Smooth the data with a filter so pca doesn't get a bunch of noise
for i=1:length(pitch_data)
    smoothed(i).pitches=jc_wfilt(pitch_data(i).pitches);
end

%Align the notes and cut out the part of the note that you are interested in (the stack)
thresh=jc_getthresh(smoothed(1).pitches,100);
a=smoothed(1).pitches(thresh+50:thresh+400); % initialize the matrix
for n=2:length(pitch_data)
    thresh=jc_getthresh(smoothed(n).pitches,100);
    a(n,:)=smoothed(n).pitches(thresh+50:thresh+400);   %cut so you get useful info.
end

%Center examples to make E[a]=0
%http://www.dtek.chalmers.se/~d95danb/ocr/pca_new.html
%dimension=351;   %Because we're going from thresh+100 to thresh + 500 and dimensions must agree
%b=mean(a,1);
%c=ones(dimension,1);
%d=b*c;
%a=a-d;

% PCA
% PC = the principalcomponents (as columns)
% SCORE = the transformed matrix a.
% LATENT = the eigenvalues
[PC, SCORE, LATENT]=princomp(a);