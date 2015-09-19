function [PC,SCORE,LATENT]=jc_pca430(pitch_data)
thresh=jc_stdplot(pitch_data(1).pitches,100);
a=pitch_data(1).pitches(thresh+100:thresh+500); % initialize the matrix
for n=2:length(pitch_data)
    thresh=jc_stdplot(pitch_data(n).pitches,100);
    a(n,:)=pitch_data(n).pitches(thresh+100:thresh+500);   %cut so you get useful info.
end
%Center examples to make E[a]=0
%http://www.dtek.chalmers.se/~d95danb/ocr/pca_new.html
dimension=401;   %Because we're going from 15 to 115
b=mean(a,1);
c=ones(dimension,1);
d=b*c;
a=a-d;

% PCA
% PC = the principalcomponents (as columns)
% SCORE = the transformed matrix a.
% LATENT = the eigenvalues
[PC, SCORE, LATENT]=princomp(a);