function [PC,SCORE,LATENT]=jc_pca616(input)


%put the data into a matrix for pca

a=input';


% PCA
% PC = the principalcomponents (as columns)
% SCORE = the transformed matrix a.
% LATENT = the eigenvalues
[PC, SCORE, LATENT]=princomp(a);





%Center examples to make E[a]=0
%http://www.dtek.chalmers.se/~d95danb/ocr/pca_new.html
%dimension=351;   %Because we're going from thresh+100 to thresh + 500 and dimensions must agree
%b=mean(a,1);
%c=ones(dimension,1);
%d=b*c;
%a=a-d;
