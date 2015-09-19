function [PC,SCORE,LATENT]=jc_pca519(arrayfile,f)

%Smooth the data with a filter so pca doesn't get a bunch of noise
for i=1:length(arrayfile)
    arrayfile(i).pitches=jc_wfilt(arrayfile(i).pitches);
end

%Align
for i=1:length(arrayfile)
    [holder]=SmoothData(f(i).datt,32000,1);
    smooth(i).smoothed=holder;
end
for i=1:length(smooth(1).smoothed)
    sumnote(i)=0;
    for j=1:length(smooth)
        sumnote(i)=sumnote(i)+smooth(j).smoothed(i);
    end
    avgnote(i)=sumnote(i)/length(smooth);
end
for i=1:length(arrayfile)
    h=xcov(avgnote,smooth(i).smoothed);
    [peak,index]=max(h);
    shift=(index-6400)/4;
    starting=round(280-shift);
    ending=starting+400;
    forpca(i).pitches=arrayfile(i).pitches(starting:ending);
end



%put the data into a matrix for pca
for i=1:length(arrayfile)
    a(i,:)=forpca(i).pitches;   
end

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
