function [PC,SCORE,LATENT,varPshift,varSlope]=jc_pcaX(norm,first,last)

varPshift=0;
varSlope=0;
%put the data into a matrix for pca
for i=1:25 %size(norm,2)
    a=norm(first:last,i)';
    varPshift=varPshift+abs(mean(a))*length(a); %area under the pitch shift
    a=a-mean(a);
    slope1=(a(length(a))-a(1))/(length(a)-1);
    varSlope=varSlope+0.5*abs(a(length(a))-a(1))*length(a);
    for j=1:length(a)
        b(i,j)=a(j)-(j-1)*slope1;
    end
    b(i,:)=b(i,:)-b(i,1);
end

% PCA
% PC = the principalcomponents (as columns)
% SCORE = the transformed matrix a.
% LATENT = the eigenvalues
[PC, SCORE, LATENT]=princomp(b);


