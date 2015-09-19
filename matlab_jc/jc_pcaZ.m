function [PC,SCORE,LATENT,varPshift,varSin]=jc_pcaZ(norm)

varPshift=0;
varSin=0;
%put the data into a matrix for pca
for i=1:25 %size(norm,2)

    a=norm(300:650,i)';
    for kk=1:length(a)
        t(kk)=kk;
    end
    varPshift=varPshift+abs(mean(a))*length(a); %area under the pitch shift
    b=a-mean(a);
    [mb,resnorm,residual] = lsqcurvefit(@jc_sinfit,[0 0.01 0 0],t,b);
    c(i,:)=residual;
    ff=jc_sinfit(mb,t);
    varSin=varSin+sum(abs(ff));
end
    
% PCA
% PC = the principalcomponents (as columns)
% SCORE = the transformed matrix a.
% LATENT = the eigenvalues
[PC, SCORE, LATENT]=princomp(c);


