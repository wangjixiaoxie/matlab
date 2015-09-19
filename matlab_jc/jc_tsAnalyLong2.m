function [xax,store]=jc_tsAnalyLong2(norm,firstpoint,lastpoint)

varPshift=0;
varSin=0;
%put the data into a matrix for pca
store=zeros(30);


for i=1:15 %25 %size(norm,2)

    a=norm(firstpoint:lastpoint,i)';
    for kk=1:length(a)
        t(kk)=kk;
    end
    varPshift=varPshift+abs(mean(a))*length(a); %area under the pitch shift
    b(i,:)=a-mean(a);
    [mb,resnorm,residual] = lsqcurvefit(@jc_sinfit,[0 0.001 0 0],t,b(i,:));
    c=residual;

    hh=jc_fft(xcorr(c));
    
    for k=1:30
        store(k)=store(k)+hh(k);
    end
end

xax=zeros(1,30);
for i=1:30
    xax(i)=(1/(i/2))*length(norm(:,1))/8;
end

store=store(2:29);
xax=xax(2:29);