function [xax,ts]=jc_tsAnalyLong(norm,firstpoint,lastpoint,filetype,bins)

varPshift=0;
varSin=0;
%put the data into a matrix for pca
Xlength=45;
store=zeros(Xlength);


for i=1:size(norm,2) %25 %size(norm,2)

    a=norm(firstpoint:lastpoint,i)';
    for kk=1:length(a)
        t(kk)=kk;
    end
    varPshift=varPshift+abs(mean(a))*length(a); %area under the pitch shift-this allows quantification of the pitch shift whereas doing xcov would just remove it
    b(i,:)=a-mean(a);
    hh=jc_fft(xcorr(b(i,:)));
    for k=1:Xlength
        store(k)=store(k)+hh(k);
    end
end

xax=zeros(1,Xlength);
if filetype=='w'
    c=1/((44100/bins)*(1/1000));
else 
    c=1/((32000/bins)*(1/1000));
end

for i=1:Xlength
    xax(i)=(2/i)*(lastpoint-firstpoint+1)*(c);
    store(i)=store(i)/xax(i);
end


ts=store(1:Xlength);
xax=xax(1:Xlength);