function [xax,store]=jc_tsAnalyMK(norm,firstpoint,lastpoint,filetype,bins)

varPshift=0;
varSin=0;
%put the data into a matrix for pca

Xlength=200;
store=zeros(Xlength);

for i=1:size(norm,2)

    a=norm(firstpoint:lastpoint,i)';
    for kk=1:length(a)
        t(kk)=kk;
    end
    %varPshift=varPshift+abs(mean(a))*length(a); %area under the pitch shift
    %b(i,:)=a-mean(a);
    hh=jc_fft(xcorr(a));
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

%store=store(1:55);
%xax=xax(1:45);