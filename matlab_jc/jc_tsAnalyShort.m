function [output]=jc_tsAnalyShort(norm,firstpoint,lastpoint,filetype,bins)

varPshift=0;
varSin=0;
Xlength=45;

%%%%%%%%%%%%xaxis%%%%%%%%%%%
xax=zeros(1,Xlength);
if filetype=='w'
    c=1/((44100/bins)*(1/1000));
else 
    c=1/((32000/bins)*(1/1000));
end
for i=1:Xlength
    xax(i)=(2/i)*(lastpoint-firstpoint+1)*(c);
end

%%%%%yaxis%%%%%%%%%%%%%%%%%%%

%period of SinFit in points
BumpWidthMS=300;
if filetype=='w'
    c=(44100/bins)*(1/1000);
else 
    c=(32000/bins)*(1/1000);
end
SineWidth=BumpWidthMS*2*c;
SinePeriod=(2*pi)/SineWidth;

store=zeros(1,Xlength);
for i=1:size(norm,2)
    a=norm(firstpoint:lastpoint,i)';
    for kk=1:length(a)
        t(kk)=kk;
    end
    varPshift=varPshift+abs(mean(a))*length(a); %area under the pitch shift
    b=a-mean(a);
    [mb,resnorm,residual] = lsqcurvefit(@jc_sinfit,[0 SinePeriod 0 0],t,b);
    residuals(i,:)=residual;
    ff=jc_sinfit(mb,t);
    hh=jc_fft(xcorr(residuals(i,:)));
    for k=1:Xlength
        store(k)=store(k)+hh(k);
    end
    varSin=varSin+sum(abs(ff));
    periodSin(i)=(2*pi)/mb(2); 
    ampSin(i)=mb(1);
end
normampSin=(abs(ampSin)/max(abs(ampSin)));
normampSin=normampSin*(1/mean(normampSin));
perSin=periodSin.*normampSin;
perSin=mean(perSin);
perSin=perSin/(2*c); %Convert to ms

timescale=store;

output.timescale=timescale;
output.xax=xax;
output.varPshift=varPshift;
output.varSin=varSin;
output.perSin=perSin;