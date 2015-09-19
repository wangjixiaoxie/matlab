function [output]=jc_tsAnaly(norm,firstpoint,lastpoint,filetype,bins)

varPshift=0;
varSin=0;
%put the data into a matrix for pca

%period of SinFit in points
BumpWidthMS=300;
if filetype=='w'
    c=(44100/bins)*(1/1000);
else 
    c=(32000/bins)*(1/1000);
end
SineWidth=BumpWidthMS*2*c;
SinePeriod=(2*pi)/SineWidth;


for i=1:13 %25 %size(norm,2)

    a=norm(firstpoint:lastpoint,i)';
    for kk=1:length(a)
        t(kk)=kk;
    end
    varPshift=varPshift+abs(mean(a))*length(a); %area under the pitch shift
    b(i,:)=a-mean(a);
    [mb,resnorm,residual] = lsqcurvefit(@jc_sinfit,[0 SinePeriod 0 0],t,b(i,:));
    c(i,:)=residual;
    ff=jc_sinfit(mb,t);
    varSin=varSin+sum(abs(ff));
    periodSin(i)=(2*pi)/mb(2);
    ampSin(i)=mb(1);
end
normampSin=(abs(ampSin)/max(abs(ampSin)));
normampSin=normampSin*(1/mean(normampSin));
perSin=periodSin.*normampSin;
perSin=mean(perSin);
perSin=perSin/16; %Convert to ms
[PC, SCORE, LATENT]=princomp(c);
[timescaleSin,xax]=jc_tscaleNOsmooth(PC,LATENT);


%%%Linear%%%%%
varLin=0;
cc=b;
clear b;
clear a;
for i=1:13 %25
    a=cc(i,:);
    slope1=(a(length(a))-a(1))/(length(a)-1);
    varLin=varLin+0.5*abs(a(length(a))-a(1))*length(a);
    for j=1:length(a)
        b(i,j)=a(j)-(j-1)*slope1;
    end
    b(i,:)=b(i,:)-b(i,1);
end
[PC2, SCORE2, LATENT2]=princomp(b);
[timescaleLin,xax2]=jc_tscaleNOsmooth2(PC2,LATENT2);

output.timescaleSin=timescaleSin;
output.timescaleLin=timescaleLin;
output.xax=xax;
output.varPshift=varPshift;
output.varLin=varLin;
output.varSin=varSin;
output.perSin=perSin;