function jc_contAlldataAV(ZFnorm,BFnorm,stds,contingency,chunk)
%aa=250-chunk(1);
longer=550;
averZ=zeros(1,longer);
nZ=ones(1,longer);
averB=zeros(1,longer);
nB=ones(1,longer);
for i=1:length(ZFnorm)
    [av]=jc_contingencyAV1(ZFnorm(i).data,1,stds,contingency,chunk);
    for j=1:length(av)-60 %-extra 10 to avoid crap in BF (resampling)
        averZ(j)=averZ(j)+av(50+j);
        nZ(j)=nZ(j)+1;
    end
end
averZt=averZ./nZ;


for i=1:length(BFnorm)
    if i==4 || i==5 || i==6
        chspec=0;
    else
        chspec=1;
    end
    [av]=jc_contingencyAV1(BFnorm(i).data,chspec,stds,contingency,chunk);
    for j=1:length(av)-60 %-extra 10 to avoid crap in BF (resampling)
        averB(j)=averB(j)+av(50+j);
        nB(j)=nB(j)+1;
    end
end
averBt=averB./nB;
%center=round(mean(chunk));
%diff=averB(center-50)-averZ(center-50);
for i=1:longer
    xax(i)=i*(4/32);
end
hold on; plot(xax,averBt,'k');plot(xax,averZt,'r','LineWidth',2)

%figure;hold on;xlim([0 3]);
%plot(1,falloffBF,'*');
%plot(2,falloffZF,'*');
%M(1)=mean(falloffBF);
%M(2)=mean(falloffZF);
%E(1)=std(falloffBF)/sqrt(length(falloffBF));
%E(2)=std(falloffZF)/sqrt(length(falloffZF));
%figure;errorbar(M,E);
%[h,p]=ttest2(falloffBF,falloffZF);


