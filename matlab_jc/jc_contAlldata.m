function [p]=jc_contAlldata(ZFnorm,BFnorm,stds,contingency,chunk)
for i=1:length(ZFnorm)
    falloffZF(i)=jc_contingency(ZFnorm(i).data,1,stds,contingency,chunk);
end
for i=1:length(BFnorm)
    if i==4 || i==5 || i==6
        chspec=0;
    else
        chspec=1;
    end
    falloffBF(i)=jc_contingency(BFnorm(i).data,chspec,stds,contingency,chunk);
end

%figure;hold on;xlim([0 3]);
%plot(1,falloffBF,'*');
%plot(2,falloffZF,'*');
M(1)=mean(falloffBF);
M(2)=mean(falloffZF);
E(1)=std(falloffBF)/sqrt(length(falloffBF));
E(2)=std(falloffZF)/sqrt(length(falloffZF));
%figure;errorbar(M,E);
[h,p]=ttest2(falloffBF,falloffZF);

