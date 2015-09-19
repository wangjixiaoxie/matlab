for i=31:37
    days=Exp1(i).predays;
    mnpre(i-30)=mean(Exp1(i).replong(days(find(days<5))));
    mnpost(i-30)=mean(Exp1(i).replong(days(find(days>5))));
end

% last index before gap
lastind=[898 53 91 116 113 109 296];
for i=1:7
    [h,p]=ttest2(Exp(i+30).rplength{1}(1:lastind(i)),Exp(i+30).rplength{1}(lastind(i)+1:end))
end

