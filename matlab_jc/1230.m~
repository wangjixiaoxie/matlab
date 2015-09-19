%%%% Compare actual vs. predicted

for i=1:17
[actshift(i).dat,normshift(i).dat,pcnormshift(i).dat]=jc_actualshift(Alldata2(i),1,0.8,24);
end

up=[1 1 0 0 1 0 1 0 1 1 1 1 0 1 1 0 1];
for n=1:17
    if up(n)==1
        prc=70;
    else
        prc=30;
    end
    figure;plot(normshift(n).dat,'r')
    %change 30 to 20
    avf=ContingSim((tfinal(n).data),Alldata2(n).baselineAC,prc);
    ga=avf-mean(Alldata2(n).baselineAC');
    % change min to max
    if up(n)==1
        nga(n).dat=ga./max(ga(Alldata2(n).startnote:Alldata2(n).endnote));
    else
        nga(n).dat=ga./min(ga(Alldata2(n).startnote:Alldata2(n).endnote));
    end
    goob(n).data=nga(n).dat(Alldata2(n).startnote:Alldata2(n).endnote);
    
    hold on;plot(goob(n).data,'b')
end



% for i=13:17
%     tfinal(i).data=toffsets(i).data((find(toffsets(i).data>Alldata2(i).startnote & toffsets(i).data<Alldata2(i).endnote)));
% end


%%%%%
% ds4 - ac622wnon_3 (exp5)
% ds3 - ac6408-2(exp8)