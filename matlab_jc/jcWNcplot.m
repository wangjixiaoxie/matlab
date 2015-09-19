%read me - plot of targeting vs. peak change
%100808
%D:\jcharles\ContingencyTest\bk20
%structure BK20BK45 in 1008_BK20BK45contingencydata.mat - saved

%



%Find the file of all catch trials in the wnon period
%Get the offsets
[fvalsT,FFstatsT,avgT,toffsetT,pitchT]=summary_statsTW('batchcomb','batchcombnotes');
%How to align the original pitch curves with toffsets:
    %Problem because using findwnoteJC
figure;plot(avgT(1:1320));hold on;plot(avgbk7,'r')
figure;plot(xcov(avgT,avgbk7)) %Take peak of this curve -- 'g'
g=74;
figure;plot(avgbk1(g+1:g+1000))
hold on;plot(toffsetT,2450,'*')




%%%%Plot the stuff as standard deviations%%%%
%1. Look at get_avn to determine note onset and offset times
%Actually the best way is to look at where the pitch curve levels off
B.notonset=150;
B.notoffset=800;
B.notlength=B.notoffset-B.notonset+1;
%2. Choose the shifted data and the baseline data
B.avgBaseline=avgbk1+avgbk2;
B.avgBaseline=B.avgBaseline./2;
B.pitchBaseline=[pitchbk1,pitchbk2];
B.stdBaseline=std(B.pitchBaseline');
B.avgShifted=avgbk7;
%avgShifted=avgShifted./6;
%Calculate and plot the trace
B.SDtrace=(B.avgShifted-B.avgBaseline)./B.stdBaseline;
B.SDTplot=B.SDtrace(B.notonset:B.notoffset);
for i=1:B.notlength
    B.xax(i)=i/8;
end
figure;plot(B.xax,B.SDTplot)
%%%%%%
B.toffset=(toffsetT-B.notonset)/8; %puts it in ms
hold on;plot(B.toffset,-2,'*')