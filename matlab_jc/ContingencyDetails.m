%read me - plot of targeting vs. peak change
colormapper.m
%All saved in D:\jcharles\ContingencyTest
%%%%%pu34 (in folder E:) -- 2000 to 2700
%Data location E:\pu34
%Baseline: ac711\batchJC
%shifted: ac716\batch16.catch
%Targeting: wnonpu34\batchTOFFS (batch.catch)
%           ac714\batchTOFFS (batch14.catch)
%All data analysis done on the fly.

%%%%%pk32bk28exp2 (Tim's)%%%%%%%   --- 2200 to 2800
REDONE
%All data has been transferred into
%D:\jcharles\ContingencyTest\pk32bk28_exp2
%Baseline: ac427\batchJC2 AND screen\batch25notes 
%Shifted: ac502\batch04JC (batch02 had yet to shift beyond 1sd
%Targeting (from oriole4 Tim): wnon428\batchTOFFS (batch.catch)
fvalsBaseline=findwnote4('batchJCnotes','a','','',0,[2800 3500],8500,1,'obs0',0);
[fvalsTOFFSwn,FFstatsTOFFSwn,toffsetTOFFSwn]=summary_statsTW2('batchTOFFS','batchTOFFSnotes');
%%%%%


PK20R49 
%D:\jcharles\ContingencyTest\pk20r49
%Baseline: \ac20908\batchJCnotes
%Shifted: \ac21608\batchJCnotes
%Targeting: \ac21208ampon\batchJCnotes
%Limits: 2800 to 3500
%Note: 'a'
%


PK32BK28exp2
%Baseline: 427ac (values from data file 902pk32ds2.mat)
%Shifted: 502 (values from the same data file)
%Targeting: 428 and 429 (all data) - from wnon428 folder
%To align the targeting with the average pitch curves, I took the pitch
%curve from a few notes in the wnon428 folder (30ms in findwnoteJC -- 
%240 points added to toffset in summary_statsTW1).  Alignment suggested to
%shift the baseline and shifted averages over by 154pts x=x(155:1500);
%ReferenceDatapk32bk28 was modified today (100908) and might be useful.

%%%%%pk37 end of note%%%%%
PK37BK19end
%Baseline for pk37END: 913-915 all data
%Shifted for pk37END: 920 and 921 all data
%Targeting for pk37END: 919 all data - because final day with much shifting

%%%%%bk50 end of note%%%%%
BK50W18end
%Baseline: 927 through 929 (all data)
%Shifted: 1002 and 1003 (all data)
%Targeting: 1001A (when most of the shift occurred)

%%%%%bk50 beginning of note%%%%
%Baseline: Sunday PM (1005C) and all day monday (1006A,B,C)
%Shifted: 1010all, 1011all, 1012A,B


%%%%%bk20 (Tim) note B %%%%%%%%%%%
%Done 100808
%See 90108contingency.txt
%BASELINE - 1 and 2
%SHIFTED - 7
%TARGETING - wnon713
%1. preAC (ac710)
%2. preAC (ac710-2)
%3. preMUSC (500mu711)
%4. shiftday1 (712)
%5. shiftday2 (713)
%6. shiftday3 (714)
%7. ac715

%%%%%bk20 (Tim) note A %%%%%%%%%%%
%All segmentation at 30;5;500;2
%a note is first a in iaaaaa repeats
%b note is first b in bbaaaa repeats
%Note a: [2900 3700]
%Note b: [2100 2800]
%Targeting: D:\jcharles\ContingencyTest\bk20\wnon713\batchWN
    %This batch includes all of 714 and early 715 data and is all of the
    %data that Tim has specified by *.keep
%Baseline: D:\..\ac711\batchBASE ---most recent song that has the proper
%amplitude calibration and thus can be segmented properly --- it is all the
%song from the afternoon of 711
%Shifted: D:\...\ac715\batchSHIFTED (all the catch trials that have notes i.e. .keep
%from the 15th and the 16th


%Done 100908
%SHIFTED:
%D:\jcharles\ContingencyTest\bk20\ac715 - post data for note 'a'
    %I went in and labeled c's
    %batchJCconting
    %'cAaaaaaa'
%BASELINE:
%D:\jcharles\ContingencyTest\bk20\ac711
    %I went in and labeled again
%The notes in the different folders had different segmentation, so it was
%necessary to calibrate for this by comparing the average notes. Ultimately
%I decided to shift the avg711 pitch curve over 7.5ms (60pts) to match the
%avg715 pitch curve.
%TARGETING DATA:
%Same as baseline data

%100808
%D:\jcharles\ContingencyTest\bk20
%structure BK20BK45 in 1008_BK20BK45contingencydata.mat - saved

%100908
%


fvalsPRE=findwnoteJC('batchJCcontingnotes','a','c','a',0,[2000 2700],8500,1,'obs0',0);
for i=1:length(fvalsPRE)
    shifted(i,:)=fvalsPRE(i).datt;
end
[pitch1,avg1]=jc_pitchmat1024(shifted(:,1:5000),1024,1020,1,2800,3600,[1 2 3],'obs0',1);
%Find the file of all catch trials in the wnon period
%Get the offsets
[fvalsTOFFS,FFstatsTOFFS,toffsetTOFFS]=summary_statsTW2('batchJC','batchJCnotes');
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
B.notonset=250;
B.notoffset=1000;
B.notlength=B.notoffset-B.notonset+1;
%2. Choose the shifted data and the baseline data
B.avgBaseline=avg;
B.avgBaseline=B.avgBaseline;
%B.pitchBaseline=[pitch1];
B.stdBaseline=std(B.pitchBaseline(155:1500,:)');
B.avstd=mean(B.stdBaseline);
B.avgShifted=B.avg7prime;
%avgShifted=avgShifted./6;
%Calculate and plot the trace
B.SDtrace=(B.avgShifted-B.avgBaseline)./B.stdBaseline;
B.SDTplot=B.SDtrace(B.notonset:B.notoffset);
clear B.xax;
for i=1:B.notlength
    B.xax(i)=i/8;
end
figure;plot(B.xax,B.SDTplot)
%%%%%%
B.toffset=(toffsetT-B.notonset)/8; %puts it in ms
hold on;plot(B.toffset,2,'*')

%Calculate center of mass of shift and compare to mean and std of targeting
B.meanTarget=mean(B.toffset); %more informative to have median?
B.stdTarget=std(B.toffset);
%B.cmassShift=sum((B.xax.*B.SDTplot))./sum(B.SDTplot);
