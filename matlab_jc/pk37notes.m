%pk37bk19 perturbations%
300;1;100;2

%  Age on September 12 ~129 days





%  Assay 1 (September 13 - September 26)
%       Pushed middle of song up.
%       Baseline: 913all,914all,915all
%  Assay 2 (October 15 - October 21 at 3pm)
%       Pushed front of song up (failed).
%       Baseline: 1015all
%  Assay 3 (October 21 afternoon - October 28)
%       Pushed middle of song down
%       Baseline: 1015all
%  Assay 4 (November 3 - recovery starts on November 10 - 
%       Tested the role of the white noise on efference copy
%       Baseline: all day 1103 and morning until 3pm on 1104




%September 12 (Fri) - started the template -- templates.mat

%September 15 (Mon) at ~2pm to nightfall - messed around with white noise
%feedback so it is likely to see changes in song here.  Messed around with
%notched and un-notched white noise files.
%September 15 (Mon) at 9pm - started the contingency of 2380Hz (to actually
%begin at 7am the next day).

%September 16 (Tues) at 10:30am-12:30pm - messed up and failed to record 
% data from the bird
%September 16 (Tues) at 4:40pm - changed escape contingency from 2380 to
%2430Hz

%September 17 (Wed) at 1:40-2:00pm - stopped recording song and changed:
%  1. escape contingency from 2430 to 2450 Hz
%  2. template --- changed from original template to OR logic between two
%  templates (first=original template; second=new template that is much
%  more precise and also is tuned to the new avg signal --- i.e. 90th pt in 
%  average of the morning of Sept 17). The result is much more precise with
%  the goal to drive the localized pitch shift.  The std(toff)=7.7ms in the
%  test data set!!! ---- template2.mat
% Turns out that 2:00pm to 4:00pm ---> all song was chunked into the same
% file which I can go back and split up if I need to. For now, I'll just
% use 10:30am-1:30pm as the second time point for sept 17.

%September 18 (Thu) at 11:07amfvals16AMa=findwnoteJC('batch16dawnnotes','a','','',0,[3000 4500],6000,1,'obs0')
%Template 2 is still good, but contingency is too low.  I raised it to
%2500Hz and changed the output folder to 2500conting
%September 18 (Thu) at 5:22pm
%Template 2 is still doing well.  Contingency is possibly too high but I'm
%leaving it where it is for a little while.  Another check at 6:15pm shows
%that it is still a good contingency.

%September 19 (Fri) at 10:54am
%Lowered contingency to 2480Hz b/c bird wasn't improving and few notes were
%above 2500Hz.  Changed output folder to 2480conting.
%5:18pm - made it 2495Hz - template still looks good - output folder
%2495conting

%%%ASILOMAR%%%%

%September 23 (Tues) at 4:35pm
%At 4:35pm, I reduced the contingency to 2440Hz and started the new output
%folder 2440conting92308. 
%At 5:25pm, I lowered it to 2425Hz.
%At 8:13pm, I lowered it to 2400Hz.

%September 24 (Wed) at 10:50am
%Trending back upwards - moved contingency up to 2430Hz.
%Move further up to 2445Hz.
%At around 8pm, moved down to 2335Hz due to low hit rate.

%September 25 (Thu)
%11am - Looks good - kept it at the same contingency - which is ~exactly 
%80th prctile.
%TSANALY - compare 917 with 914 (same time of day --- 11am-12:30pm)
%2pm - there is still great targeting and it is right on the 80th prctile
%9pm - set contingency to 2445 in homes that his circadian rhythm would
%continue and he would learn in the morning.

%September 26 (Fri)
%pk37bk19: Bird seems to be completely not learning despite good targeting.  Thus, I 
%transferred him to SBox 3 in Kris's room and am recording his song in SAP.
%This will allow determination of timecourse of return to baseline.  

%New bird: I would like to set up a new bird and discussed which birds to use with
%Evren.  I have cleared the soundbox used for pk37bk19 for use with this
%putative bird.

%TSA: Great success!  919 data matches with 917 data regarding both
%development of a correlation and the reason for the development of this
%correlation (namely learning ~3 notes after escapes).  Also, I decided to
%do the control by mimicking proportion of hits and using getvals to
%mimic the difference between these values and my pitch contour values.



%October 14
%around 8pm I started recording song in the evtaf-enabled soundbox

%October 15
%1:34pm - tested templa1front
%9pm - decided to set at 2370Hz (around 65th percentile) - targeting was
%around 89 (near front of note) whereas the original was around 112 - this
%is a difference of 23ms.  The false +/- rates are good.

%October 16
A,B-templa1front,ctnrng1front
C-templa2front,cntrng2front
%11:35am - raised to 2380Hz.
%12:10pm - raised to 2390Hz.
%12:20pm - raised to 2400Hz.
%5:15pm - new template (templa2front+ cntrng2front) that is based on the
    %1016Bfiles - hoping to move targeting closer to note onset.
%6:10pm - set to 2410Hz (same folder)

%October 17
A,B,C-templa2front,cntrng2front
%10:30am - lowered to 2400Hz
%12:45am - raised to 2420Hz
%12:57am - raised to 2430Hz
%6:15pm - raised to 2440Hz

%October 21 - despite many efforts it seems that this specific bird cannot
%change the front of its song. Thus I decided to try to drive the song
%negative.  I turned the amp off at 3pm.
%3:30pm - turned the amp back on and started actively pushing down(templa1down)
%8pm - templa2down (templa1down + current low level[batchU])
%pushed it so a good proportion (~30-40% are below 2430)

%October 22 (wednesday) - went huge pushing the bird down and whistling
%until he was at least a standard deviation down by the end of the day.
%Left the contingency set at 2290Hz and went to Princeton Homecoming until
%Monday - the average hit rate was around 50-60% during this time.

%October 27 (monday) 
A - templa2
B,C - templa3front,cntrng3front
%returned to town and the bird was shifted around 2.5
%sigma down (success).  Naturally the targeting had gotten pretty bad by
%then so I made a new template (templa3) starting at around 12pm.  At
%approximately this time, I also lowered the contingency from 2290Hz to
%2250Hz.

%October 28 (tuesday)
B - templa3
%Computer is now on standard time instead of daylight savings time.
%The program crashed overnight, thus between 8am and 11am there was no
%recording.  Luckily the bird appeared to stay fairly well shifted down.

%4:30pm - decided there was no longer any reason to try to shift him down.
%Turned it to closed logic to allow recovery without white noise.

clear LMANdependency
clear deltaVariation
for i=1:3
    LMANdependency(i)=mean(abs(PCnorm(1:83,i)))-mean(abs(PCnorm(84:106,i)));
    deltaVariation(i)=mean(abs(PCnorm(107:162,i)))-mean(abs(PCnorm(1:83,i)));
end
figure;plot(LMANdependency,deltaVariation,'*')





%%%%%%%%%%%%%%   Assay 4 - hipass white noise %%%%%%%%%%%%%%
% highpasswn.m
%
% November 4 - 2pm - made template based on baseline (1103,1004) data
%       at 3:30pm, I started hitting it upwards at 2360Hz (70th prctile)
% November 5 - 12pm - upped it to 2400Hz
%          4pm - lowered back to 2380Hz 
%          4:30pm - new template -- templa2, cntrng2 - better targeting
% November 8 - 12pm - new template -- templa3,cntrng3
%  Problem with targeting this note because it's so oddly shaped.
%
% November 10 - turned amp off (actually closed logic) at 10:30 am because
% the previous day's data looked great.
% Recovery to baseline occurs.