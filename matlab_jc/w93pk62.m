% w93pk62
% Evren worked with this bird previously - 900 days old

% 12.08 - noon - placed in soundbox (wn was accidentally on, but few hits)
    % wn hit 6 and up for another bird - hit occassionally
% 12.09 - 7:30pm - wn off
    
% 12.12 - dawn - wn on
%        - 2:36pm - someone turned amp off accidentally?
%        - 5:03pm - I turned amp back on        
% 12.13 - wn on
% 12.14 - wn off at 10:38am

% 12.18 - dawn - hit repeats 1-7 (65%) to increase.
%       - 11:17am - hit repeats 1-5, lower wn volume a bit (make sure he sings through and gets escapes)

% 12.20 - dawn - hit 1-5 to increase repeat number
%   used template1219reps.mat (MIN=2, MAX=3, TH=2.5, synshifts)
% --- looks like learning - but is also doing branch point learning...
% 

% 03.05.11 - wn on at dawn to cause pitch shift - hit below 3540 Hz (hit first A in repeat)
% 03.06.11 - wn off at dark - great pitch shift

% 03.28.11 - POSTLESION
    % There was an appearance of subharmonics and increase in pitch of high stack following
    % the lesion.  There was no change in the length of the targeted repeat.

    % 4:30pm - new template for pitch shift - 1sec refractory period
   % 9pm - adjusted cntrng settings, hit below 3710Hz (mean+20Hz), 1second
   % refractory period
% 03.29.11 - adjusted to 3730Hz around noon
% 03.30.11 - adjusted to 3715Hz around 11:15am
%   9pm - wn off

% 03.31.11 - around 4pm - transferred to new computer
% Looks like repeat increase (possibly from hitting first A in pitch shift) - wait to let it stabilize...

% 04.01.11 - around 7pm - new template to hit repeats 5 to 25 - ready for WN
    % template adjusted to recognize notes with lower pitch (which often
    % occur later in the train of repeats)
% 04.03.11 - 13/36 (36%) have six or fewer repeats - plan to hit repeats 7 to 25 with WN
%   WN on at 9pm
% 04.04.11 - WN on - repeat reduce - hit repeats 7 to 25 with WN
% 11am - Looking good


%%% pitch shift pre
%%%
figure;hold on;
win1=180:240;
plot(runningaverage(timing3(fvalsPRE),15),runningmedian(mean(pitchPRE(win1,:)),15),'.')
plot(runningaverage(timing3(fvalsDUR),15),runningmedian(mean(pitchDUR(win1,:)),15),'r.')
plot(runningaverage(timing3(fvalsPOST),15),runningmedian(mean(pitchPOST(win1,:)),15),'.')

[timevals1210pre,distnt1210pre] = jcrepeatdist('a','batch1210notes');
[timevals1211pre,distnt1211pre] = jcrepeatdist('a','batch1211notes');
%
[timevals1212wn,distnt1212wn] = jcrepeatdist('a','batch1212notes');
[timevals1213wn,distnt1213wn] = jcrepeatdist('a','batch1213notes');
[timevals1214wn,distnt1214wn] = jcrepeatdist('a','batch1214notes');
%
[timevals1214post,distnt1214post] = jcrepeatdist('a','batch1214notes');
[timevals1215post,distnt1215post] = jcrepeatdist('a','batch1215notes');
[timevals1216post,distnt1216post] = jcrepeatdist('a','batch1216notes');
[timevals1217post,distnt1217post] = jcrepeatdist('a','batch1217notes');
%
[timevals1218wn,distnt1218wn] = jcrepeatdist('a','batchnotes');
%
[timevals1220wn,distnt1220wn] = jcrepeatdist('a','batch1220notes');
[timevals1221wn,distnt1221wn] = jcrepeatdist('a','batch1221notes');
[timevals1222post,distnt1222post] = jcrepeatdist('a','batch1222notes');
[timevals1223post,distnt1223post] = jcrepeatdist('a','batch1223notes');


figure;hold on;
window1=10;
plot(runningaverage(timevals1210pre,window1),runningaverage(distnt1210pre,window1))
plot(runningaverage(timevals1211pre,window1),runningaverage(distnt1211pre,window1))
plot(runningaverage(timevals1212wn,window1),runningaverage(distnt1212wn,window1),'r')
plot(runningaverage(timevals1213wn,window1),runningaverage(distnt1213wn,window1),'r')
%plot(runningaverage(timevals1214wn,window1),runningaverage(distnt1214wn,window1),'r')
plot(runningaverage(timevals1214post,window1),runningaverage(distnt1214post,window1))
plot(runningaverage(timevals1215post,window1),runningaverage(distnt1215post,window1))
plot(runningaverage(timevals1216post,window1),runningaverage(distnt1216post,window1))
plot(runningaverage(timevals1217post,window1),runningaverage(distnt1217post,window1))
plot(runningaverage(timevals1218wn,window1),runningaverage(distnt1218wn,window1),'r')
plot(runningaverage(timevalsX,window1),runningaverage(distntX,window1),'b')
plot(runningaverage(timevals1220wn,window1),runningaverage(distnt1220wn,window1),'r')
plot(runningaverage(timevals1221wn,window1),runningaverage(distnt1221wn,window1),'r')
plot(runningaverage(timevals1222post,window1),runningaverage(distnt1222post,window1),'b')
plot(runningaverage(timevals1223post,window1),runningaverage(distnt1223post,window1),'b')


