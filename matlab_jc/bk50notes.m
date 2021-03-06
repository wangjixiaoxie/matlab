%bk50w18 notes
%    160 days old on Sept 27

% Assay 1 (September 27 - October 4)
%   Goal: shift up rear of note
%   Baseline: September 27 - September 29
%   Shifting: September 30 - October 4

%October 5a,b - return to baseline

% Assay 2 (October 5 - October 12B)
%   Goal: shift up front of note
%   Baseline: October 5c and all of October 6
%   Shifting: October 7 - October 12B

% Assay 3 (October 12-13
%   Goal: directed TSA experiment
%   Applied white noise while song returned to baseline



%%%%BASELINE%%%%
%2 per day
%A - 7am-around 1pm or 2pm
%B - 2pm until 9pm (except monday - stopped at 4pm when started wn testing)
%%%%SHIFTING%%%
%3 per day
%A - 7am until around noon
%B - noon until around 4pm
%C - 4pm until nightfall
%exceptions when I change the contingency or template

%Friday 092608
    %took bk50w18 at 7:10pm at put into cage 4 and started recording song in
    %folder initscreen_092608
    
%Sunday 092808
    %3pm - changed name on config file to bk50w18 and saved new config file
    %Labeled at 30;5;100;2
    %templa=mean(avB(1:2:256,[90]),2);
    %4pm - loaded template, made config file and changed folder to
    %92808_tmptest
    
%Monday 092908
    %More baseline recording - targeting is very good.
    %11:30am - set number of matches to 3 because targeting still will look
    %good and it looks (based on ContingSim) that the early-middle part of
    %the note is more pliable - ContingSim also suggests that I should
    %shift the note down instead of up. conting_3hits_92908;
    %092908_tmptest2 folder
    %12:14pm - back to original (4 hits)
    %12:23pm - upped it to 5hits
    %It's going to be interesting to see what happens with this bird
    %because the original (baseline) pattern of variation was to raise or
    %lower the entire note.
    %Getting the third harmonic looks very comparable to getting the 2nd
    %and 3rd harmonics, so I'm just going to make a white noise file that
    %is notched at 5000 Hz.
    %5pm-6pm - ran notched white noise to calibrate amplitude that gives a
    %good notch. --- chose 54
    %7pm - did a new white noise at lower amplitude and shorter time (30ms)
    %because the bird was stopping its song after the wn burst
    %'wn30low5000'
    %7:33pm - started at 2450Hz
    
    %Note that I chose to push the note upwards instead of downwards
    %because it appeared that the note was more flexible in this direction
    %on short timescales.  I concluded this by running a similation
    %(ContingSim) on some of the baseline data).
    
    %30;5;1e-6;2 are the segmenting specs for the wav data for the
    %deterioration of pk37
%Tuesday 093008
    %Today I segmented all of the previous song, finished virus updates,
    %formatted new hard drive, ordered new credit card, e-mailed Francesca.
    
    %I need to: feed bird, label pk37 deterioration, etc. ---transfer more
    %files
    
    %notch is okay but not great - but I'm keeping the amplitude at 54
    %12:45 - upped it to 2460Hz
    %Things look good! By setting the hit rate low, it appears that the
    %bird is now willing to sing.
    %There was a definite increase in pitch over the course of the morning.
    %5:25pm - upped it again to 2475Hz - unclear whether to raise or lower.
    %
%Wednesday 100108
    %10:35am - upped it again to 2495Hz - the note has gone up quickly in
    %the morning=great!
    %10:53am - shifted from 5consecutive hits to 6consecutive hits because
    %the bird was shifting the whole note up --- new targeting is later in
    %the note.
    %NOTE THAT all of 1001Afiles occur before the 5-6 switch to allow using
    %the first cntrng on them. 
    %Everything else uses cntrng2.
    %12:47pm - post-labmeeting -- song is rocketing upward - I changed
    %contingency to 2515Hz (keeping steady at 6consecutive hits).
    %2:15pm - upped it to 2525Hz
    %2:55pm - added a high template
    %3:12pm - shifted back to old template (and old folder) because too
    %many false positives - and the bird gets pissed off by white noise
%Thursday 100208
    %2:25pm - upped contingency to 2540Hz
    %3:25pm - brought it back down to 2535Hz (hit rate was a bit high)
    %3:35pm - added new template - it was the 90th point in the average
    %note of 1002Anotes - it was critical to set the second threshold at 3
    %instead of 4 to avoid false positives --- also upped to 2540Hz
    %3:40pm - upped it to 2550Hz
    %3:45pm - back to 2540Hz
    %3:55pm - up to 2545
    
%Friday 100308
    %10am - lowered contingency to 2530Hz
    %11am - increased white noise amplitude and made it un-notched
    %At some point in the mid-afternoon, he stopped singing as much, so I
    %reduced the amplitude.
%Saturday 100408
    %Kept the relatively low contingency, but there wasn't any sign that
    %would indicate ability to shift higher - so at night (just after 9pm),
    %I switched off the white noise.
    
%%%%%%% Experiment 2: Front_Up
%   Baseline: 1005c,1006a,1006b,1006c
%Monday 100608
    %At around 7pm it appeared that the mean pitch curve over the last four
    %segments (1005c through all of 1006) was the same as baseline and the
    %standard dev was also the same as baseline.  Thus I made a template
    %(templa3) and started running it on the same --- running on previous
    %song suggested a 1% false neg rate and around an 8% false positive
    %rate and a standard deviation of 5.6.  The template was made by taking
    %the average of points 85 to 90 in the average note from 1006Bfiles.
    %This targets the very beginning of the note (first threshold
    %crossing).
    %Based on simulation of the day's song, I decided to set the
    %contingency at 2440Hz.
    %Started running it at 8:30pm
%Tuesday 100708
    %11:00 - upped it to 2450Hz.  The original template is working well and
    %the original threshold is still fairly optimal. There were no false
    %negatives thus far in catch trials but song is quickly moving upwards.
    %11:35 - upped it to 2480Hz.
    %11:40 - decided 2470Hz was better - more conservative.
    %11:45 - Upped to 2480Hz
    %1:40pm - new template (tmp4) - same as old template but shifted
    %upwards (right) by one frequency bin.
    %
    %1:57pm - upped it to 2500Hz
    %4pm - realized that I should have both templates so bird can't escape
    %by going low - so I used OR logic on both templates === templa5
    
    %4:40 - discovered the problem - since I had been using OR logic the
    %bird could escape by going too low -- but correcting it allowed false
    %positives that caused false negatives due to their refractory period.
    %After correcting this, things seem better.
    %4:55 - bumped up to 2465
    %5:50 - bumped up to 2475 (hit rate had been 60%)
    %6pm - bumped up to 2485
    %8pm - new template (templa6 w/ cntrng3) - set to 2465 at beginning of 
    %note because avgnote appeared like the middle was being favored.  
    
%Wednesday 10/8/08
%A and B use templa6 and cntrng3
%C uses templa6 and cntrng6
    %10:20am - bumped up to 2475 because of low hit rate
    %10:40am - bumped further up to 2490 based on FFstats and low hit rate
    %10:50am - set to 2500
    %Rapidly trending upward - changed targeting to second hit because
    %it's hitting sometimes before the note starts.  Changed contingency to
    %2520Hz at 2pm.       templa6; cntrng6
    %Changed contingency to 2510.
    %Staying at 75% hit rate, not improving.
    %Contingency at 2500 overnight.
%Thursday 10/9/08
%A uses templa6 and cntrng6
%B uses templa8 and cntrng8
%C uses templa8 and cntrng8
    %10:40am - bumped up to 2510 based on low hit rate
    %11am - up to 2520 based on FFstats
    %11:50am - new template (tmp8) with cntrng8 - TH 2.4 for both, 2nd
    %crossing for both  --tmp8 is tmp6 plus a new, higher template
    %Anything before 11:50 is 1009A and uses cntrng6 and templa6.
    %Anything after is 1009B and uses cntrng8 and templa8
    %12:30pm - up to 2530 and right back down to 2520
    %2515 overnight
%Friday 10/10/08
    %A uses templa8 and cntrng8
    %B uses templa9 and cntrng9
    %C uses templa9 and cntrng9
    %
    %10am - bumped up to 2525 based on hit rate
    %11:13am (end of Afiles) - shifted to templa9 = combination of template
    %based on this morning and templa6. - better hit rate and targeting
    
%Sunday 10/12/08
%All data uses templa9 and cntrng9
%The data from the front shift experiment looks really good.




%%%% TSA experiment - October 12C until October 13

%New experiments starts with C (~4:27pm)
    %Started the new experiment at 4:27pm - I set the catch trial fraction
    %at 0.5 and the catch song fraction at 0.  I upped the contingency to
    %2540Hz (then immediately back down to 2535 at 4:28pm - after the first
    %8 detections).
%For the new experiment, I should use the program evtaf_freq to get the
%actual trigger values (neccessary for the non-catch trials).
    %At 5pm when I left, everything was working well.

%Monday 10/13/08
%2:18pm - song was remaining high - thus I was not driving specific
%learning, thus I wasn't tested what I wanted to test.  Based on the
%current version of evtaf_amp.vi, I couldn't figure out a way to do what I
%really want to do (that is, to HIT all notes below while hitting 20% of
%the notes above and compare the hits above to the escapes above).
%Instead, I decided to take shifted song and hit 50% of notes, thus
%removing any useful information and presumably causing new learning of the
%absence of association between white noise and neg/pos reinforcement).  I
%plan on analyzing the deterioration of song to see what happens.

%Tuesday 10/14/08
%TSA notcatch analysis - analyzed batch13pm (650 notes)
%TSA1014_3.m (does it all)
%sim.m (determines significance)
%saved as TSA_1014.mat in bk50w18 folder
dirf('*.cbin','batch')
evsonganaly
mk_tempf('batch1800',templa9,2,'obs0');
get_trigt2('batch1800',cntrng9,0.3,128,1,1);
vals1800=jctaf_freq('batch1800',[2200 2750],'a',128,'obs0',1,0);
%Appears that lots of recovery has occurred
%AllTSAdata.mat in bk50w18 folder has vals around 10am,1pm(1300)and 6pm(1800)

%Made figures - onsets and offsets - twice of minimum of SDbaseline.
%Subtracted 16 from mean to get the targdist
%Note - for bk50w18end I averaged the two peaks that were adjacent

%830pm --- removed the bird and put in pk37bk19 --- bk50 is now in Kris's
%room in soundbox 3
%838pm --- started collecting song



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Assay 3 - shifting the front of the note downwards.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% October 29 - placed in soundbox in Sam's room and began recording.
% October 30 and 31 - tested templates for targeting the front of the note
% October 31  - decided to target the first crossing of
% templa1down and began running the bird with a crappy speaker that didn't deliver a
% very high volume of white noise.
%   3:50pm - started with Evren's 80ms stimulus at 2440Hz

%   5:20pm - changed to my old 50ms stimulus at 2430Hz

% November 2 - came in and adjusted settings.  Preoccupied with file
% transfers and other things.
% November 4 (tuesday) at 1pm - got aggressive about shifting this bird.
% Changed to a good speaker with high volume white noise and decided that
% targeting was too early in the note.  Made a new template based on recent
% song (templa2down) and targeted the second template crossing.
% WHOOPS - I was hitting below - problem corrected at 6pm - set to 2400Hz.
% November 5 - 11am - set up to 2430Hz b/c hit rate was too low
%           - 5pm - set down to 2420Hz - looking better; targeting is good
% 
% November 8 - 12pm - very precise targeting and effective shifting. If
% this keeps up, I'm pretty much done.

