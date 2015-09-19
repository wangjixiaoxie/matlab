% bk80w28

% September 2011 - not a good directed singer

% Did something (hit/miss) with him in August 2009 - see /cardinal
% September 2010 -  NatNeuro ctls

% 11.16.10 - bilateral LMAN stim array implantation
% 11.17.10 - song looks good - decent singing rate
cleandir4('batch',50000,500,6,10);

% 11.23.10 - make stim template

%**************************%
% Created new stim rig - using one stimulator with current split (shorted)
%**************************%
% Stim template recognizes noisy note that precedes the low stack note
% Stim in premotor window of low stack note
% 12.13.10 - testing stim configurations
    % 11:20am --- 20ms_delay, 80ms_duration, R14, L23, 80uA - seems to reduce s.d. w/o much offset
    % 12:30pm --- 20ms_delay, 80ms_duration, R23, L14, 80uA - seems to not reduce s.d., possible offset
    % 3:00pm --- 20ms_delay, 80ms_duration, R14, L23, 80uA 
    % 9:00pm --- 0ms_delay, 80ms_duration, R14, L23, 100uA -
% 12.14.10
    % 12pm - 0ms delay, 80ms duration, R13, L23, 100uA - upward deflection
    % 4pm - 0ms delay, 80ms duration, R23, L23, 100uA - upward deflection, great s.d. reduction
    % 7pm - 0ms delay, 80ms duration, R23, L34, 100uA - upward deflection, okay s.d. reduction
% 12.15.10
    % 10:40am - 0ms delay, 80ms duration, R23, L23, 100uA - seeemd to be best at reducing s.d.
    % 12:05pm - hit below 2310Hz(median of non-stimmed trials +25Hz) - 10% stim
    % 6:45pm - hit below 2325Hz (should have done this originally- poor initial estimate of median)
% 12.16.10 - thursday
    % good learning, partial reversion
    % 12:05pm - wn off, stim 50%
    % 2pm - stim 20%, higher uA (150uA)
% 12.17.10
    % 12:36pm - new stim template - better precision]
% 12.18.10 (Experiment(12:15)
    % 100% experiment completely failed to block learning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  --- START with Experiment #16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% using rig with two stimulators 
% 1.05.11 - screening
% 1.06.11 - 10am - 80uA, R23, L14, 0ms delay, 80ms duration
                % Big upward offset, especially later in note, no s.d. reduction
%         - 3:15pm - switch to R23, L23 
% 1.07.11 - 10:27am - same settings (R23,L23,0msdel,80msdur,80uA) - 100% stim
%   % median as 2310Hz (might be slightly higher)
    % hit above 2310-50/2=2285Hz
    % - 12:20pm - WN ON - hit above 2285Hz, no pic chip

% 1.08.11 
%   - 10:30am - 100% stim - still not singing much
%   - 1:13pm - 2320Hz+50/2=2345Hz - hit below 2345Hz, no pic chip

% 1.09.11 - not singing much - note that stimulators produce small deltaV
% 1.10.11
%   - 10:10am - wn off - failed to block learning, despite large reduction in singing rate
%   - tried a few configs and increased voltage,  problem with singing rate
%   reduction 
% - increased to 200uA at R23,L23 (what I had been using) had only subtle effects suggesting not near LMAN
% - WHEREAS increasing to 200uA at R24,L14 caused major disruption confirming placement near LMAN
% - big reduction in singing rate but promising chocie of stim wires
% 1.11.11 
%  5pm - start 20% stim R23,L23 at 20prct - 100uA
% 1.12.11
%  9pm - start WN on, hit above 2310Hz (20th prctile) - shift down 
% 1.13.11 (DAY ONE - intermittent stim) - 100uA
% 12:20pm - hit above 2280Hz
% 1.14.11
% 10am - now at 2245Hz
% 1.15.11
% noon - lots of learning, not much reversion - increase to 200uA per side
% switch configs if that doesn't work
% 1.16.11
% the 200uA causes a great (50Hz) reversion
% 11:10am - hit above 2185Hz
% ~6:30pm - hit above 2160Hz
% 1.17.11
% lots of learning and consolidation has occurred
% noon - stim effect still looks good               
% noon - hit above 2130Hz

% 1.18.11 - nightfall - wn off

% 1.25.11 (DAY ONE - intermittent stim) 
    % 200uA on each of the stimulators
    % ~0.75 Volt offset has remained consistent during pvs 10 days
    % dawn - hit above 2250Hz (~25th percentile) - shift down
    % 5pm - hit above 2230Hz
% 1.26.11 (DAY TWO) - reversion looks great!
    % 11am - hit above 2110Hz
    % 7pm - hit above 2180Hz
% 1.27.11 (DAY THREE)
    % 11am - hit above 2160Hz
% 1.28.11 (DAY FOUR)
    % noon - hit above 2140Hz
    % 5pm - voltage looks good - 0.75V on right, 0.7V on left
% 1.29.11
    % 7pm - hit above 2110Hz
% 1.30.11
    % 11:30am - hit above 2120Hz (maintain this - asymptote)

% 2.01.11 - 9pm WN OFF

% 2.08.11 - 10:30am - stim 100% - note that voltage is less (0.5V on right,0.6V on left)
%         - 12:10pm - stim 100%, wn on, hit below 2370Hz (since stim causes upward offset)
% problem - oscilloscope tester left attached - experiment stopped at 2pm

% 2.09.11 
    % 9:30am - 250uV
    % 11am - stim 100% 
    % 12:30pm - stim 100%, wn on, hit below 2370Hz
% 2.10.11
    % learning occurred with stim 100% - no block of expression
% 2pm - leads out - transferred to different soundbox
%


% Did something (hit/miss) with him in August 2009 - see /cardinal
% September 2010
% 09.03 - put in launchpad soundbox - for NatNeuro ctls
%   wntest from 12-2pm - never sang through but good initiation of song

% 09.07 - wn test for notch - too loud
%   template test - first threshold crossing
%   hit below 2285Hz at 8pm
% 09.08 
%   hit below 2300Hz at 10:06am
%   hit below 2320Hz at 12:55pm
%   hit below 2335Hz at 3:25pm
%   hit below 2350Hz at 4:20pm
%   hit below 2365Hz at 5pm
%   hit below 2380Hz at 8:15pm
%  looks like he learned a lot in the morning and then it has leveled off

% 09.09
%   hit below 2420Hz
% 09.10
%   WN off at 12:25pm - allow recovery and then repeat the same experiment

% 09.12 - median for entire day is 2300Hz
%   10:30pm - wn on, hit below 2300Hz
% 09.13
%   10:38am - hit below 2330Hz
%   10:46am - hit below 2380Hz
%   11:05am - adjusted template slightly to recognize BTAF better
%   12:23pm - hit below 2400Hz
% WN OFF at dark

% 09.15
%   dark - hit below 2320Hz
% 09.16
%   hit below 2350Hz
%   4:25pm - increased WN volume
%   9:00pm - decreased WN volume - too loud for notch to work
%  hit below 2380Hz - keep this as the "asymptote" over the weekend
% 09.17

fvBaseline=findwnoteJC('batchnotes','b','','',0,[6000 8100],8500,1,'obs0',0);
for i=1:length(fvBaseline)
    shiftedBaseline(i,:)=fvBaseline(i).datt;
end
pitchBaseline=jc_pitchmat1024(shiftedBaseline,1024,1020,1,6000,8100,[1],'obs0',1);
fvWNon=findwnoteNN('batchnotes','b','','',0,[6000 8100],8500,1,'obs0',0);
        recognized=[];
        hitrange=[];
        catchtrial=[];
        for i=1:length(fvWNon)
            shiftedWNon(i,:)=fvWNon(i).datt;
            recognized(i)=fvWNon(i).CATCH>-1;
            hitrange(i)=fvWNon(i).TEMPLATE==1;
            catchtrial(i)=fvWNon(i).CATCH==1;
        end
        indEscapeAbove=recognized & hitrange & catchtrial;
        indHitAbove=recognized & hitrange & ~catchtrial;
pitchWNon=jc_pitchmat1024(shiftedWNon,1024,1020,1,6000,8100,[1],'obs0',1);    
figure;plot(median(pitchWNon(:,indEscapeAbove)'))
hold on;plot(median(pitchWNon(:,indHitAbove)'),'r')

for i=1:size(pitchWNon,2)
    residWNon(:,i)=pitchWNon(:,i)-mean(pitchBaseline')';
end
b=zeros(1363,1363);
for i=1:1363
    [b(i,1:i),bint(i).data]=regress(residWNon(:,i),residWNon(:,1:i));
    i
end    
regcoefs=b;
clear mb
for i=1:length(regcoefs)
    mb(i)=mean(regcoefs(i+1:length(regcoefs),i));
end



% load /cardinal3/bk80w28/Analyzed910.mat
    % Beautiful!!!!!
        figure;plot(runningaverage(t(indEscapeAbove),50),runningaverage(mb(indEscapeAbove),50),'r')
        hold on;plot(runningaverage(t(indHitAbove),50),runningaverage(mb(indHitAbove),50),'b')












fvals80w7Ahit=findwnoteJC('batch807Anotes','a','','',0,[2000 2700],8500,1,'obs0');
fvals807Amiss=findwnoteJC('batch807Anotes','b','','',0,[2000 2700],8500,1,'obs0');
for i=1:length(fvals807Ahit)
shifted807Ahit(i,:)=fvals807Ahit(i).datt;
end
for i=1:length(fvals807Amiss)
shifted807Amiss(i,:)=fvals807Amiss(i).datt;
end
pitch807Ahit=jc_pitchmat1024(shifted807Ahit,1024,1020,1,2000,2700,[1],'obs0',1);
pitch807Amiss=jc_pitchmat1024(shifted807Amiss,1024,1020,1,2000,2700,[1],'obs0',1);
pitch807A=[pitch807Ahit,pitch807Amiss];
figure;plot(mean(pitchBase'))
hold on;plot(mean(pitch807A'),'r')
hold on;plot(mean(pitch807Amiss'),'g')
hold on;plot(mean(pitch807Ahit'),'k')
