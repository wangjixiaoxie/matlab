% pu87bk30
% Used as control for covert - see Covertpositivecontrols.m
        % 12.16.09 at 3:15pm - amp on
        %   WN hits the lowest 70% (3170Hz) - attempt to drive pitch up
        % 12.16.09 at dark - amp off
        % increased by ~30Hz (~70-80% escapes)
        % rapid decay the next morning - within an hour or so
  % LEARNING THE NEXT MORNING WAS NOT PRESENT OR NOT SIGNIFICANT      
        
        
        % 12.18.09 at 2:01pm - amp on
        %   WN hits the highest 70% (3100Hz) - attempt to drive pitch down
        % pre increases just before WN (using the end of the pre as baseline, it
        % looks like ~30Hz decrease and then rapid recovery on the next morning
        % 220:300
   % LEARNING THE NEXT MORNING IS UNCONVINCING
   
% surgery on 12.21.09 - bilateral RA cannulae implantation using new
    % coordinates - 10 degrees at 0.3mm caudal - localized very well
% singing in the morning of 12.22.09 - on computer in Sam's room

% good bird for syntax and syllable structure

% Implanted on 01.06.10 at ~5pm 
%       The first implants were quickly torn out by the bird.  Careful
%       investigation revealed that they were not damaged, thus I soaked
%       them in ethanol and then PBS and re-inserted, adding a small amount
%       of epoxy to the probe on the right side of the bird's head.  This
%       appeared to resolve the matter.  There was a slight amount of blood
%       and possibly also saline??? on the right side, which I cleaned out
%       by reinserting the place holders and then dabbing with a kimwipe.
% As of 01.07 at 2pm, the bird has chirped but not sang.

% 01.08.10 at 1:45pm - 4mM apv test run - 1.0uL/min
%          at 2:01pm - new template - hits all A's (i.e. high stack notes)
%                                   includes synshift
%                                   does note hit calls
%          at 6:48pm - acsf on at 1.5uL/min
% 01.09.10 at 11:21am - 4mM apv on - 1.0uL/min
%          at 12:25pm - wn on hit above 3100Hz - 0.8uL/min

% 1.10.10 - disk full at ~5pm
% 1.11 at ~10:30am - back on acquisition
%          11:45am - wn on all A transitions

% Syntax Stats
% 108acsf - 83.1% A
% 108apv - 74.5% A
% 109acsf - 81.3% A (small sample size)
% 109apv - 76.7% A
% 110acsf - 74.0% A
% 111am - 77.1% A

% 1.13 from 4-5:30pm --> 68% A (~70% A for the day)
% 1.13 at 6:15pm --> 4mMapv on at 1.0uL/min at 6:09pm
    % apv to ~71% - n.s.
% 1.13 at lights out (9:15pm) -acsf on
% 1.14 morning - 72.7% A
% 1.14 - probes out around 1pm - clogged

% Syntax quantification - WN on high stack (A) to increase p(noisy note B)
tv110acsfA=timing3(fv110acsfA);
tv110acsfB=timing3(fv110acsfB);

tvall110acsf=[tv110acsfA tv110acsfB];
numsall2=[ones(1,length([tv110acsfA])) zeros(1,length([tv110acsfB]))];
[b2,ind2]=sort(tvall110acsf);
numsorted110acsf=numsall2(ind2);
tvsorted110acsf=tvall110acsf(ind2);
% probability of B (escaping transition) as a function of time
figure;hold on;subplot(212);hold on;
ylim([0 0.5])
plot(runningaverage(tvsorted109acsf,40),1-runningaverage(numsorted109acsf,40),'*','Color','k')
plot(runningaverage(tvsorted109apv,40),1-runningaverage(numsorted109apv,40),'*','Color','b')
plot(runningaverage(tvsorted110acsf,40),1-runningaverage(numsorted110acsf,40),'*','Color','k')
% wn on
plot(runningaverage(tvsorted111acsf,40),1-runningaverage(numsorted111acsf,40),'*','Color','r')
plot(runningaverage(tvsorted113apv,40),1-runningaverage(numsorted113apv,40),'*','Color','g')
%%%%
% spectrogram
subplot(211);
imagesc(t,f,log(avZ));syn;ylim([0,1e4]);xlim([-1.3 0.1])
% variability reduction as metric of effectiveness - GREAT0.
mean(std(pitch113apvA(200:300,:)'))/mean(mean(pitch113apvA(200:300,:)'))
mean(std(pitch109apvA(200:300,:)'))/mean(mean(pitch109apvA(200:300,:)'))
mean(std(pitchPRE(200:300,:)'))/mean(mean(pitchPRE(200:300,:)'))
mean(std(pitch111acsfA(200:300,:)'))/mean(mean(pitch111acsfA(200:300,:)'))