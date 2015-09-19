%%% Conclusion - odd - in both cases there are upshifts in pitch during wn
%%% playback, possibly due to notched WN???

%%% Experiment 1 - upshift %%%%
% 6/22/09 - baseline song - pitchBase
% 6/23/09 - 1:00pm chicktime - wnlowpass5000 hits all song
%                            - wnlowpass5000 hits bottom 70% of song from
%                               pitchBase dataset
% No playback overnight
% 6/24/09 - ~11am-3pm - lights out, playback remains on
% 6/24/09 - ~6pm - ampoff - i.e (WN off, playback off)

%%% Experiment 2 - downshift %%%%
% 6/25/09 - decay to baseline

% 6/26/09 - 11:50am realtime - wnlowpass6000 hits all song
%                            - wnlowpass6000 hits top 70% of song from
%                               pitch626ampoff dataset 
% 6/26/09 - 3pm - changed notch of real wn playback to 32000 to be correct
% 6/26/09 - 3:30pm - changed cntrng to hit later in the note to be more
% consistent with the playback version
% 6/27/09 - 1:35pm - lights off, amp on

pitchAll=[pitch623pm,pitch624A,pitch624B,pitch624C,pitch625ampoff,pitch626ampoff];
tvalsAll=[tvals623pm tvals624A tvals624B tvals624C tvals625ampoff,tvals626ampoff];
figure;plot(tvalsAll,median(pitchAll(800:860,:)),'*','Color','r')
hold on;plot(tvalsAll(1:end-172),median(pitchAll(800:860,1:end-172)),'*','Color','k')
plot(tvalsAll(1:136),median(pitchAll(800:860,1:136)),'*')

wc=20
for i=1:floor(size(pitchAll,2)/wc)
    start=i*wc-wc+1;
    last=start+wc-1;
    val1=[mean(median(pitchAll(800:860,start:last)))];
    plot([tvalsAll(start) tvalsAll(last)],[val1 val1],'Color','g','LineWidth',8)
end

figure;plot([tvals626ampoff tvals626on tvals627A],[median(pitch626ampoff(