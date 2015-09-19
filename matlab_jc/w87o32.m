% w87o32

% MMAN lesion candidate - REPEATS
% segment at 500
% no need to cleandir
% label repeated note 'b'
% 6.03.11 - start screening
% 6.04.11
    [timevals0604,distnt0604] = jcrepeatdist('b','batchnotes'); 
    % mean=9.0, s.d.=2.06
% 6.07.11 - T.Dembo lesion

% 6.09.11 - female present from 4:15pm to 8:55pm
    % sang one song - at 4:50pm
% 6.11.11 - started singing a lot at dawn
    % No apparent changes to song
    % Can he learn???
    
% 6.12
% problem with refractory period on the repeat counter (0.1 is too short)
% 4:10pm - adjusted to 0.2

% 6.13 - WN day one - reduce repeat number (hit 8 and up)
    % template hits 7 and up (but often misses the first one)
% 6.15 - looks like starting to learn
% 6.16 - definite learning - wn off at dark
% 6.21 - still no recovery ?!
% 6.24 - finally recovered - indicates that it was learning

% 7.11 - did it really learn - TRY AGAIN with lower threshold
    % Do one more experiment - hit 5 and up (can it learn to only sing 3 or 4)?
        % adjust so TH say MIN=4, but that really means hitting 5 and up
        % because the template misses the first one
    % WN on at dawn on 7.15
    % clearly learned
    % recovered


figure;hold on;
plot(runningaverage(timevalsPre,20),runningaverage(distntPre,20),'b')
plot(runningaverage(timevalsWN,20),runningaverage(distntWN,20),'r')
plot(runningaverage(timevalsPost,20),runningaverage(distntPost,20),'b')
%
plot(runningaverage(timevalsPre2,20),runningaverage(distntPre2,20),'b')
plot(runningaverage(timevalsWN2,20),runningaverage(distntWN2,20),'r')
plot(runningaverage(timevalsPost2,20),runningaverage(distntPost2,20),'b')
