% A1-->A2-->A3-->A4 ...
%   p1   p2   p3    ...

% transitionTW has 9 experiments, all of which have pre/dur/post except #4,
% which lacks post

% dur corresponds to the afternoon of the last day (usually 2nd day) of WN

% transitionTW(i).pstay is the important field
% the jth column of pstay corresponds to pj 
     % (i.e. the transition from the jth syllable to the j+1th syllable - see above)
% the rows of pstay correspond to pre, dur, post
