% pu90bk42
% Screened 05.15.11 - data on /cardinal4
cleandir4('batch',500,500,5,5);

% decent candidate for syntax learning
% very solid variability, possible to target 

% (A)(B)(C)DEFFFF
% (A)(B)CDEGHI
% (A)(B)CDGHI
% A is rare
% B is about 50/50
% Most common are BCDEFFFF vs. BCDEGHI vs. BCDGHI
% CDEFFFF vs. CDEGHJ vs. CDGHJ
% Target CDEFFFF vs. CDEGHJ
%          71%         29%
% Seems targetable - should be pretty easy to hit the first F - worth screening
% Also has targetable high stack notes

%  6.1.11 - transferred to launchpad sound box for screening/learning
%  6.2.11 - 5:40pm - changed filter settings - segment at 500
    labels=notestats
    length(findstr(labels','cdef')) % 56 - 72%
    length(findstr(labels','cdeg')) % 22 - 28%
% 6.3.11 - 1:45pm - loaded template - but had the third tmp Evtaf instead
    % of BirdTAF
% 6.4.11 - 12:45pm - corrected template

%******** Data (tmp, etc) on /cardinal9/SyntaxScreening/pu90bk42 *****%