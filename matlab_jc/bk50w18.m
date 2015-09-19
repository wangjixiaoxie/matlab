bk50w18
%
% 100 percent hit rate
% WN on at dawn on March 23, 2010
% Starting singing frequently at dawn on March 24th

% Shifted up - significantly as of the 25th
% 1. at 8:45am computer time, there were escapes
% 2. possible effect of notched WN
% 3. pushing the WN back along the note

    % 11:20am - new template (more permissive)
    % 11:40am - new WN file (non-notched)
% Looks like these changes were successful

% March 26, 2010 - WN off at 12:07pm

% March 28, 2010 - WN on at dawn

% March 29, 2010 - still looks good

% March 31, 2010 
            % After being very stable, pitch began to decrease
         % AT ~ 9pm (after 4 days of 100%, changed to 90% of pvs day):
            % Hit below 2470Hz, this threshold was fairly low because of decrease 
            % in pitch during the final day of 100% WN
% April 1, 2010 - noon - hit below 2490Hz (around 90% of pvs four days)

figure;hold on;
plot(tvals100prct,pitch100prct(625,:),'*')
plot(tvals90prct,pitch90prct(625,:),'*','Color','r')
plot(tvals322preA,pitch322preA(625,:),'*','Color','k')
plot(tvals327postA,pitch327postA(625,:),'*','Color','k')

