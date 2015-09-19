% LesionAnalysisRotation
% Look through /swift 5, /swift6, /swift7, /swift8, /swift9
% In each hard drive partition (e.g. /swift5), there are a bunch of folders
        % with birds' names (e.g. 'bk72w64') each containing a bunch of folders
        % containing song files. 
% For each bird, find all complete days (around 7am to around 8pm, or at
    % least 9am to 7pm), that are not within 48hrs after the end of an ampon
    % session (check this by looking at time stamps in ampon folders).  
    % Look at prelesion data and postlesion data (which is in an explicitly
    % labeled 'postlesion' folder)
    % I have done this step for two of the birds (see below).
% Load data files with frequency and time of day information for pvsly
% labeled song.
        load FREQDATA.mat % Should be in birds' name direction but may be a separate file in /postlesion directory in some cases
        figure;plot(vals(:,1),vals(:,2),'*') % plots time (minutes) on x-axis and frequency (Hz) on y-axis
        % Note the gaps in between days 
        datestr(vals(i,1)) % allows you to convert from minutes (from arbitrary starting point)
                                % to day of the year (to align freq data with pvsly
                                % identified complete days of interest)
% 1. Find the complete days--- DONE
% 2. For each of these dates (e.g - May 5) how many days after beginning of recording?
%         datestr(LESIONDAY) % e.g May 17
%         figure;plot(vals(:,1)-LESIONDAY,vals(:,2),'*') % plots time (days) on x-axis and frequency (Hz) on y-axis
    % e.g. for bk13w63, day 1 is n=-15.2, or floor(n)=16 days before the lesion day (May 17th) - so it is May 1st.
    % So now you have enough information to determine which dates
    % correspond to which vals.
    % check with: datestr(vals(i,1))
% 2. (explicitly)
% Segment vals into days and stick data into an array
%figure;plot(diff(sort(vals(:,1))-LESIONDAY)*24) % look where this is > 10
[b,sortedtvals]=sort(vals(:,1));
dayends=find(diff(b-LESIONDAY)*24>10); % tells you the last point in each day
timesorted=vals(sortedtvals,1);
FFsorted=vals(sortedtvals,2);
bk13w63.indices=[]; %
daystart=1;
for i=1:length(dayends)
    bk13w63(i).indices=[daystart:1:dayends(i)];
    daystart=dayends(i)+1;
    bk13w63(i).timevals=timesorted(bk13w63(i).indices);
    bk13w63(i).FFvals=FFsorted(bk13w63(i).indices);
end
figure;hold on; % Make sure this looks the same as the original
for i=1:length(dayends)
    plot(bk13w63(i).timevals,bk13w63(i).FFvals,'*')
end

% Now all you need to do is determine which days you care about
    % For each value of m, is this one of the days I'm interested in?
    % Look at a few values of n ---> e.g. 2, 11, 20
        datestr(vals(bk13w63(m).indices(n),1));
 
% Now all you need to do is determine which days you care about
% FOR EACH m
    % For each value of m, is this one of the days I'm interested in?
    % Look at a few values of n to verify ---> e.g. 2, 11, 20
        datestr(vals(bk13w63(m).indices(n),1));
    % Then look at the first value of the day
        datestr(vals(bk13w63(m).indices(1),1));
       firsttime=2    % number of minutes the first value of the day is after 7:00 am
       firsttime=bk13w63(m).timevals(1)*60*24-firsttime;
       for i=1:length(bk13w63(m).timevals)
           bk13w63(m).timevals(i)=bk13w63(m).timevals(i)*60*24-firsttime;
       end
        
        
    % And...make a new field:
    bk13w63(m).keep=1; % if you want to keep it or
    bk13w63(m).keep=0; % if you want to ignore it
% END

               
% If any of the non-ignored days have less than 20 notes labeled, you will
% need to find the folder corresponding to those days, label more song, and
% calculate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SWIFT5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%
% bk13w63
%%%%%%%%%
    % pre-lesion days (May 2-5) can be found in
        /swift5/bk13w63/screen
        /swift5/bk13w63/trigtest
        /swift5/bk13w63/trigtest2
    % post-lesion days (May 21-30 and June 21-22) can be found in
        /swift5/bk13w63/postlesion/screen
        /swift5/bk13w63/postlesion/trigtest
        /swift5/bk13w63/postlesion/trigtest2
        /swift5/bk13w63/postlesion/ampoff2
%%%%%%%%%
% bk72w64
%%%%%%%%%
    % pre-lesion days (June 14-18) can be found in
        /swift5/bk72w64/screen
        /swift5/bk72w64/trigtest
    % post-lesion days (June 28-July 7 and July 18-21 and August 3-11 and August 30-September 5) can be found in
        /swift5/bk72w64/postlesion/screen % note that song files for last days of June are in alpha order after first song files for July
        /swift5/bk72w64/postlesion/trigtest
        /swift5/bk72w64/postlesion/trigtest2
        /swift5/bk72w64/postlesion2/trigtest       
        /swift5/bk72w64/postlesion2/trigtest2        
        /swift5/bk72w64/postlesion2/ampoff        
%%%%%%%%%
% bk35w47
%%%%%%%%%
    % pre-lesion days (ampon: Dec14-18)
        /swift5/bk35w47/trigtest            % Dec11-13
        /swift5/bk35w47/trigtest2           % finishes Dec13
        /swift5/bk35w47/ampoff              % Dec21-28
    % post-lesion days (ampon: Jan13-18, Jan25-30, Feb11-15)
        /swift5/bk35w47/postlesion/screen               % Dec31, jan1, begin Jan2 
        /swift5/bk35w47/postlesion/screen2              % end Jan2, Jan3-4, Jan5-7, Jan9
        /swift5/bk35w47/postlesion/trigtest             % end Jan9, Jan10
        /swift5/bk35w47/postlesion/trigtest2            % end Jan10, begin Jan11
        /swift5/bk35w47/postlesion/trigtest3            % end Jan11, Jan12
        /swift5/bk35w47/postlesion/trigtest4            % Jan22, begin jan23...
        /swift5/bk35w47/postlesion/trigtest5            % Jan24
        /swift5/bk35w47/postlesion/trigtest_rep         % Feb5-7, begin Feb8
        /swift5/bk35w47/postlesion/trigtest_rep_shortpb % end Feb8, Feb9-10
        /swift5/bk35w47/postlesion/ampoff               % Jan21, begin Jan22
        /swift5/bk35w47/postlesion/ampoff2              % Jan31, Feb1-3
        
        
%%%%%%%%%
% bk24w14 NO LESION DATA
%%%%%%%%%

%%%%%%%%%
% o50bk72
%%%%%%%%%
    % pre-lesion days (ampon: Nov4-10)
        /swift5/o50bk72/pretrain                    % Oct30-31, Nov1-3
        /swift5/o50bk72/ampoff                      % Nov13-16
    % post-lesion days (ampon: Dec2-10, Jan8-14)
        /swift5/o50bk72/postlesion/ampoff           % Dec13-15
        /swift5/o50bk72/postlesion/ampoff2          % Jan15-25
        /swift5/o50bk72/postlesion/morescreen       % Jan3-4, Dec17-20
        /swift5/o50bk72/postlesion/screen           % Dec1, Nov21, Nov24-30
        /swift5/o50bk72/postlesion/trigtest         % Jan5-7 
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SWIFT6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       

%%%%%%%%%
% pk20bk78
%%%%%%%%%
    % pre-lesion days (ampon: Feb14-21)
        /swift6/pk20bk78/ampoff                     % Feb22-28, Mar1-13
        /swift6/pk20bk78/pretrain                   % Feb8-12
    % post-lesion days (ampon: Mar26-Apr4)
        /swift6/pk20bk78/afterlesion/ampoff         % Apr7-10
        /swift6/pk20bk78/afterlesion/screen         % Mar17-21
        /swift6/pk20bk78/afterlesion/trigtest       % end Mar21, Mar22-25
        
%%%%%%%%%
% pk208k66
%%%%%%%%%
    % pre-lesion days (ampon: Sept18-22)
        /swift6/pk28bk66/ampoff                     % Sept25
        /swift6/pk28bk66/ampoff_b                   % Sept26-29, Oct1-2
        /swift6/pk28bk66/screen                     % Sept14
        /swift6/pk28bk66/trigtest                   % end Sept14, Sept15-17
        /swift6/pk28bk66/newscreen                  % Oct22-31, Nov1-2
        /swift6/pk28bk66/newscreen2                 % Nov4-6,
        % post-lesion days (ampon: Dec5-10,Dec22-27)
        /swift6/pk28bk66/postlesion/ampoff          % Dec13-21 
        /swift6/pk28bk66/postlesion/ampoff2         % Dec30-31, Jan1
        /swift6/pk28bk66/postlesion/screen          % Nov10-11, Dec1-2, begin Dec3
        /swift6/pk28bk66/postlesion/trigtest        % end Dec3 
        /swift6/pk28bk66/postlesion/trigtest2       % end Dec3, Dec14
      
%%%%%%%%%
% pk29bk36
%%%%%%%%%
    % pre-lesion days (ampon: Feb10-19)
        /swift6/pk29bk36/ampoff                     % Feb22-28
        /swift6/pk29bk36/pretrain                   % Jan26-31, Feb1-3, Feb6-9
    % post-lesion days (ampon:Mar21-29, Apr11-17)
        /swift6/pk29bk36/afterlesion/ampoff         % Apr1-10
        /swift6/pk29bk36/afterlesion/ampoff2        % Apr17-23
        /swift6/pk29bk36/afterlesion/screen         % Mar4-13
        /swift6/pk29bk36/afterlesion/trigtest       % Mar14-20
        /swift6/pk29bk36/afterlesion/trigtest2      % Apr10
    
%%%%%%%%%
% pk30bk79
%%%%%%%%%
    % pre-lesion days (ampon: Dec13-24)
        /swift6/pk30bk79/prelesion/ampoff           % Dec25-Jan4
        /swift6/pk30bk79/prelesion/screen           % Dec9-12
    % post-lesion days (ampon: Jan29-Feb3)
        /swift6/pk30bk79/postlesion/ampoff          % Feb6-7
        /swift6/pk30bk79/postlesion/ampoff2         % Feb17-22
        /swift6/pk30bk79/postlesion/ampoff3         % Feb24-26
        /swift6/pk30bk79/postlesion/screen          % Jan15-21
        /swift6/pk30bk79/postlesion/trigtest        % Jan23
        /swift6/pk30bk79/postlesion/trigtest2       % Jan25-28
        /swift6/pk30bk79/postlesion/trigtest3       % Feb9-10
        
%%%%%%%%%
% pk42bk73
%%%%%%%%%
    % pre-lesion days (ampon: Jun15-17)
        /swift6/pk42bk73/ampoff                     % Jun20-29 
        /swift6/pk42bk73/screen                     % Jun13
        /swift6/pk42bk73/trigtest1                  % Jun14
        /swift6/pk42bk73/trigtest2                  % end Jun14
    % post-lesion days (ampon:)
        /swift6/pk42bk73/afterlesion/ampoff         % July14-20
        /swift6/pk42bk73/afterlesion/screen         % July21-26
        /swift6/pk42bk73/afterlesion/trigtest2      % July 13 
   
        
%%%%%%%%%
% pk80r6
%%%%%%%%%
    % pre-lesion days (ampon: Nov11-16)
        /swift6/pk80r6/prelesion/wnoff1             % Nov19
        /swift6/pk80r6/prelesion/wnoff2-b           % Dec10-13
        /swift6/pk80r6/prelesion/screen2            % Nov8
        /swift6/pk80r6/prelesion/temptest3          % Nov10
    % post-lesion days (ampon: Dec28-Jan2)
        /swift6/pk30bk79/postlesion/ampoff          % Jan5-13
        /swift6/pk30bk79/postlesion/screen          % Dec17-20, Dec23-24
        /swift6/pk30bk79/postlesion/morescreen      % Jan14
        /swift6/pk30bk79/postlesion/trigtest        % Dec27
        /swift6/pk30bk79/postlesion/trigtest2       % end Dec27
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SWIFT7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%%
% o7bk77
%%%%%%%%%
    % pre-lesion days (ampon: ?)
        /swift6/o7bk77/prelesion/screen             % 
    % post-lesion days (ampon: )
        /swift6/o7bk77/afterlesion/ampoff           % 