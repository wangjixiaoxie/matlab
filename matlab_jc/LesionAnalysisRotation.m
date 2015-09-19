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
%         datestr(vals(i,1)) % allows you to convert from minutes (from arbitrary starting point)
                                % to day of the year (to align freq data with pvsly
                                % identified complete days of interest)
% 1. Find the complete days--- DONE
% 2. For each of these dates (e.g - May 5) how many days after beginning of recording?
       datestr(LESIONDAY) % e.g May 17
       figure;plot(vals(:,1)-LESIONDAY,vals(:,2),'*') % plots time (days) on x-axis and frequency (Hz) on y-axis
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
birdname.indices=[]; %
daystart=1;
for i=1:length(dayends)
    birdname(i).indices=[daystart:1:dayends(i)];
    daystart=dayends(i)+1;
    birdname(i).timevals=timesorted(birdname(i).indices);
    birdname(i).FFvals=FFsorted(birdname(i).indices);
end
figure;hold on; % Make sure this looks the same as the original
for i=1:length(dayends)
    plot(birdname(i).timevals,birdname(i).FFvals,'*')
end

% Now all you need to do is determine which days you care about
    % For each value of m, is this one of the days I'm interested in?
    % Look at a few values of n ---> e.g. 2, 11, 20
%         datestr(vals(birdname(m).indices(n),1));
 
% Now all you need to do is determine which days you care about
% FOR EACH m
    % For each value of m, is this one of the days I'm interested in?
    % Look at a few values of n to verify ---> e.g. 2, 11, 20
%         datestr(vals(birdname(m).indices(n),1));
    % Then look at the first value of the day
%         datestr(vals(birdname(m).indices(1),1));

%%%%%%%%% changes for each bird %manually find this for all days you care about.
% keepmatrix=[x y z]; %indices of days you want to keep
% firsttimematrix = [xx yy zz]; % corresponding to |keepmatrix|
                %% number of minutes the first value of the day is after 7:00 am 


% when is the lesion day? 
% datestr(LESIONDAY)


 % spits out ALL days and their corresponding index. go through manually and setup the |keepmatrix|  
 sortedvalues=vals(sortedtvals,:);
 for  m = 1:length(birdname)
     datestr(sortedvalues(birdname(m).indices(1),1)) 
     m % use to make |keepmatrix|
 end

%%%%%%%%%%%%%% 
% keepmatrix = [];

% after figuring out |keepmatrix| - create the |firsttimematrix| easily
 for  k = keepmatrix
     datestr(sortedvalues(birdname(k).indices(1),1)) 
     k % use to make |firsttimematrix|
 end

%%%%%%%%%%%%%%
% firsttimematrix = [];


% double check that dates and minutes matrice are same length 
length(keepmatrix)
length(firsttimematrix)


%% pull out dates you care about and change the timevals to "minutes after 7am" so can analyze later
for m = 1:length(birdname)
    birdname(m).keep=0; % default
    if ~isempty(find(keepmatrix==m)) % 0 or 1 depending on value of keepmatrix
        birdname(m).keep=1; %change default to 1
        firsttime=birdname(m).timevals(1)*60*24-firsttimematrix(find(keepmatrix==m));
        
        % for all kept days, how many days before/after lesion. negative
        % values = pre-lesion. positive values = postlesion
        % days are rounded to nearest integer
        birdname(m).dayspostlesion = round(sortedvalues(birdname(m).indices(1),1) - LESIONDAY);
        
        
        for i=1:length(birdname(m).timevals)
           birdname(m).timevals(i)=birdname(m).timevals(i)*60*24-firsttime; % change time vals to "minutes after 7am"
        end
    end
end

%% at end - change |birdname| to actual birdname then save as .mat
o7bk77  = birdname;
save o7bk77ANALYZED o7bk77

% If any of the non-ignored days have less than 20 notes labeled, you will
% need to find the folder corresponding to those days, label more song, and
% calculate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SWIFT5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%
% bk13w63 ** DONE
%%%%%%%%%
% keepmatrix=[2 3 5 19:27]; 
% firsttimematrix = [0 4 2 25 10 7 5 19 3 5 8 3]; 

    %% pre-lesion days 
    % May 2-3, 5 [2 3 5]
    % Unlabeled: May 4
        /swift5/bk13w63/screen
        /swift5/bk13w63/trigtest
        /swift5/bk13w63/trigtest2
    
    %% post-lesion days 
    % May 21-25, May27-30 [19 27] 
    % Unlabeled: June 21-22 
        /swift5/bk13w63/postlesion/screen
        /swift5/bk13w63/postlesion/trigtest
        /swift5/bk13w63/postlesion/trigtest2
        /swift5/bk13w63/postlesion/ampoff2

        
%%%%%%%%%
% bk72w64 ** DONE
%%%%%%%%%

% keepmatrix=[1:5 14 16:21 32:36 38:40];
% firsttimematrix = [1 32 8 1 60 126 10 53 36 163 48 93 6 190 69 36 -5 48 -6 -4];

    %% pre-lesion days 
    %June 14-18 [1:5]
        /swift5/bk72w64/screen
        /swift5/bk72w64/trigtest
        
    %% post-lesion days 
    % June 28-30 [14-16]    
    % July1-2 [17-18]
    % July5-7 [19-21]
    % July18-20 [32-34]
    % Aug6-7 [35-36]
    % Aug 9-11 [38-40] 
    % Not Labeled: July3-4, July21, Aug3-5, Aug8, Aug30-Sept5
        /swift5/bk72w64/postlesion/screen % note that song files for last days of June are in alpha order after first song files for July
        /swift5/bk72w64/postlesion/trigtest
        /swift5/bk72w64/postlesion/trigtest2
        /swift5/bk72w64/postlesion2/trigtest       
        /swift5/bk72w64/postlesion2/trigtest2        
        /swift5/bk72w64/postlesion2/ampoff       
        
        
        
%%%%%%%%%
% bk35w47 ** DONE
%%%%%%%%%

% keepmatrix=[2:4 12:18 20:23 32:35 41:43 45];
% firsttimematrix = [167 135 136 133 181 128 95 57 38 29 65 28 36 39 22 23 18 170 10 58 81 19];

    %% pre-lesion days (ampon: Dec14-18)
        /swift5/bk35w47/trigtest            % Dec11-13
        /swift5/bk35w47/trigtest2           % finishes Dec13
        /swift5/bk35w47/ampoff              % Dec21-28
    % Dec11-13 [2:4]
    % Dec21 [12]
    % Dec 27-28 [13:14]
    % Unlabeled: Dec22-26
    
    % post-lesion days (ampon: Jan13-18, Jan25-30, Feb11-15)
        /swift5/bk35w47/postlesion/screen               % Dec31, jan1, begin Jan2 
        /swift5/bk35w47/postlesion/screen2              % end Jan2, Jan3-7, Jan9
        /swift5/bk35w47/postlesion/trigtest             % end Jan9, Jan10
        /swift5/bk35w47/postlesion/trigtest2            % end Jan10, begin Jan11
        /swift5/bk35w47/postlesion/trigtest3            % end Jan11, Jan12
        /swift5/bk35w47/postlesion/trigtest4            % Jan22, begin jan23...
        /swift5/bk35w47/postlesion/trigtest5            % Jan24
        /swift5/bk35w47/postlesion/trigtest_rep         % Feb5-7, begin Feb8
        /swift5/bk35w47/postlesion/trigtest_rep_shortpb % end Feb8, Feb9-10
        /swift5/bk35w47/postlesion/ampoff               % Jan21, begin Jan22
        /swift5/bk35w47/postlesion/ampoff2              % Jan31-Feb3
    % Dec31-Jan3 [15:18]
    % Jan9-12 [20:23]
    % Jan21-24 [32:35]
    % Jan30-Feb1 [41:43]
    % Feb5 [45]
    % Unlabeled: Jan4-7, Feb2, Feb6-10
        
        
%%%%%%%%%
% bk24w14 NO LESION DATA
%%%%%%%%%

%%%%%%%%%
% o50bk72 ** DONE
%%%%%%%%%

% keepmatrix = [1:5 15:18 22 25:32 44:46 48:50 52:56];
% firsttimematrix = [-38 -40 -38 -39 4 1 1 -2 -1 52 20 10 2 7 2 0 3 6 -1 0 237 -2 1 2 3 1 0 -2 -1];

    %% pre-lesion days (ampon: Nov4-10)
        /swift5/o50bk72/pretrain                    % Oct30-Nov3
        /swift5/o50bk72/ampoff                      % Nov13-16
    % Oct30-Nov3 [1:5]
    % Nov13-16 [15:18]
    % Unlabeled:
    
    %% post-lesion days (ampon: Dec2-10, Jan8-14)
        /swift5/o50bk72/postlesion/ampoff           % Dec13-15
        /swift5/o50bk72/postlesion/ampoff2          % Jan15-25
        /swift5/o50bk72/postlesion/morescreen       % Jan3-4, Dec17-20
        /swift5/o50bk72/postlesion/screen           % Dec1, Nov21, Nov24-30
        /swift5/o50bk72/postlesion/trigtest         % Jan5-7 
    % Nov21 [22]
    % Nov24-Dec1 [25:32]
    % Dec13-15 [44:46]
    % Dec17-19 [48:50]
    % Jan3-7 [52:56]
    % Unlabeled: Dec20 
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SWIFT6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       

%%%%%%%%%
% pk20bk78 ** DONE
%%%%%%%%%

% keepmatrix = [1:5 15:26 30:38 51:53];
% firsttimematrix = [12 10 26 16 14 14 23 62 0 20 12 13 27 17 20 -1 24 82 52 115 58 30 131 0 0 0 42 114 23];

    %% pre-lesion days (ampon: Feb14-21)
        /swift6/pk20bk78/ampoff                     % Feb22-28, Mar1-13
        /swift6/pk20bk78/pretrain                   % Feb8-12
    % Feb8-12 [1-5]
    % Feb22-28 [15-21]
    % Mar1-2, 5-6, 11 [22-26]
    % Unlabeled: Mar3-4, Mar7-10, Mar12-13
        
    %% post-lesion days (ampon: Mar26-Apr4)
        /swift6/pk20bk78/afterlesion/ampoff         % Apr7-10
        /swift6/pk20bk78/afterlesion/screen         % Mar17-21
        /swift6/pk20bk78/afterlesion/trigtest       % end Mar21, Mar22-25
    % Mar17-25 [30-38]
    % Apr7-9 [51-53]
    % Unlabeled: Apr10

    
%%%%%%%%%
% pk28k66 ** DONE  
%%%%%%%%%
    
% keepmatrix = [2:5 12:15 16:18 20 21 24:26 36:39 40 49:50];
% firsttimematrix = [27 32 160 26 101 86 57 10 -43 -53 -53 32 131 4 113 23 3 42 126 23 183 11 3];

        %% pre-lesion days (ampon: Sept18-22)
        /swift6/pk28bk66/ampoff                     % Sept25
        /swift6/pk28bk66/ampoff_b                   % Sept26-29, Oct1-2
        /swift6/pk28bk66/screen                     % Sept14
        /swift6/pk28bk66/trigtest                   % end Sept14, Sept15-17
        /swift6/pk28bk66/newscreen                  % Oct22-31, Nov1-2
        /swift6/pk28bk66/newscreen2                 % Nov4-6,
        % Sept14-17 [2-5]
        % Sept25-28 [12-15]
        % Nov4-6 [16-18]
        % Unlabeled: Sept 29, Nov1-2 , Oct1-2, Oct22-31
        
        %% post-lesion days (ampon: Dec5-10,Dec22-27)
        /swift6/pk28bk66/postlesion/ampoff          % Dec13-21 
        /swift6/pk28bk66/postlesion/ampoff2         % Dec30-31, Jan1
        /swift6/pk28bk66/postlesion/screen          % Nov10-11, Dec1-2, begin Dec3
        /swift6/pk28bk66/postlesion/trigtest        % end Dec3 
        /swift6/pk28bk66/postlesion/trigtest2       % end Dec3, Dec14
        % Nov10-11 [20-21]
        % Dec1-3 [24-26]
        % Dec13-16 [36-39]
        % Dec20 [40]
        % Dec30-31 [49-50]
        % Unlabeled: Dec17-19, Dec21
        
%%%%%%%%%
% pk29bk36 ** DONE
%%%%%%%%%

% keepmatrix = [2:12 25:31 33 35:38 40:49 61:69 75:79];
% firsttimematrix = [89 -17 5 -19 17 -18 186 164 15 -4 -8 -7 22 -4 30 1 -5 24 91 58 16 15 -12 26 123 19 46 26 44 13 16 19 12 -13 61 56 5 -24 -27 15 -24 -18 22 -17 -24 -6 -20];

    %% pre-lesion days (ampon: Feb10-19)
        /swift6/pk29bk36/ampoff                     % Feb22-28
        /swift6/pk29bk36/pretrain                   % Jan26-31, Feb1-3, Feb6-9
    % Jan26-Feb1 [2-8]
    % Feb6-9 [9-12]
    % Feb22-28 [25-31]
    % Unlabeled: Feb2-3, 
    
    %% post-lesion days (ampon:Mar21-29, Apr11-17)
        /swift6/pk29bk36/afterlesion/ampoff         % Apr1-10
        /swift6/pk29bk36/afterlesion/ampoff2        % Apr17-23
        /swift6/pk29bk36/afterlesion/screen         % Mar4-13
        /swift6/pk29bk36/afterlesion/trigtest       % Mar14-20
        /swift6/pk29bk36/afterlesion/trigtest2      % Apr10
    % Mar4-20 [33-49]
    % Apr1-10 [61-69]
    % Apr17-19, 21-22 [75-79]
    % Unlabeled: Apr20, Apr23
    
    
%%%%%%%%%
% pk30bk79  ** DONE
%%%%%%%%%

% keepmatrix = [2:5 18:23 31:37 39 41:44 53:56 63:67];
% firsttimematrix = [-2 -1 5 6 -3 10 8 24 -3 5 -3 13 34 12 2 19 149 19 28 42 12 2 33 19 40 47 106 100 24 26 35];

    %% pre-lesion days (ampon: Dec13-24)
        /swift6/pk30bk79/prelesion/ampoff           % Dec25-Jan4
        /swift6/pk30bk79/prelesion/screen           % Dec9-12
    % Dec9-12 [2:5]
    % Dec25-Dec28 [18:21]
    % Jan3-4 [22:23]
    % Unlabeled: Dec29-31
        
    %% post-lesion days (ampon: Jan29-Feb3)
        /swift6/pk30bk79/postlesion/ampoff          % Feb6-7
        /swift6/pk30bk79/postlesion/ampoff2         % Feb17-22
        /swift6/pk30bk79/postlesion/ampoff3         % Feb24-26
        /swift6/pk30bk79/postlesion/screen          % Jan15-21
        /swift6/pk30bk79/postlesion/trigtest        % Jan23
        /swift6/pk30bk79/postlesion/trigtest2       % Jan25-28
        /swift6/pk30bk79/postlesion/trigtest3       % Feb9-10
    % Jan15-21 [31:37]
    % Jan23 [39]
    % Jan25-28 [41:44]
    % Feb6-7 [53:54]
    % Feb9-10 [55:56]
    % Feb17-21 [63:67]
    % Unlabeled: Feb22, Feb24-26
        
%%%%%%%%%
% pk42bk73  ** DONE
%%%%%%%%%

% keepmatrix = [2:3 9:10 12:14 16:22 24:25 35:38];
% firsttimematrix = [3 -7 10 -12 -15 -10 -20 59 62 0 7 -3 -9 -8 -8 -3 1 -3 86 2];

    %% pre-lesion days (ampon: Jun15-17)
        /swift6/pk42bk73/ampoff                     % Jun20-29 
        /swift6/pk42bk73/screen                     % Jun13
        /swift6/pk42bk73/trigtest1                  % Jun14
        /swift6/pk42bk73/trigtest2                  % end Jun14
    % Jun13-14 [2:3]
    % Jun20, 21, 23 [9 10 12]
    % Jun27-29 [13-15]
    % Unlabeled: Jun 22, Jun24-26, Jun29
    
    %% post-lesion days (ampon: July14-20)
        /swift6/pk42bk73/afterlesion/ampoff         % July 23-26
        /swift6/pk42bk73/afterlesion/screen         % July 2-12
        /swift6/pk42bk73/afterlesion/trigtest2      % July 13 
    % July2, 4-10, 12-13 [16:22 24:25]
    % July 23-26 [35:38]
    % Unlabeled: July3, July 11
        
%%%%%%%%%
% pk80r6  ** DONE
%%%%%%%%%

% keepmatrix = [1:2 10:13 15 21 23:25];
% firsttimematrix = [13 62 16 30 0 19 170 16 18 179 146];

    %% pre-lesion days (ampon: Nov11-16)
        /swift6/pk80r6/prelesion/wnoff1             % Nov19
        /swift6/pk80r6/prelesion/wnoff2-b           % Dec10-13
        /swift6/pk80r6/prelesion/screen2            % Nov8
        /swift6/pk80r6/prelesion/temptest3          % Nov10
    % Nov8, 10 [1:2]
    % Nov19 [10]
    % Dec12-13 [11:12]
    % Unlabeled: Dec10-11
        
    %% post-lesion days (ampon: Dec28-Jan2)
        /swift6/pk30bk79/postlesion/ampoff          % Jan5-13
        /swift6/pk30bk79/postlesion/screen          % Dec17-20, Dec23-24
        /swift6/pk30bk79/postlesion/morescreen      % Jan14
        /swift6/pk30bk79/postlesion/trigtest        % Dec27
        /swift6/pk30bk79/postlesion/trigtest2       % end Dec27
    % Dec24 [13]
    % Dec27 [15]
    % Jan5-7, 9, 13 [21:25]
    % Unlabeled: Dec17-20, Dec23, Jan6, Jan8, Jan10-12, Jan14
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SWIFT7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%%
% o7bk77
%%%%%%%%%

% keepmatrix = [2:5 9 22:27];
% firsttimematrix = [36 31 47 49 45 77 43 37 38 117 41];

    %% pre-lesion days (ampon: ?...can't find folder)
        /swift7/o7bk77/prelesion/screen             % July21-24
    % July21-24 [2:5]
    
    %% post-lesion days (ampon: Aug1-2, Aug5-10)
        /swift7/o7bk77/afterlesion/ampoff           %  Aug14-20
        /swift7/o7bk77/afterlesion/trigtest         %  July31
    % July31 [9]
    % Aug14-19 [22:27]
    % Unlabeled: Aug20    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        