%To know the time range that you analyzed
try
    timerange{1} = datestr(min(check_stuff.(sofinterest).freq.vals(:,1)));
    timerange{2} = datestr(max(check_stuff.(sofinterest).freq.vals(:,1)));
    
    display(timerange)
catch err
end

%Summary information about hit rate
try
    display(' ')
    display(['syllable: ' parameters.sofinterest ' hit rate information'])
    display(check_stuff.(sofinterest).hit)
    display(' ')
catch err
    disp('did not analyze');
end

%Summary information about match rate
try % maybe did not analyze
    display(' ')
    display(['syllable: ' parameters.sofinterest ' match rate information'])
    display(check_stuff.(sofinterest).match)
    display(' ')
catch err
    disp('did not analyze');
end

%Summary information about FF
try % i.e. maybe did not analyze FF
    display(' ')
    display(['syllable: ' parameters.sofinterest ' FF information'])
    display(check_stuff.(sofinterest).freq)
    display(' ')
catch err
    disp('did not analyze');
end

% QUICK SUMMARY
display(' ')
display(['syllable: ' parameters.sofinterest ' QUICK SUMMARY'])
disp('# Labeled')
disp(check_stuff.(sofinterest).match.sum(2));
try
    disp('# WN hit (on target)')
    disp(check_stuff.(sofinterest).hit.sum(1));
    disp('# WN hit (all)')
    disp(check_stuff.(sofinterest).hit.sum(3));
catch err
end
try
    disp('# Offline match (on target, ignoring pitch contingency)')
    disp(check_stuff.(sofinterest).match.sum(1));
    disp('# Offline match (all)')
    disp(check_stuff.(sofinterest).match.sum(3));
catch err
end


% FOR TEMPLATE MATCH ANALYSIS
% Find the song filenames that had labels that were missed
TotMissed=0;
disp('OFFLINE MATCHING:')
for k=1:size(check_stuff.(sofinterest).match.trigger,2);
    NumMissed=check_stuff.(sofinterest).match.values(k,2)-check_stuff.(sofinterest).match.values(k,1);
    if NumMissed>0;
        TotMissed=TotMissed+NumMissed;
        disp([check_stuff.(sofinterest).match.trigger(k).fn ' has ' num2str(NumMissed) ' labels that were not matched offline.']);
        
        % REPORT the inds of the misses:
        inds=find(check_stuff.(sofinterest).match.trigger(k).Labeled_HitsEscapes==0); % finds targets that were missed
        disp(['     Misses (rendition # within song): ' num2str(inds')])
        
    end
    
    % get stats on missed stuff
end

disp(['TOTAL num of missed labels found: ' num2str(TotMissed)]);

% Find false positive songs

syl=[parameters.sofinterest_pre parameters.sofinterest parameters.sofinterest_post];
FalsePos=check_stuff.(sofinterest).match.values(:,3)-check_stuff.(sofinterest).match.values(:,1); % fp equals (hits total) - (hits on targ)
for k=1:length(FalsePos);
    if FalsePos(k)>0;
        disp([check_stuff.(sofinterest).match.trigger(k).fn ' has ' num2str(FalsePos(k)) ' false positive matches.']);
        
        % disp facts about each FP
        tmp=[];
        for i=1:size(check_stuff.(sofinterest).match.trigger(k).tnotes,1); % figure out location of the FP, within song
            if ~strcmp(check_stuff.(sofinterest).match.trigger(k).tnotes(i,:),syl); % if doesn't match, then is FP.
                tmp=[tmp i]; % gives indices of the FPs
            end
        end
        
        
        for kk=1:length(tmp); % for all FPs, get data
            ind=tmp(kk);
            disp(['     FP# ' num2str(kk) ': ' check_stuff.(sofinterest).match.trigger(k).tnotes(ind,:) ' at ' ...
                num2str(check_stuff.(sofinterest).match.trigger(k).ttimes(ind)) 'sec mark'])
        end
    end
end
disp(['TOTAL num of false positives: ' num2str(sum(FalsePos))]);


% FOR ACTUAL WN HIT ANALYSIS.
% (SKIPPED, as this is not informative.  miss likely meens passed pitch
% contingency)
% % Find the song filenames that had labels that were missed
% TotMissed=0;
% disp('WN TARGETING:')
% for k=1:size(check_stuff.(sofinterest).hit.trigger,2);
%     NumMissed=check_stuff.(sofinterest).hit.values(k,2)-check_stuff.(sofinterest).hit.values(k,1);
%     if NumMissed>0;
%         TotMissed=TotMissed+NumMissed;
%         disp([check_stuff.(sofinterest).hit.trigger(k).fn ' has ' num2str(NumMissed) ' labels that were not hit offline.']);
%     end
% end
% disp(['Total num of missed labels found: ' num2str(TotMissed)]);

% Find false positive songs
if get_WN_hits==1;
    disp(' ');
    disp('ONLINE HITS:')
    FalsePos=check_stuff.(sofinterest).hit.values(:,3)-check_stuff.(sofinterest).hit.values(:,2);
    for k=1:length(FalsePos);
        if FalsePos(k)>0;
            disp([check_stuff.(sofinterest).hit.trigger(k).fn ' has ' num2str(FalsePos(k)) ' false positive Hits.']);
            
            % REPORT facts about each FP
            tmp=[];
            for i=1:size(check_stuff.(sofinterest).hit.trigger(k).tnotes,1); % figure out location of the FP, within song
                if ~strcmp(check_stuff.(sofinterest).hit.trigger(k).tnotes(i,:),syl); % if doesn't match, then is FP.
                    tmp=[tmp i]; % gives indices of the FPs
                end
            end
            
            
            for kk=1:length(tmp); % for all FPs, get data
                ind=tmp(kk);
                disp(['FP# ' num2str(kk) ': ' check_stuff.(sofinterest).hit.trigger(k).tnotes(ind,:) ' at ' ...
                    num2str(check_stuff.(sofinterest).hit.trigger(k).ttimes(ind)) 'sec mark'])
            end
            
        end
    end
end

