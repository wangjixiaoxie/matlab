function [islearning, LearnSummary, switchtime] = lt_neural_v2_QUICK_islearning(birdname, exptname, extractlearninfo)

% extractlearninfo =1, then gets summary fo learning for this expt. if 0, then just tells you if islearning.
% [note: this only runs if actually is laernig]

if ~exist('extractlearninfo', 'var')
    extractlearninfo = 0;
end

%% lt - decides if this was a learning expt

LearnStruct = lt_neural_v2_LoadLearnMetadat;

% is this learning expt?
birdindtmp = strcmp({LearnStruct.bird.birdname}, birdname);
if ~any(birdindtmp)
   islearning = [];
   LearnSummary = [];
   switchtime = [];
    return
    
end
    
islearning = any(strcmp([LearnStruct.bird(birdindtmp).info(1,:)], exptname));

% --- optional: output info about this learn expt
if extractlearninfo==1
        LearnSummary = struct;
        
    if islearning==1
        
        exptinds = strcmp([LearnStruct.bird(birdindtmp).info(1,:)], exptname);
        
        LearnDat = LearnStruct.bird(birdindtmp).info(:, exptinds);
        
        % --- extract vector of datenums (switches) and statuses (pre
        % and post)
        numtargs = size(LearnDat, 2);
        
        for i =1:numtargs
            
            targsyl = LearnDat{2, i};
            LearnSummary.targnum(i).targsyl = targsyl;
            
            switches = LearnDat(3:end, i);
            
            switches = switches(~cellfun('isempty', switches));
            
            % - for each switch extract datenum and statuses
            numswitches = length(switches);
            
            for j=1:numswitches
                
                dnum = datenum(switches{j}(1:14), 'ddmmmyyyy-HHMM');
                statuspre = switches{j}(16:17);
                statuspost = switches{j}(19:20);
                
                LearnSummary.targnum(i).switches(j).datenum = dnum;
                LearnSummary.targnum(i).switches(j).statuspre = statuspre;
                LearnSummary.targnum(i).switches(j).statuspost = statuspost;
                
            end
        end
    end
end

%% ===== output time of WN onset (earliest time across all targs)
if extractlearninfo==1
switchtime = [];
if islearning==1
        numtargs = length(LearnSummary.targnum);
        for jj=1:numtargs
            
            
            if find(strcmp({LearnSummary.targnum(jj).switches.statuspre}, 'Of'), 1, 'first') ~=1
                % then this xperiments started with something other than WN
                % ,,, throw out all data
                indtmp = [];
            else
                indtmp = find(strcmp({LearnSummary.targnum(jj).switches.statuspre}, 'Of'), 1, 'last');
            end
            
            if isempty(indtmp)
                % then WN was always on
                swtimethis = 0; % make this 0 so this will always be earliest, and therefore must throw out all data
            else
                swtimethis = LearnSummary.targnum(jj).switches(indtmp).datenum;
            end
            
            
            % time when WN began
            switchtime = min([switchtime swtimethis]); % get earliest time across all targets
        end
else
end
end