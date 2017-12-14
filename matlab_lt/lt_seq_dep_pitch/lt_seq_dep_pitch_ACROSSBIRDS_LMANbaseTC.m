function lt_seq_dep_pitch_ACROSSBIRDS_LMANbaseTC(DATSTRUCT)

% For each expt, for each class of syl, plot all contexts, across days,
% PBS and MUSC

datstruct = DATSTRUCT.singlesyls;
%%
ListOfBirds = unique(datstruct.Birdnames);
ListOfExpts = unique(datstruct.Exptnames);
ListOfSingleSyls = unique(datstruct.SingleSylAll);

AllAFPbias = [];
AllPitchShiftTommorowPBS = [];
for bname = ListOfBirds
    for ename = ListOfExpts
        indtmp = strcmp(datstruct.Birdnames, bname) & strcmp(datstruct.Exptnames, ename);
        
        if ~any(indtmp)
            % then this bird/expt combo doesn't exist
            continue
        end
        
        
        
        figcount=1;
        subplotrows=4;
        subplotcols=2;
        fignums_alreadyused=[];
        hfigs=[];

        
        
        ListOfSingleSyls = unique(datstruct.SingleSylAll(indtmp));
        
        for ssyl = ListOfSingleSyls
           
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            
% -- find all contexts for this syl
            inddat = find(strcmp(datstruct.Birdnames, bname) & strcmp(datstruct.Exptnames, ename) ...
                & strcmp(datstruct.SingleSylAll, ssyl));
            
            % --- for each context, plot running dat
            for i = 1:length(inddat)
                ind = inddat(i);
                
                % ========= PBS
                plotcol = 'k';
                tvals = datstruct.PitchTvalsPBS{ind};
                fvals = datstruct.PitchFFvalsPBS{ind};
                
                % -- get mean and sem
                [fmean, fsem] = grpstats(fvals, floor(tvals), {'mean', 'sem'});
                tmean = unique(floor(tvals));
                lt_plot(tmean, fmean, {'Errors', fsem, 'Color', plotcol, ...
                    'LineStyle', '-'});
                
                % --- save
                tmeanPBS = tmean;
                fmeanPBS = fmean;
                
                
                % ====== MUSC
                plotcol = 'r';
                tvals = datstruct.PitchTvalsMUSC{ind};
                fvals = datstruct.PitchFFvalsMUSC{ind};
                
                % -- get mean and sem
                [fmean, fsem] = grpstats(fvals, floor(tvals), {'mean', 'sem'});
                tmean = unique(floor(tvals));
                lt_plot(tmean+0.3, fmean, {'Errors', fsem, 'Color', plotcol});
                % -- save
                tmeanMUSC = tmean;
                fmeanMUSC = fmean;
                
                % ====== lines between PBS and MUSC
                for j=1:length(tmeanMUSC)
                   
                    t1 = tmeanMUSC(j);
                    t2 = tmeanMUSC(j)+0.3;
                    
                    f1 = fmeanPBS(tmeanPBS==t1);
                    f2 = fmeanMUSC(j);
                    
                    plot([t1 t2], [f1 f2], '-', 'LineWidth', 2);
                end
                
                % ====== plot text showing FF change
                lt_plot_text(tmeanPBS(end)+0.5, fmeanPBS(end), ...
                    ['shift ' num2str(datstruct.Pitch_MUSC_minus_PBS(ind))], 'r');
                
                
                
                % ############################ COLLECT DATA
                % -- for each musc day, collect 1) AFP bias, 2) pitch shift
                % tomorrow
                for j=1:length(tmeanMUSC)
                   
                    t1 = tmeanMUSC(j);
                    f1 = fmeanPBS(tmeanPBS==t1); % pbs
                    f2 = fmeanMUSC(j); % musc
                    
                    afpbias = f1-f2;
                    
                    ftomorrow = fmeanPBS(tmeanPBS==(t1+1)); % pbs pitch, tomorrow
                    
                    
                    if isempty(ftomorrow)
                        % THEN DON'T KEEP (this is last base day)
                        continue
                    else
                        % == save datapoint
                        AllAFPbias = [AllAFPbias afpbias];
                        AllPitchShiftTommorowPBS = [AllPitchShiftTommorowPBS ...
                            ftomorrow-f1];
  
                    end
                end
                
                
%                 %% ================= paired data [IGNORE FOR NOW]...
%                 for ii=i+1:length(inddat)
%                    ind1 = ind;
%                    ind2 = inddat(ii);
%                     
%                    % ============================== COLLECT DATA FOR EACH
%                    % DAY
%                    for j=1:length(tmeanMUSC)
%                        
%                        % -- for every musc day
%                        datstruct.PitchTvalsPBS
%                        
%                    end
%                    
%                    
%                    
%                 end
                
                
            end
            
            if length(inddat)>1
            title(ssyl, 'Color', 'b');
            else
                title(ssyl);
            end
            
        end
        
        lt_subtitle([bname{1} '-' ename{1}]);
    end
end

%% ==================== PLOT (afp bias vs shift tomorrow)
lt_figure; hold on
xlabel('afp bias');
ylabel('PBS pitch (tomorrow - today)');

plot(AllAFPbias, AllPitchShiftTommorowPBS, 'ok');