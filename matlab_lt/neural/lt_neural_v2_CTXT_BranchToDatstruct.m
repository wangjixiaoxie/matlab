function [DATSTRUCT, Diagnostics] = lt_neural_v2_CTXT_BranchToDatstruct(ALLBRANCH, birdstoexclude, ...
    LMANorX, RemoveRepeats, durThreshOmega, dattoplot)
%% lt 10/23/17 - converts branch struct to dat vectors (and stretches using autocorrelation function)

% LMANorX = 0; % 0, both; 1, LMAN; 2, X
% birdstoexclude = {};
% % birdstoexclude = {'bk7', 'bu77wh13', 'or74bk35', 'wh6pk36', 'br92br54'};
%
% % durThreshOmega.syl = 0.15; % omega2 (will only keep if lower) [leave empty to ignore]
% % durThreshOmega.gappre= 0.5;
% % durThreshOmega.gappost= 0.2;
% durThreshOmega.syl = []; % omega2 (will only keep if lower) [leave empty to ignore]
% durThreshOmega.gappre= [];
% durThreshOmega.gappost= [];
%
% RemoveRepeats=0; % if 1, then removes any branch with a class with token preceded by same syl (e.g. a(a)bc or a(a)ab);


%%
numalign = length(ALLBRANCH.alignpos);

DATSTRUCT = struct;

Allalign = [];
Allbirdnum = [];
Allbranchnum = [];
Allneuron = [];

NumRemovedDueToThresh = 0;
NumKeptDueToThresh = 0;
for i=1:numalign
    alignsyl = ALLBRANCH.alignpos(i).alignsyl;
    alignons = ALLBRANCH.alignpos(i).alignonset;
    
    numbirds = length(ALLBRANCH.alignpos(i).bird);
    motifpredur = ALLBRANCH.alignpos(i).ParamsFirstIter.motifpredur;
    
    % === one for each branch/syl
    Xcell = {};
    Ycell = {};
    Xcontcell = {};
    Ycontcell = {};
    
    Xcell_pos = {};
    Ycell_pos = {};
    
    Xcell_neg = {};
    Ycell_neg = {};
    
    BirdnumAll = [];
    NeuronnumAll = [];
    BranchmumAll = [];
    
    Datstruct = struct;
    
    for ii=1:numbirds
        
        
        birdname = ALLBRANCH.SummaryStruct.birds(ii).birdname;
        if any(strcmp(birdstoexclude, birdname))
            disp(['skipping ' birdname]);
            continue
        end
        
        numbranch = length(ALLBRANCH.alignpos(i).bird(ii).branch);
        
        for j=1:numbranch
            numneuron = length(ALLBRANCH.alignpos(i).bird(ii).branch(j).neuron);
            
            for nn=1:numneuron
                
                datneur = ALLBRANCH.alignpos(i).bird(ii).branch(j).neuron(nn);
                
                
                %% if care about location
                if LMANorX==1
                    location = ALLBRANCH.SummaryStruct.birds(ii).neurons(nn).NOTE_Location;
                    % LMAN
                    if ~strcmp(location, 'LMAN')
                        continue
                    end
                elseif LMANorX==2
                    location = ALLBRANCH.SummaryStruct.birds(ii).neurons(nn).NOTE_Location;
                    if ~strcmp(location, 'X')
                        continue
                    end
                end
                
                %% === COLLECT FOR THIS BRANCH/ALIGNMENT
                
                if isempty(datneur.yvals)
                    continue
                end
                if length(datneur.xtimes)==1
                    continue
                end
                
                %% ==== remove anything with repeats?
                if RemoveRepeats ==1
                    isrepeat=0;
                    for k =1:length(datneur.prms_regexpstrlist)
                        thisstr = datneur.prms_regexpstrlist{k};
                        
                        tokensyl = thisstr(strfind(thisstr, '(')+1);
                        presyl = thisstr(strfind(thisstr, '(')-1);
                        
                        if tokensyl == presyl
                            isrepeat =1;
                        end
                        
                        % NOTE: ad hoc, change this to separate section.
                        if presyl =='i'
                            isrepeat=1;
                        end
                    end
                    
                    if isrepeat==1
                        continue
                    end
                end
                
                %% ===== if want to filter by syl/gap duration differences across classes (within context)
                if datneur.DurAnovas.syl_omega>durThreshOmega.syl | ...
                        datneur.DurAnovas.gappre_omega > durThreshOmega.gappre | ...
                        datneur.DurAnovas.gappost_omega > durThreshOmega.gappost
                    NumRemovedDueToThresh = NumRemovedDueToThresh+1;
                    continue
                else
                    NumKeptDueToThresh = NumKeptDueToThresh+1;
                end
                
                
                
                %% ========================================== DATA
                
                DatN = nan(1,3);
                if strcmp(dattoplot, 'classperform')
                    Xcell = [Xcell datneur.xtimes];
                    Ycell = [Ycell datneur.yvals];
                elseif strcmp(dattoplot, 'dprime');
                    % take average for each branch
                    if (1)
                        Ycell = [Ycell nanmean(datneur.DprimeAllPairwise,2)];
                        xtimes = -motifpredur + (1:size(datneur.DprimeAllPairwise,1))./1000;
                        Xcell = [Xcell xtimes'];
                    else
                        % don't take average (each pairwise is one val)
                        Ycell = [Ycell datneur.DprimeAllPairwise];
                        xtimes = -motifpredur + (1:size(datneur.DprimeAllPairwise,1))./1000;
                        xtimes = repmat(xtimes', 1, size(datneur.DprimeAllPairwise,2));
                        Xcell = [Xcell xtimes];
                    end
                elseif strcmp(dattoplot, 'frmean')
                    % take average for each branch
                    xtimes = -motifpredur + (1:size([datneur.FR.classnum.frmean],1))./1000;
                    Xcell = [Xcell xtimes];
                    tmp = nanmean([datneur.FR.classnum.frmean],2);
                    if size(tmp,1)>1
                        tmp = tmp';
                    end
                    Ycell = [Ycell tmp];
                    
                end
                % syl contour
                %                 plot(datneur.sylcontours_x, datneur.sylcontours_mean, '-', 'Color', [0.6 0.2 0.2]);
                Xcontcell = [Xcontcell datneur.sylcontours_x];
%                 datneur.sylcontours_mean = datneur.sylcontours_mean(1:length(datneur.sylcontours_x));
                Ycontcell = [Ycontcell datneur.sylcontours_mean];
                
                DatN(1) = length(Ycell{end});
                
                
                
                % =================================== POS CONTROL
                if strcmp(dattoplot, 'classperform')
                    Xcell_pos = [Xcell_pos datneur.xtimes];
                    Ycell_pos = [Ycell_pos datneur.yvals_pos];
                elseif strcmp(dattoplot, 'dprime');
                    % take average for each branch
                    if (1)
                        Ycell_pos = [Ycell_pos nanmean(datneur.DprimeAllPairwise_Pos,2)];
                        xtimes = -motifpredur + (1:size(datneur.DprimeAllPairwise_Pos,1))./1000;
                        Xcell_pos = [Xcell_pos xtimes'];
                    else
                        % don't take average (each pairwise is one val)
                        Ycell_pos = [Ycell_pos datneur.DprimeAllPairwise_Pos];
                        xtimes = -motifpredur + (1:size(datneur.DprimeAllPairwise_Pos,1))./1000;
                        xtimes = repmat(xtimes', 1, size(datneur.DprimeAllPairwise_Pos,2));
                        Xcell_pos = [Xcell_pos xtimes];
                    end
                elseif strcmp(dattoplot, 'frmean')
                    % take average for each branch
                    xtimes = -motifpredur + (1:size([datneur.FR_POSCONTR.classnum.frmean],1))./1000;
                    Xcell_pos = [Xcell_pos xtimes];
                    
                    tmp = nanmean([datneur.FR_POSCONTR.classnum.frmean],2);
                    if size(tmp,1)>1
                        tmp = tmp';
                    end
                    Ycell_pos = [Ycell_pos tmp];
                    
                end
                DatN(2) = length(Ycell_pos{end});
                
                
                % ========================================= NEG CONTROL
                if strcmp(dattoplot, 'classperform')
                    Xcell_neg = [Xcell_neg datneur.xtimes];
                    Ycell_neg = [Ycell_neg datneur.yvals_neg];
                elseif strcmp(dattoplot, 'dprime');
                    % take average for each branch
                    if (1)
                        Ycell_neg = [Ycell_neg nanmean(datneur.DprimeAllPairwise_Neg,2)];
                        xtimes = -motifpredur + (1:size(datneur.DprimeAllPairwise_Neg,1))./1000;
                        Xcell_neg = [Xcell_neg xtimes'];
                    else
                        % don't take average (each pairwise is one val)
                        Ycell_neg = [Ycell_neg datneur.DprimeAllPairwise_Neg];
                        xtimes = -motifpredur + (1:size(datneur.DprimeAllPairwise_Neg,1))./1000;
                        xtimes = repmat(xtimes', 1, size(datneur.DprimeAllPairwise_Neg,2));
                        Xcell_neg = [Xcell_neg xtimes];
                    end
                elseif strcmp(dattoplot, 'frmean')
                    % take average for each branch [IDENTICAL TO DATA]
                    xtimes = -motifpredur + (1:size([datneur.FR.classnum.frmean],1))./1000;
                    Xcell_neg = [Xcell_neg xtimes];
                    
                                        
                    tmp = nanmean([datneur.FR.classnum.frmean],2);
                    if size(tmp,1)>1
                        tmp = tmp';
                    end
                    Ycell_neg = [Ycell_neg tmp];
                end
                DatN(3) = length(Ycell_neg{end});
                
                
                % =================================== SAMPLE SIZE TALLY
                Allalign = [Allalign i];
                Allbirdnum = [Allbirdnum ii];
                Allbranchnum = [Allbranchnum j];
                Allneuron = [Allneuron nn];
                
                BirdnumAll = [BirdnumAll ii];
                NeuronnumAll = [NeuronnumAll nn];
                BranchmumAll = [BranchmumAll j];

                % ======= confirm that is paired - i.e. each datapoint has dat and both controls
                %                 disp(DatN);
                assert(length(unique(DatN))==1, 'dat and controls have diff lengths ...');
            end
        end
    end
    
    assert(length(Ycell) == length(Ycell_pos), 'asfds');
    assert(length(Ycell) == length(Ycell_neg), 'asdf');
    
    % ================== PUT INTO STRUCT
    Datstruct.Dat.Xcell = Xcell;
    Datstruct.Dat.Ycell = Ycell;
    Datstruct.Dat.Xcontcell = Xcontcell;
    Datstruct.Dat.Ycontcell = Ycontcell;
    
    Datstruct.PosContr.Xcell = Xcell_pos;
    Datstruct.PosContr.Ycell = Ycell_pos;
    
    Datstruct.NegContr.Xcell = Xcell_neg;
    Datstruct.NegContr.Ycell = Ycell_neg;
    
    Datstruct.Dat.BirdnumAll = BirdnumAll;
    Datstruct.Dat.NeuronnumAll = NeuronnumAll;
    Datstruct.Dat.BranchmumAll = BranchmumAll;
    
    %% ====================== SUBTRACT CONTROLS
    
    for nn = 1:length(Datstruct.Dat.Ycell);
        % ============================ VS. NEG
        % --- confirm that xvalues match
        assert(all(Datstruct.Dat.Xcell{nn} - Datstruct.NegContr.Xcell{nn} < 0.001), 'safasd');
        Datstruct.DatMinusNeg.Xcell{nn} = Datstruct.Dat.Xcell{nn};
        
        % -- get diff
        ydiff = Datstruct.Dat.Ycell{nn} - Datstruct.NegContr.Ycell{nn};
        Datstruct.DatMinusNeg.Ycell{nn} = ydiff;
        
        
        % ============================== VS POS
        % --- confirm that xvalues match
        assert(all(Datstruct.Dat.Xcell{nn} - Datstruct.PosContr.Xcell{nn} < 0.001), 'safasd');
        Datstruct.DatMinusPos.Xcell{nn} = Datstruct.Dat.Xcell{nn};
        
        % -- get diff
        ydiff = Datstruct.Dat.Ycell{nn} - Datstruct.PosContr.Ycell{nn};
        Datstruct.DatMinusPos.Ycell{nn} = ydiff;
    end
    
    
    
    %% ========= time warping, align all neur by warping syl contours
    % ============= METHOD 1, linear stretching, based on autocorrelation
    % width between peaks
    
    % --------- Actual dat
    dattype = 'Dat';
    Datstruct = warp_stretch(Datstruct, dattype);
    
    dattype = 'NegContr';
    Datstruct = warp_stretch(Datstruct, dattype);
    
    dattype = 'PosContr';
    Datstruct = warp_stretch(Datstruct, dattype);
    
    dattype = 'DatMinusNeg';
    Datstruct = warp_stretch(Datstruct, dattype);
    
    dattype = 'DatMinusPos';
    Datstruct = warp_stretch(Datstruct, dattype);
    
    
    %%
    
    DATSTRUCT.numalign(i).Datstruct = Datstruct;
    
    
end

Diagnostics.Allalign =Allalign;
Diagnostics.Allbirdnum = Allbirdnum;
Diagnostics.Allbranchnum = Allbranchnum;
Diagnostics.Allneuron = Allneuron;

Diagnostics.NumRemovedDueToThresh = NumRemovedDueToThresh;
Diagnostics.NumKeptDueToThresh = NumKeptDueToThresh;


end



function Datstruct = warp_stretch(Datstruct, dattype)
%  WILL STRETCH RELATIVE TO ACTUAL DAT FOR THAT BRANCH

% -- MODIFY
Ycell = Datstruct.(dattype).Ycell; % will be stretched
Xcell = Datstruct.(dattype).Xcell;

% -------------- HARD PARAMS
stretchtempl = 'Dat';
Ycontcell = Datstruct.(stretchtempl).Ycontcell; % will determing how much to stretch
Xcontcell = Datstruct.(stretchtempl).Xcontcell;

pkwidthtarg = 100; % 120 ms

% ----- RUN
numsamps = length(Ycell);
Ycell_WARP = cell(size(Ycell));
Ycontcell_WARP = cell(size(Ycontcell));

for j=1:numsamps
    
    tmp = xcorr(Ycontcell{j}(~isnan(Ycontcell{j})));
    tmp = tmp(floor(end/2):end);
    [~, pklocs] = findpeaks(double(tmp), 'sortstr', 'descend', 'npeaks', 2, 'minpeakdistance', 20);
    pkwidth = pklocs(2)-pklocs(1);
    
    % ====== stretch based on pkwidth compared to target
    % --- performance
    xtimes = Xcell{j};
    yvals = Ycell{j};
    if length(xtimes)<3
        Ycell_WARP{j} = Ycell{j};
        Ycontcell_WARP{j} = Ycontcell{j};
        continue
    end
    xtimes_new = xtimes.*pkwidthtarg/pkwidth;
    % resample at the original xtimes
    yvals_new = interp1(xtimes_new, yvals, xtimes);
    
    % -- contour
    xtimes_cont = Xcontcell{j};
    ycont = Ycontcell{j};
    xtimes_cont_new = xtimes_cont.*pkwidthtarg/pkwidth;
    ycont_new = interp1(xtimes_cont_new, ycont, xtimes_cont);
    
    if (0)
        if rand<0.1
            lt_figure; hold on;
            
            subplot(311); hold on;
            title('class performance (bk = original)');
            plot(xtimes, yvals, '-ok');
            plot(xtimes_new, yvals, '-r');
            plot(xtimes, yvals_new, '-or');
            
            subplot(312); hold on;
            title('syl contours');
            plot(xtimes_cont, ycont, '-xk');
            plot(xtimes_cont_new, ycont, '-r');
            plot(xtimes_cont, ycont_new, '-xr');
            
            subplot(313);  hold on;
            title('autocor')
            plot(tmp);
            line([pkwidthtarg pkwidthtarg], ylim, 'Color','k');
            lt_plot_annotation(1, ['wdth:' num2str(pkwidth)], 'r')
            pause; close all;
        end
    end
    
    % ==== OUTPUT (overwrite)
    Ycell_WARP{j} = yvals_new;
    Ycontcell_WARP{j} = ycont_new;
end

% ---- dat
Datstruct.(dattype).Ycell_WARP = Ycell_WARP;
% -- template
Datstruct.(stretchtempl).Ycontcell_WARP = Ycontcell_WARP;

end
