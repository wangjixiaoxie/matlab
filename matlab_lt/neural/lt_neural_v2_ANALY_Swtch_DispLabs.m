function lt_neural_v2_ANALY_Swtch_DispLabs(MOTIFSTATS_Compiled, SwitchStruct, onlyWNonset)
%%

%%
% ---
numbirds = length(MOTIFSTATS_Compiled.birds);

for bb = 1:numbirds
    
    numexpts = length(MOTIFSTATS_Compiled.birds(bb).exptnum);
    birdname = MOTIFSTATS_Compiled.birds(bb).birdname;
    for ee =1:numexpts
        
        numswitches = length(SwitchStruct.bird(bb).exptnum(ee).switchlist);
        motifstats = MOTIFSTATS_Compiled.birds(bb).exptnum(ee).MOTIFSTATS;
        exptname = SwitchStruct.bird(bb).exptnum(ee).exptname;
        
        for ss = 1:numswitches
            
            % ===================== only continue if is onset of experiment
           if onlyWNonset==1
            tmptmp = SwitchStruct.bird(bb).exptnum(ee).switchlist(ss).learningContingencies(2:2:end);
            tmptmp = cell2mat(tmptmp');
            if ~all(tmptmp(:,1)==0)
                disp('SKIPPED - not starting from WN off')
                continue
            end
           end
           
            % ======= note down experiment type
            singlesyls = motifstats.params.SingleSyls_inorder;
            learnconting = SwitchStruct.bird(bb).exptnum(ee).switchlist(ss).learningContingencies;
            
            expttype = [];
            numtargs = [];
            if length(learnconting) ==2
                numtargs = 1;
                % then only one targ - other contexts exist?
                targsinglesyl = learnconting{1}(findstr(learnconting{1}, '(')+1);
                if sum(strcmp(singlesyls, targsinglesyl))>1
                    expttype =   'one targ context';
                elseif sum(strcmp(singlesyls, targsinglesyl)) ==1
                    expttype = 'one targ - stereotyped';
                else
                    asdfhasdlhihildsaf;
                end
            elseif length(learnconting)>2
                numtargs = length(learnconting)/2;
                % then multiple targets - are they same direction?
                if length(unique(cell2mat(SwitchStruct.bird(bb).exptnum(ee).switchlist(ss).learningDirs(2:2:end)))) ==1
                    expttype = 'mult targ context - samedir';
                else
                    expttype = 'mult targ context - diff dir';
                end
            end
            
            % ########################## display things about this expt
            disp(' ');
            disp(' ############################################## ');
            disp([birdname '-' exptname '-sw' num2str(ss)]);
            tmp = datestr(SwitchStruct.bird(bb).exptnum(ee).switchlist(ss).switchdnum, 'ddmmmyyyy-HHMM');
            disp(['time of switch: ' tmp]);
            disp(['expt type: ' expttype]);
            disp(['targets: ' SwitchStruct.bird(bb).exptnum(ee).switchlist(ss).learningDirs]);
            disp(['same type syls: ' motifstats.params.SameTypeSyls]);
            
            % #####################
            
            % ==========================================================
            % collect data for each neuron
            goodneurons = SwitchStruct.bird(bb).exptnum(ee).switchlist(ss).goodneurons;
            
            for nn=goodneurons
                
                motiflist = motifstats.params.motif_regexpr_str;
                nummotifs = length(SwitchStruct.bird(bb).exptnum(ee).switchlist(ss).neuron(nn).DATA.motif);
                mcount = 0;
                %                 nummotifs = length(motifstats.params.motif_regexpr_str);
                
                % ################################# display song labeling
                % progress for this neuron
                
                cd(MOTIFSTATS_Compiled.birds(bb).exptnum(ee).SummaryStruct.birds(1).neurons(nn).dirname);
                batchfile = MOTIFSTATS_Compiled.birds(bb).exptnum(ee).SummaryStruct.birds(1).neurons(nn).batchfilename;
                eval(['!cp ' batchfile ' ..']); % copy batch file one dir up. [as this is most accurate batch file]
                cd .. % labeled files are up one.
                disp(['neuron ' num2str(nn) ' out of ' num2str(max(goodneurons))]);
                disp(MOTIFSTATS_Compiled.birds(bb).exptnum(ee).SummaryStruct.birds(1).neurons(nn).dirname)
                
                lt_batch_disp_all_labels(batchfile);
                
            end
        end
    end
    
    disp('enter press anything to continue');
    input(' ');
    
    
end