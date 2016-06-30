%% LT 11/25/15 - a) keeps only repeats, b) throws out repeats
function SeqDepPitch_AcrossBirds=lt_seq_dep_pitch_ACROSSBIRDS_FilterRepeats(SeqDepPitch_AcrossBirds, remove_repeats, extract_repeats)
% SHOULD DO: [not yet]
%

% NOTE ON WHAT DOES:
% remove_repeats=2; [syl pos 2 +] a) throws out all experiments if target was syl2 or higher in repeat, and b) For all experiments, throw out syls that are at repeat position 2+ [EVEN IF IS DIFF TYPE SYL] [CAN END UP WITH EXPTS WITH NO/FEW NONTARGS]
% remove_repeats=3; [same as above, but for repeats with syls 3+];
% remove_repeats=9; if targ is in repeat, remove syls in that repeat (save this for direct analysis of experiments targeting repeats), but keep all other syls
% for all repeats across experiments, remove syls in repeats past certain number (pos 3)
% remove_repeats=10; if targ is in repeat, remove syls in that repeat coming after targ, but keep all other syls
% remove_repeats=0; does not alter SeqDepPitch_AcrossBirds

% extract_repeats=1; CAN ONLY DO ONE of extract or remove in one go (code is fool proof) ---> Extract list of repeat experiments - i.e. if target is in repeat, what is 1) position of targ in repeat, list of syls in that repeat

% %% == save original
%
% SeqDepPitch_AcrossBirds_ORIG=SeqDepPitch_AcrossBirds;
%
%% ++++++++++++++++++++++ REMOVING REPEATS [2+]
if remove_repeats==2;
    
    %% ==== 1) Experiment automatically thrown out if the target was syl 2 or higher in repeat
    
    filter = 'ThrowOutIfTargInPosTwoOrHigherInRepeat';
    [SeqDepPitch_AcrossBirds, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds, filter);
    
    
    
    
    %% ===== 2) For all experiments, throw out syls that are at repeat position 2+ [EVEN IF IS DIFF TYPE SYL]
    for i=1:length(SeqDepPitch_AcrossBirds.birds)
        
        birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
        numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        
        for ii=1:numexpts
            exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
            SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
            
            % GO THROUGH ALL SYLS AND REMOVE THEM IF THEY ARE POS 2+ IN REPEATS
            %         disp([i ii]);
            SylsUnique_temp={};
            disp(' ');
            disp([birdname '-' exptname ' removed:']);
            for j=1:length(SylsUnique);
                syl=SylsUnique{j};
                
                
                presyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).preceding_syl;
                presyl=lower(presyl);
                
                % IF PRESYL IS SAME AS SYL, THEN THIS IS POS 2+ IN A REPEAT
                syl_lower=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).single_syl;
                
                
                if strcmp(syl_lower, presyl);
                    % then remove this syl
                    disp(syl);
                else
                    % keep this syl
                    SylsUnique_temp=[SylsUnique_temp syl];
                    
                    
                end
            end
            
            disp(['kept: ' SylsUnique_temp]);
            
            % ===== OUTPUT
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique=SylsUnique_temp;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.SylFields_Unique=SylsUnique_temp;
            
        end
    end
    
    
    
end



%% ++++++++++++++++++++++ REMOVING REPEATS [3+]
if remove_repeats==3;
    
    %% ==== 1) Experiment automatically thrown out if the target was syl 3 or higher in repeat
    
    filter = 'ThrowOutIfTargInPosThreeOrHigherInRepeat';
    [SeqDepPitch_AcrossBirds, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds, filter);
    
    
    
    
    %% ===== 2) For all experiments, throw out syls that are at repeat position 3+ [EVEN IF IS DIFF TYPE SYL]
    for i=1:length(SeqDepPitch_AcrossBirds.birds)
        
        birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
        numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        
        for ii=1:numexpts
            exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
            SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
            
            % GO THROUGH ALL SYLS AND REMOVE THEM IF THEY ARE POS 2+ IN REPEATS
            %         disp([i ii]);
            SylsUnique_temp={};
            disp(' ');
            disp([birdname '-' exptname ' removed:']);
            for j=1:length(SylsUnique);
                syl=SylsUnique{j};
                
                twosylback=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).two_syl_back;
                twosylback=lower(twosylback);
                
                presyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).preceding_syl;
                presyl=lower(presyl);
                
                syl_lower=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).single_syl;
                
                % if this syl, presyl, and 2sylback are same, then this is
                % pos 3 in repeat
                if strcmp(twosylback, presyl) & strcmp(presyl, syl_lower)
                    disp(syl);
                else
                    % keep this syl
                    SylsUnique_temp=[SylsUnique_temp syl];
                end
                
            end
            
            disp(['kept: ' SylsUnique_temp]);
            
            % ===== OUTPUT
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique=SylsUnique_temp;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.SylFields_Unique=SylsUnique_temp;
            
        end
    end
end

if remove_repeats==9;
    %% 1) if the target is in a repeat, then remove all same-type syls
    % that directly precede or follow it
    disp('------------ REMOVING SYLS THAT ARE IN REPEATS, IN SAME MOTIF AS TARGET, IF TARGET WAS IN REPEAT')
    for i=1:length(SeqDepPitch_AcrossBirds.birds);
        birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
        NumExperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        
        for ii=1:NumExperiments;
            exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
            
            targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
            targsyl_lower=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).single_syl;
            
            MotifNum=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).motif_num;
            MotifSyls=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.SylLists.FieldsInOrder{MotifNum};
            
            % === position of targ in motif
            TargPosInMotif=find(strcmp(MotifSyls, targsyl));
            
            % === convert all in motif to single syls
            MotifSyls_lower={};
            for j=1:length(MotifSyls);
                syltmp=MotifSyls{j};
                
                if length(syltmp)>1;
                    syltmp=regexp(syltmp, '[A-Z]', 'match');
                end
                syltmp=lower(syltmp);
                
                MotifSyls_lower=[MotifSyls_lower syltmp];
            end
            
            % === go backwards, if find similar syl, remove. stop once find
            % diff syl
            SylIndsToRemove=[];
            for k=1:25
                if k>=TargPosInMotif
                    continue
                end
                
                SylToCheck_lower=MotifSyls_lower{TargPosInMotif-k};
                
                if strcmp(SylToCheck_lower, targsyl_lower);
                    % then remove this syl and keep going
                    SylIndsToRemove=[SylIndsToRemove TargPosInMotif-k];
                    
                else
                    % then reached edge of repeat, stop
                    break
                end
                
            end
            
            % === go forward, if find similar syl, remove. stop once find
            % diff syl
            for k=1:25
                if k+TargPosInMotif>length(MotifSyls_lower)
                    continue
                end
                
                SylToCheck_lower=MotifSyls_lower{TargPosInMotif+k};
                
                if strcmp(SylToCheck_lower, targsyl_lower);
                    % then remove this syl and keep going
                    SylIndsToRemove=[SylIndsToRemove TargPosInMotif+k];
                    
                else
                    % then reached edge of repeat, stop
                    break
                end
                
            end
            
            
            % ===== WHAT SYLS TO REMOVE?
            SylsToRemove_syls=MotifSyls(SylIndsToRemove);
            
            
            
            % ===== ACTUALLY REMOVE
            SylsUnique_tmp={};
            SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
            
            for k=1:length(SylsUnique);
                syl=SylsUnique{k};
                
                if any(strcmp(SylsToRemove_syls, syl));
                    % then don't keep this syl
                else
                    SylsUnique_tmp=[SylsUnique_tmp syl];
                end
            end
            
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique=SylsUnique_tmp;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.SylFields_Unique=SylsUnique_tmp;
            
            
            % ===== display
            disp([birdname '-' exptname ':']);
            disp(['target: ' targsyl , ' in motif: ' MotifSyls]);
            disp(['removed: ' SylsToRemove_syls]);
            disp(['kept: ' SylsUnique_tmp]);
        end
    end
    
    
    
    %% ===== 2) For all experiments, throw out syls that are at repeat position 3+ [EVEN IF IS DIFF TYPE SYL]
    disp(' ---------------------- REMOVING ALL SYLS AT POS 3+')
    for i=1:length(SeqDepPitch_AcrossBirds.birds)
        
        birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
        numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        
        for ii=1:numexpts
            exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
            SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
            
            % GO THROUGH ALL SYLS AND REMOVE THEM IF THEY ARE POS 2+ IN REPEATS
            %         disp([i ii]);
            SylsUnique_temp={};
            disp(' ');
            disp([birdname '-' exptname ' removed:']);
            for j=1:length(SylsUnique);
                syl=SylsUnique{j};
                
                twosylback=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).two_syl_back;
                twosylback=lower(twosylback);
                
                presyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).preceding_syl;
                presyl=lower(presyl);
                
                syl_lower=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).single_syl;
                
                % if this syl, presyl, and 2sylback are same, then this is
                % pos 3 in repeat
                if strcmp(twosylback, presyl) & strcmp(presyl, syl_lower)
                    disp(syl);
                else
                    % keep this syl
                    SylsUnique_temp=[SylsUnique_temp syl];
                end
                
            end
            
            disp(['kept: ' SylsUnique_temp]);
            
            % ===== OUTPUT
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique=SylsUnique_temp;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.SylFields_Unique=SylsUnique_temp;
            
        end
    end
end


if remove_repeats==10;
    %% 1) if the target is in a repeat, then remove all same-type syls that come after
    % that directly precede or follow it
    disp('------------ REMOVING SYLS THAT ARE IN REPEATS, IN SAME MOTIF AS TARGET, IF TARGET WAS IN REPEAT')
    for i=1:length(SeqDepPitch_AcrossBirds.birds);
        birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
        NumExperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        
        for ii=1:NumExperiments;
            exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
            
            targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
            targsyl_lower=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).single_syl;
            
            MotifNum=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).motif_num;
            MotifSyls=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.SylLists.FieldsInOrder{MotifNum};
            
            % === position of targ in motif
            TargPosInMotif=find(strcmp(MotifSyls, targsyl));
            
            % === convert all in motif to single syls
            MotifSyls_lower={};
            for j=1:length(MotifSyls);
                syltmp=MotifSyls{j};
                
                if length(syltmp)>1;
                    syltmp=regexp(syltmp, '[A-Z]', 'match');
                end
                syltmp=lower(syltmp);
                
                MotifSyls_lower=[MotifSyls_lower syltmp];
            end
            
            SylIndsToRemove=[];
            %             % === go backwards, if find similar syl, remove. stop once find
            %             % diff syl
            %             for k=1:25
            %                 if k>=TargPosInMotif
            %                     continue
            %                 end
            %
            %                 SylToCheck_lower=MotifSyls_lower{TargPosInMotif-k};
            %
            %                 if strcmp(SylToCheck_lower, targsyl_lower);
            %                     % then remove this syl and keep going
            %                     SylIndsToRemove=[SylIndsToRemove TargPosInMotif-k];
            %
            %                 else
            %                     % then reached edge of repeat, stop
            %                     break
            %                 end
            %
            %             end
            
            % === go forward, if find similar syl, remove. stop once find
            % diff syl
            for k=1:25
                
                if k+TargPosInMotif>length(MotifSyls_lower)
                    continue
                end
                
                SylToCheck_lower=MotifSyls_lower{TargPosInMotif+k};
                
                if strcmp(SylToCheck_lower, targsyl_lower);
                    % then remove this syl and keep going
                    SylIndsToRemove=[SylIndsToRemove TargPosInMotif+k];
                    
                else
                    % then reached edge of repeat, stop
                    break
                end
                
            end
            
            
            % ===== WHAT SYLS TO REMOVE?
            SylsToRemove_syls=MotifSyls(SylIndsToRemove);
            
            
            % ===== ACTUALLY REMOVE
            SylsUnique_tmp={};
            SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
            
            for k=1:length(SylsUnique);
                syl=SylsUnique{k};
                
                if any(strcmp(SylsToRemove_syls, syl));
                    % then don't keep this syl
                else
                    SylsUnique_tmp=[SylsUnique_tmp syl];
                end
            end
            
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique=SylsUnique_tmp;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.SylFields_Unique=SylsUnique_tmp;
            
            
            % ===== display
            disp([birdname '-' exptname ':']);
            disp(['target: ' targsyl , ' in motif: ' MotifSyls]);
            disp(['removed: ' SylsToRemove_syls]);
            disp(['kept: ' SylsUnique_tmp]);
        end
    end
    
end

%% +++++++++++++++++++++++++++ EXTRACTING REPEATS
% 1b) Extract list of repeat experiments - i.e. if target is in repeat, what is 1) position of targ in repeat, list of syls in that repeat
if extract_repeats==1 && remove_repeats~=0;
    disp('----- SKIPPING EXTRACT, since remove_repeats~=0');
elseif extract_repeats==1 && remove_repeats==0;
    
    disp('------------ EXTRACTING INFOR OF EXPERIMENTS WITH TARGS IN REPEATS')
    for i=1:length(SeqDepPitch_AcrossBirds.birds);
        birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
        NumExperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        
        for ii=1:NumExperiments;
            exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
            
            targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
            targsyl_lower=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).single_syl;
            
            MotifNum=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).motif_num;
            MotifSyls=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.SylLists.FieldsInOrder{MotifNum};
            
%             % ==== REMOVE ANY MOTIF SYLS THAT ARE NOT DEFINED IN UNIQUE
%             % SYLS
%             SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
%             MotifSyls_tmp={};
%             for j=1:length(MotifSyls);
%                 motifsyl=MotifSyls{j};
%                 
%                 if any(strcmp(SylsUnique, motifsyl));
%                     % then keep
%                     MotifSyls_tmp=[MotifSyls_tmp motifsyl];
%                 end
%             end
%             MotifSyls=MotifSyls_tmp;
%                 
            
            % === position of targ in motif
            TargPosInMotif=find(strcmp(MotifSyls, targsyl));
            
            % === convert all in motif to single syls
            MotifSyls_lower={};
            for j=1:length(MotifSyls);
                syltmp=MotifSyls{j};
                
                if length(syltmp)>1;
                    syltmp=regexp(syltmp, '[A-Z]', 'match');
                end
                syltmp=lower(syltmp);
                
                MotifSyls_lower=[MotifSyls_lower syltmp];
            end
            
            % === go backwards, if find similar syl, remove. stop once find
            % diff syl
            SylIndsToRemove=[];
            for k=1:25
                if k>=TargPosInMotif
                    continue
                end
                
                SylToCheck_lower=MotifSyls_lower{TargPosInMotif-k};
                
                if strcmp(SylToCheck_lower, targsyl_lower);
                    % then remove this syl and keep going
                    SylIndsToRemove=[SylIndsToRemove TargPosInMotif-k];
                    
                else
                    % then reached edge of repeat, stop
                    break
                end
                
            end
            
            % === go forward, if find similar syl, remove. stop once find
            % diff syl
            for k=1:25
                
                if k+TargPosInMotif>length(MotifSyls_lower)
                    continue
                end
                
                SylToCheck_lower=MotifSyls_lower{TargPosInMotif+k};
                
                if strcmp(SylToCheck_lower, targsyl_lower);
                    % then remove this syl and keep going
                    SylIndsToRemove=[SylIndsToRemove TargPosInMotif+k];
                    
                else
                    % then reached edge of repeat, stop
                    break
                end
                
            end
            
            
            % ==== SKIP THIS EXPT IF NOT IN REPEAT
            if isempty(SylIndsToRemove)
                disp(' ');
                disp([birdname '-' exptname]);
                disp(['target: ' targsyl]);
                disp('NO REPEAT')
                continue;
            end
            
            % ==== CONVERT THOSE INDS INTO POSITION IN REPEAT
            IndsOfRepeat=sort([SylIndsToRemove TargPosInMotif]);
            
            PosOfTargInRepeat=find(IndsOfRepeat==TargPosInMotif);
            SylsInRepeat=MotifSyls(IndsOfRepeat);
            
            
            
            disp(' ');
            disp([birdname '-' exptname]);
            disp(['target: ' targsyl]);
            disp(['syls in repeat: ' SylsInRepeat]);
            disp(['position of target in repeat: ' num2str(PosOfTargInRepeat)]);
            
            
            % === save
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DataRepeats.PosOfTargInRepeat=PosOfTargInRepeat;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DataRepeats.SylsInRepeat=SylsInRepeat;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DataRepeats.targsyl=targsyl;
            
            
            
        end
    end
    
end

