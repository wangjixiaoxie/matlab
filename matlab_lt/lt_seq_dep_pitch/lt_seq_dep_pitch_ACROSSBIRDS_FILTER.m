function [SeqDepPitch_AcrossBirds, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_FILTER(SeqDepPitch_AcrossBirds, PARAMS)
%% LT  - filters out information for each syl in each experiment into SylID field.

%% PARAMS
NumBirds=length(SeqDepPitch_AcrossBirds.birds);

%% HIT RATE OF ALL EXPERIMENTS

for i=1:NumBirds;
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        
        syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        for j=1:length(syls_unique);
            syl=syls_unique{j};
            num_days=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).HitRateStatus);
            
            NumHits=[];
            NumTotal=[];
            HitRate=[];
            for jj=1:num_days;
                NumHits(jj)=sum(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).HitRateStatus{jj});
                NumTotal(jj)=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).HitRateStatus{jj});
            end
            HitRate=NumHits./NumTotal;
            
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed.(syl).HIT_RATE_overdays.NumHits=NumHits;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed.(syl).HIT_RATE_overdays.NumTotal=NumTotal;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PostProcessed.(syl).HIT_RATE_overdays.HitRate=HitRate;
        end
    end
end





%% FILTER ALL SYLLABLES (i.e. giving each one vector of features)
% THIS CURRENTLY ONLY WORKS FOR SINGLE TARGET EXPERIMENTS
disp('NOTE - FIlter will operate on all expts (even those with 2 targets), but is only reliable for One target expts currently (take the first targsyl as the target)!!');
disp('NOTE - Only gets data for syls that are defined in "fields in order" (e.g. a single syl will not be analyzed)');
% format as a structure.

% first get all syllable names without redundancy - use motifs that I have
% hand designated

counter1=0;
counter2=0;

for i=1:NumBirds;
    NumExperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    for ii=1:NumExperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        % ============================ what motifs are in this expt?
        motifs=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.SylLists.FieldsInOrder;
        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.motifs=motifs;
        
        
        % ======== GET INFO ABOUT TARGET FOR THIS EXPT -------------------------
        targsyl=[];
        targsyl_pre=[];
        targsyl_post=[];
        
        % targ syl
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.SylLists.TargetSyls{1};
        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl=targsyl;
        
        % - preceding syl
        try % try: sometimes targyl is just one letter. sometimes is >1 letter, but is 1st in the string. that gives error
            if length(targsyl)>1;
                tmp=regexp(targsyl,'[A-Z]'); % find caps
                if tmp==1;
                    targsyl_pre=nan;
                else
                    targsyl_pre=targsyl(tmp-1);
                end
            end
        catch err
            disp('error: targsyl is not multiletter string, or is at start of string'); % if this is triggered, use method used for post syl, below.
            keyboard
        end
        
        % - post syl
        for k=1:length(motifs);
            if any(strcmp(motifs{k},targsyl));
                % find index of target
                ind=find(strcmp(motifs{k},targsyl));
                
                % post syl is one ahead
                if ind<length(motifs{k})
                    targsyl_post=motifs{k}{ind+1};
                    if length(targsyl_post)>1; % just get the single syl
                        ind=regexp(targsyl_post,'[A-Z]');
                        targsyl_post=lower(targsyl_post(ind));
                    else
                        targsyl_post=lower(targsyl_post);
                    end
                else
                    % then this is the last syl in motif, no post syl
                    targsyl_post=nan;
                end
            end
        end
        
        % - convert targ syl to single syllable
        if length(targsyl)>1
            % where is the upper case?
            upperind=regexp(targsyl,'[A-Z]');
            single_syl_targ=lower(targsyl(upperind));
        else
            % then could be either upper or lower
            single_syl_targ=lower(targsyl);
        end
        % - check that we got everything.
        if isempty(targsyl) | isempty(targsyl_pre) | isempty(targsyl_post);
            disp('problem, see code');
            keyboard
        end
        
        % ----- PCA SCORE OF TARGET
        fv_PCAscore_target_mean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(targsyl).fv_baseline_PCAscore_mean;
        fv_zscore_target_mean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(targsyl).fv_baseline_zscore_mean;
        % ======================================================
        
        % ============================ EXTRACT INFO for each syl.
        SylFields_Unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        for j=1:length(SylFields_Unique);
            
            % -- extract individual syllables
            syl=SylFields_Unique{j};
            
            % what is the syl (ignoring context)?
            if length(syl)>1
                % where is the upper case?
                upperind=regexp(syl,'[A-Z]');
                single_syl=lower(syl(upperind));
            else
                % then could be either upper or lower
                single_syl=lower(syl);
            end
            
            
            % is this the target?
            if strcmp(targsyl,syl)==1;
                is_target=1;
            else
                is_target=0;
            end
            
            
            % is that syl similar or different to target?
            if strcmp(single_syl,single_syl_targ)==1;
                similar_to_targ=1;
            else
                similar_to_targ=0;
            end
            
            
            % === Acoustic features
            % --- Euclidian distance from target (PCA space)
            fv_PCA_mean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl).fv_baseline_PCAscore_mean;
            eucldist_from_targ_PCA=sqrt(sum((fv_PCA_mean-fv_PCAscore_target_mean).^2));
            
            fv_zscore_mean=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl).fv_baseline_zscore_mean;
            eucldist_from_targ_zscore=sqrt(sum((fv_zscore_mean-fv_zscore_target_mean).^2));
            
            
            % =======================MOTIF RELATED STATS
            % ID the motif
            motif_num=nan;
            for iii=1:length(motifs);
                if any(strcmp(motifs{iii}, syl));
                    % then this is the motif
                    motif_num=iii;
                end
            end
            
            if isnan(motif_num);
                % then this is a syl without defined motif
                % set all outputs to nan
                is_in_defined_motif=0;
                pre_syl=nan;
                presyl_similar_to_targ_presyl=0;
                post_syl=nan;
                postsyl_similar_to_targ_postsyl=nan;
                distance_from_targ=nan;
                two_syl_back=nan;
                
            else
                % then this syllable has a defined motif
                is_in_defined_motif=1;
                
                syls_in_motif=motifs{motif_num};
                
                % what is this syl's position in motif;
                syl_pos_in_motif=find(strcmp(syls_in_motif, syl));
                
                % what is preceding syl?
                % --- first try to get based on syl name
                %                 if length(syl)>1;
                %                     indCaps=regexp(syl, '[A-Z]');
                %                     preceding_syl=lower(syl(indCaps-1));
                %                 elseif syl_pos_in_motif==1;
                %                     % --- then try to extract using motif
                %                     % then this is the first syl
                %                     preceding_syl=nan;
                %                 else
                %                     preceding_syl=syls_in_motif{syl_pos_in_motif-1};
                %                 end
                %
                
                
                pre_syl=nan;
                % first choice
                indCaps=regexp(syl, '[A-Z]');
                if indCaps>1;
                    pre_syl=lower(syl(indCaps-1));
                else
                    if syl_pos_in_motif==1;
                        % second choice
                        
                        % --- then try to extract using motif
                        % then this is the first syl
                        pre_syl=nan;
                    else
                        pre_syl=syls_in_motif{syl_pos_in_motif-1};
                    end
                end
                
                %  what is 2 syls preceding?
                % --- first try using syl name
                %                 if length(syl)>2
                %                     indCaps=regexp(syl, '[A-Z]');
                %                     two_syl_back=lower(syl(indCaps-2));
                %                 elseif syl_pos_in_motif>2;
                %                     % --- else try to use motif
                %                     two_syl_back=syls_in_motif{syl_pos_in_motif-2};
                %                 else
                %                     two_syl_back=nan;
                %                 end
                
                two_syl_back=nan;
                % first choice
                indCaps=regexp(syl, '[A-Z]');
                if indCaps>2;
                    two_syl_back=lower(syl(indCaps-2));
                else
                    if syl_pos_in_motif>2;
                        % --- else try to use motif
                        two_syl_back=syls_in_motif{syl_pos_in_motif-2};
                    else
                        two_syl_back=nan;
                    end
                end
                
                
                if length(two_syl_back)>1;
                    tmp=regexp(two_syl_back,'[A-Z]');
                    two_syl_back=two_syl_back(tmp);
                else
                    two_syl_back=two_syl_back;
                end
                two_syl_back=lower(two_syl_back);
                
                
                % ============ is the preceding syl similar to preceding syl of target?
                if length(pre_syl)>1;
                    tmp=regexp(pre_syl,'[A-Z]');
                    pre_syl=pre_syl(tmp);
                else
                    pre_syl=pre_syl;
                end
                pre_syl=lower(pre_syl);
                
                if strcmpi(pre_syl,targsyl_pre); % if they are the same (this fails even if they are both nan)
                    presyl_similar_to_targ_presyl=1; % 1 means is similar
                else
                    presyl_similar_to_targ_presyl=0;
                end
                
                
                % what is the post syl?
                if syl_pos_in_motif==length(syls_in_motif);
                    % then this is last syl, no post syl
                    post_syl=nan;
                else
                    post_syl=syls_in_motif{syl_pos_in_motif+1};
                end
                
                
                % is post syl similar to target post syl?
                if length(post_syl)>1;
                    tmp=regexp(post_syl,'[A-Z]');
                    post_syl=lower(post_syl(tmp));
                else
                    post_syl=lower(post_syl);
                end
                if strcmpi(post_syl,targsyl_post); % if they are the same
                    postsyl_similar_to_targ_postsyl=1; % 1 means is similar
                else
                    postsyl_similar_to_targ_postsyl=0;
                end
                
                
                % if this is in same motif as targsyl, how many renditions
                % away is it? (i.e. -1 is preceding, +1 is post)
                if any(strcmp(syls_in_motif,targsyl));
                    % then this motif contains targ syl
                    targsyl_ind=find(strcmp(syls_in_motif,targsyl));
                    distance_from_targ=syl_pos_in_motif-targsyl_ind;
                else
                    distance_from_targ=nan;
                end
            end
            
            % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ CORRELATION STATS
            % SKIP FOR REPEAT EXPERIMENTS
            % NOTE: THIS CAN ONLY GET CORRS (motif or song) if both syls
            % are defined in regexp
            
            
            func_tmp=@(X) findstr(X, '+');
            if ~any(cell2mat(cellfun(func_tmp, SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.RegularExpressionsList, 'UniformOutput',0)));
                % then there is no repeat in the reg exp
                
                % ==================================================== To other syls in same motif
                if is_in_defined_motif==1;
                    if counter2==0;
                        counter2=1;
                        disp('NOTE: Correlations motif-to-motif: using subclass #1, assuming that is the important/only subclass');
                    end
                    
                    
                    
                    % -----------------------------------------------------------------
                    % figure out which class/subclass (strings) it is in
                    % -- compile motifs that are stereotyped.
                    motifs_potential={};
                    motifs_with_mult_subclasses={};
                    for k=1:length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.Params.RegExpr.subexpressions);
                        if length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.Params.RegExpr.subexpressions{k})>1; % skip, only works for things with one subexpression (e.g. abbccbb);
                            motifs_with_mult_subclasses=[motifs_with_mult_subclasses SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.Params.RegExpr.subexpressions{k}{1}];
                            
                            motifs_potential=[motifs_potential SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.Params.RegExpr.subexpressions{k}{1}];
                            
                            continue
                        end
                        
                        motifs_potential=[motifs_potential SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.Params.RegExpr.subexpressions{k}{1}];
                    end
                    
                    % tell user what threw out
                    for k=1:length(motifs_with_mult_subclasses);
                        
                        disp(['subclass motif correlations: had more than one subclass, took the first one for ' motifs_with_mult_subclasses{k}]);
                    end
                    % -----------------------------------------------------------------
                    
                    % Find which motif it is in, and which position in
                    % motif.
                    OUTPUT=Fn_FindSylInRegExp(motifs_potential, syl);
                    regexp_motifnum=OUTPUT.motif_of_origin;
                    regexp_PosInMotif_thissyl=OUTPUT.position_in_motif_of_origin;
                    
                    % --- if this syl is not in a motif, then skip it for
                    % motif by motif.
                    if ~isempty(regexp_motifnum)
                        
                        % -- Get list of other syls in motif
                        for k=1:length(syls_in_motif);
                            syl_other=syls_in_motif{k};
                            
                            % -- Get position of the other syl in the motif
                            OUTPUT=Fn_FindSylInRegExp(motifs_potential, syl_other);
                            if OUTPUT.motif_of_origin~=regexp_motifnum;
                                % then there is problem! they should be in same
                                % motif
                                disp('PROBLEM - other syl found to be in diff motif (regexp) than syl');
                                asdceagea; % halts
                            end
                            
                            regexp_PosInMotif_othersyl=OUTPUT.position_in_motif_of_origin;
                            
                            
                            % ====== COLLECT STATS
                            % Using motif-by-motif
                            corrcoeff=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.baseline.data_ParsedIntoSubclasses...
                                {regexp_motifnum}.sub_class{1}.CORRELATIONS.RhoMat(regexp_PosInMotif_thissyl, regexp_PosInMotif_othersyl);
                            pval=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.baseline.data_ParsedIntoSubclasses{regexp_motifnum}.sub_class{1}.CORRELATIONS.PvalMat(regexp_PosInMotif_thissyl, regexp_PosInMotif_othersyl);
                            
                            if isempty(corrcoeff)
                                corrcoeff=nan;
                            end
                            if isempty(pval)
                                pval=nan;
                            end
                            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.corrcoeff_vs.(syl_other)=corrcoeff;
                            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.p_val_vs.(syl_other)=pval;
                            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.regexp_motifnum=regexp_motifnum;
                            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.regexp_PosInMotif=regexp_PosInMotif_thissyl;
                            
                            % Using motif-by-motif, subtracting song mean
                            corrcoeff=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.baseline.data_ParsedIntoSubclasses{regexp_motifnum}.sub_class{1}.CORRELATIONS.SUBTRACT_SONG_MEAN.RhoMat(regexp_PosInMotif_thissyl, regexp_PosInMotif_othersyl);
                            pval=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.baseline.data_ParsedIntoSubclasses{regexp_motifnum}.sub_class{1}.CORRELATIONS.SUBTRACT_SONG_MEAN.PvalMat(regexp_PosInMotif_thissyl, regexp_PosInMotif_othersyl);
                            if isempty(corrcoeff)
                                corrcoeff=nan;
                            end
                            if isempty(pval)
                                pval=nan;
                            end
                            
                            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.SONG_MEAN_SUBTRACTED.corrcoeff_vs.(syl_other)=corrcoeff;
                            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.SONG_MEAN_SUBTRACTED.p_val_vs.(syl_other)=pval;
                            
                            % Using motif-by-motif, subtracting song mean, only
                            % keeping songs with >N motifs
                            corrcoeff=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.baseline.data_ParsedIntoSubclasses{regexp_motifnum}.sub_class{1}.CORRELATIONS.SUBTRACT_SONG_MEAN_NoShortSongs.RhoMat(regexp_PosInMotif_thissyl, regexp_PosInMotif_othersyl);
                            pval=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.baseline.data_ParsedIntoSubclasses{regexp_motifnum}.sub_class{1}.CORRELATIONS.SUBTRACT_SONG_MEAN_NoShortSongs.PvalMat(regexp_PosInMotif_thissyl, regexp_PosInMotif_othersyl);
                            if isempty(corrcoeff)
                                corrcoeff=nan;
                            end
                            if isempty(pval)
                                pval=nan;
                            end
                            
                            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.SONG_MEAN_SUBTRACTED_NoShortSong.corrcoeff_vs.(syl_other)=corrcoeff;
                            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.SONG_MEAN_SUBTRACTED_NoShortSong.p_val_vs.(syl_other)=pval;
                            
                        end
                    else
                        disp(['HEYHEYHEY - syl ' syl '(' birdname '-' exptname  ') is in defined motif (fields) but not in regexpr']);
                        for k=1:length(syls_in_motif);
                            syl_other=syls_in_motif{k};
                            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.corrcoeff_vs.(syl_other)=nan;
                            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.p_val_vs.(syl_other)=nan;
                            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.regexp_motifnum=nan;
                            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.regexp_PosInMotif=nan;
                            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.SONG_MEAN_SUBTRACTED.corrcoeff_vs.(syl_other)=nan;
                            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.SONG_MEAN_SUBTRACTED.p_val_vs.(syl_other)=nan;
                            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.SONG_MEAN_SUBTRACTED_NoShortSong.corrcoeff_vs.(syl_other)=nan;
                            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.SONG_MEAN_SUBTRACTED_NoShortSong.p_val_vs.(syl_other)=nan;
                        end
                    end
                    
                else
                    % not in a defined motif - note that down in the final
                    % structure
                    SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.motif_by_motif.NotInDefinedMotif=1;
                    regexp_motifnum=nan;
                    regexp_PosInMotif_thissyl=nan;
                end
                
                % ==================================================== To all syls in song (even in other motif)
                % String of order of correlations in matrix.
                if counter1==0;
                    disp('NOTE - Matching syls (e.g. dB) to position in correlation matrix only works if the regular expression strings are explicit (e.g. abbccb, and not [dj]bb, and have no subclasses');
                    counter1=1;
                end
                
                
                % ===== Find index of this syl (in the corr matrix)
                % where is the single syl? (upper case) (in motifs used for
                % regular expressions)
                motifs_potential=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.Params.RegExpr.expressions;
                
                OUTPUT=Fn_FindSylInRegExp(motifs_potential, syl);
                ind_of_this_syl=OUTPUT.ind_of_syl;
                
                % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if ~isempty(ind_of_this_syl)
                    for k=1:length(SylFields_Unique); % get corr coeff to all other syls.
                        
                        % === Find index of other syl
                        syl_other=SylFields_Unique{k};
                        
                        
                        OUTPUT=Fn_FindSylInRegExp(motifs_potential, syl_other);
                        ind_of_other_syl=OUTPUT.ind_of_syl;
                        
                        
                        % ================================= OUTPUT
                        corrcoeff=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.baseline.Correlations_Across_Classes.SONG_BY_SONG.RhoMat(ind_of_this_syl,ind_of_other_syl);
                        p_val=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.baseline.Correlations_Across_Classes.SONG_BY_SONG.PvalMat(ind_of_this_syl,ind_of_other_syl);
                        
                        if isempty(corrcoeff);
                            corrcoeff=nan;
                            p_val=nan;
                        end
                        
                        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.song_by_song.corrcoeff_vs.(syl_other)=corrcoeff;
                        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.song_by_song.p_val_vs.(syl_other)=p_val;
                        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.song_by_song.ind_of_this_syl=ind_of_this_syl;
                        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.song_by_song.ind_of_other_syl.(syl_other)=ind_of_other_syl;
                        
                    end
                    
                else
                    
                    for k=1:length(SylFields_Unique); % get corr coeff to all other syls.
                        
                        % === Find index of other syl
                        syl_other=SylFields_Unique{k};
                        
                        
                        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.song_by_song.corrcoeff_vs.(syl_other)=nan;
                        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.song_by_song.p_val_vs.(syl_other)=nan;
                        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.song_by_song.ind_of_this_syl=nan;
                        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).CORRELATIONS.song_by_song.ind_of_other_syl.(syl_other)=nan;
                        
                    end
                    
                    
                end
                % ------------------------------------------------------------------------------------------------------------
                
                
                SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).regexp_motifnum=regexp_motifnum;
                SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).regexp_PosInMotif_thissyl=regexp_PosInMotif_thissyl;
                
                
            else
                disp(['repeats in reg exp: no corr put in syl id for ' birdname '-' exptname]);
                
            end
            
            
            % ======================================= Put features into the original structure
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).single_syl=single_syl;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target=is_target;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ=similar_to_targ;
            
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).eucldist_from_targ_PCA=eucldist_from_targ_PCA;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).eucldist_from_targ_zscore=eucldist_from_targ_zscore;
            
            
            
            % motif related things (not regexp)
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).motif_num=motif_num;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).preceding_syl=pre_syl;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).two_syl_back=two_syl_back;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).post_syl=post_syl;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).distance_from_targ=distance_from_targ;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl=presyl_similar_to_targ_presyl;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).postsyl_similar_to_targ_postsyl=postsyl_similar_to_targ_postsyl;
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_in_defined_motif=is_in_defined_motif;
        end
    end
end



% ========================= TROUBLESHOOTING - looks through all syls,
% displays if not in motif, and subset, random, of syls that are in defined
% motif.
if (0)
    for i=1:NumBirds;
        expts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        for ii=1:expts;
            fields=fieldnames(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions);
            for iii=1:length(fields);
                syl=fields{iii};
                
                if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_in_defined_motif==0;
                    SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl)
                    SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}
                    
                end
                
                if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_in_defined_motif==1 && rand>0.95;
                    SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl)
                    SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}
                    
                    % SHOULD PUT - display all motifs and target syl
                    
                end
            end
        end
    end
end


%% DETERMINE IF 2 SYL BACK IS SIMILAR TO TARGET
for i=1:length(SeqDepPitch_AcrossBirds.birds);
    
    numexpt=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpt
        
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        
        two_syl_back_targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).two_syl_back;
        
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        for j=1:length(SylsUnique);
            syl=SylsUnique{j};
            
            two_syl_back=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).two_syl_back;
            
            if strcmp(two_syl_back, two_syl_back_targsyl);
                SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).two_syl_back_same_as_targ=1;
            else
                SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).two_syl_back_same_as_targ=0;
            end
        end
    end
end



%% GET INFORMATION ABOUT TARGET (i.e. learning at target) (start and end of WN)

for i=1:NumBirds;
    NumExperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:NumExperiments;
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.SylLists.TargetSyls{1};
        
        % WN start (1st 3 days)
        ind=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.AllDaysSliding.WindSize_3.(targsyl).FirstWNInd;
        
        try
            Y=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.AllDaysSliding.WindSize_3.(targsyl).meanFF_minusBaseline(ind);
        catch err
            % no WN day
            Y=nan;
        end
        % SLOT back into structure
        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Learning_by_target.WN_start_daybins=Y;
        
        % WN end
        if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.NumTargs==1;
            % only one target throughout
            ind=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.AllDaysSliding.WindSize_3.(targsyl).LastWNInd;
            
            try
                Y=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.AllDaysSliding.WindSize_3.(targsyl).meanFF_minusBaseline(ind);
            catch err
                Y=nan;
            end
        elseif SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.NumTargs==0;
            % one targ epoch, then 2;
            if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData,'Snapshot');
                tmpfield=fieldnames(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Snapshot);
                Y=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Snapshot.(tmpfield{1}).(targsyl).meanFF_minusBaseline;
            end
        end
        
        % SLOT back into structure
        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Learning_by_target.WN_end_daybins=Y;
    end
end

end



function OUTPUT=Fn_FindSylInRegExp(motifs_potential, syl)
% e.g.
% motifs_potential{1}='acbb';
% motifs_potential{2}='ddbbcc';
% syl= 'dB' (b after d)

% OUTPUT.ind_of_syl  - index if motifs_potnetial strings were concact (e.g
% in acbbddbbcc
% OUTPUT.motif_of_origin --- which ind of motifs_potential?
% OUTPUT.position_in_motif_of_origin    - ind in
% motifs_potential{OUTPUT.motif_of_origin}

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
syl_list_string=[];
for k=1:length(motifs_potential);
    syl_list_string=[syl_list_string '_' motifs_potential{k}];
end


% ===== Find index of this syl (in the corr matrix)
% where is the single syl? (upper case)
if length(syl)>1;
    location_of_upper_case=regexp(syl, '[A-Z]');
else
    location_of_upper_case=1; % just a single syl
end


% Find its position in the string that tells order of
% syls in correlation
ind_of_syl=regexp(syl_list_string, lower(syl));

% account for position of the upper case letter
% (desired syl)
ind_of_syl=ind_of_syl+location_of_upper_case-1;

% account for fact that I put a underscore for each motif in
% the syl_list_string
uscores_indices=strfind(syl_list_string(1:ind_of_syl), '_');
uscores_before_ind=length(uscores_indices);
ind_of_syl=ind_of_syl-uscores_before_ind;

% 1) which motif is this syl in? 2) what is position in
% that motif?
if uscores_before_ind==0; % then not in any regexp motif
    motif_of_origin=[];
    position_in_motif_of_origin=[];
else
    motif_of_origin=uscores_before_ind; % i.e. num uscores = which motif you are from.
    position_in_motif_of_origin=ind_of_syl+uscores_before_ind-uscores_indices(motif_of_origin);
end

% SANITY CHECK
if length(uscores_before_ind)>1;
    disp('PROBLEM, matched too many locations in correlations');
elseif isempty(uscores_before_ind);
    disp('PROBLEM, failed to match to correaltion string');
end

% -- output
OUTPUT.ind_of_syl=ind_of_syl;
OUTPUT.motif_of_origin=motif_of_origin;
OUTPUT.position_in_motif_of_origin=position_in_motif_of_origin;


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end

