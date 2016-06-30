function [SeqDepPitch_AcrossBirds] = lt_seq_dep_pitch_ACROSSBIRDS_RemoveSyls(SeqDepPitch_AcrossBirds, Params);
%% LT 9/7/15 - Remove syls by eye


NumBirds=length(SeqDepPitch_AcrossBirds.birds);

%%
for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        
        % Check whether this birdname expt pair has syl that should be removed.
        % If so, remove them from unique syls
        
        inds=find(strcmp(Params.global.SylsToRemove, birdname));
        
        for j=1:length(inds);
            
            expt_toremove=Params.global.SylsToRemove{inds(j)+1};
            syls_toremove=Params.global.SylsToRemove{inds(j)+2};
            
            % IF CURRENT EXPERIMENT IS THE ONE TO REMOVE, THEN DO SO
            if strcmp(exptname, expt_toremove);
                
                for k=1:length(syls_toremove);
                    
                    tmp_sylremove=syls_toremove{k};
                    
                    syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
                    
                    ind_to_remove=strcmp(syls_unique, tmp_sylremove);
                    
                    if isempty(ind_to_remove)
                        disp([birdname '-' exptname '-' tmp_sylremove '(syl to remove) is not actually in SylsUnique']);
                    end
                    
                    syls_unique(ind_to_remove)=[];
                    
                    % Put back into structure
                    SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique=syls_unique;
%                     SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.SylFields_Unique=syls_unique;
                    
                    
                    % tell user this is done
                    disp([birdname '-' exptname ': removed ' tmp_sylremove]);
                end
            end
        end
    end
end
