function [SeqDepPitch_AcrossBirds, Params]=lt_seq_dep_pitch_ACROSSBIRDS_AdHocChanges(SeqDepPitch_AcrossBirds, Params);


%% ===== 6/22/16 [ relate to Sober analysis] AD HOC, EXPT I HAVE NOT YET FINISHED LABELING, BUT DELAY IS
        % WRONG
 NumBirds =length(SeqDepPitch_AcrossBirds.birds);
 
for i=1:NumBirds
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpts
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
            
        if strcmp(birdname, 'or100pu10') & strcmp(exptname, 'SeqDepPitchLMAN2')
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.NumEmptyDays_StartWN_FromZscoreCode=0;
                disp(['daysDelay to 0 for: ' birdname '-' exptname]);
        end
            
        if strcmp(birdname, 'gr87bu18') & strcmp(exptname, 'SeqDepPitch2')
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.NumEmptyDays_StartWN_FromZscoreCode=0;
                disp(['daysDelay to 0 for: ' birdname '-' exptname]);
        end
        
    end
end
        
