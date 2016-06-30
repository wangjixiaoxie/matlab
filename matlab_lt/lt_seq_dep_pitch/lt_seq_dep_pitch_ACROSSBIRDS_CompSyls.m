function SeqDepPitch_AcrossBirds=lt_seq_dep_pitch_ACROSSBIRDS_CompSyls(SeqDepPitch_AcrossBirds)
%% LT 1/12/16 - compiles 4 classes (type x presim) into cells

for i=1:length(SeqDepPitch_AcrossBirds.birds);
    
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpts;
    
        
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        % === collect 
        STSS={};
        STDS={};
        DTSS={};
        DTDS={};
        
        for j=1:length(SylsUnique);
            syl=SylsUnique{j};
            
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==1;
                continue;
            end
            
            similar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
            presimilar=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl;
                
            if similar & presimilar
                STSS=[STSS syl];
            elseif similar & ~presimilar
                STDS=[STDS syl];
            elseif ~similar & presimilar
                DTSS=[DTSS syl];
            elseif ~similar & ~presimilar
                DTDS=[DTDS syl];
            else
                disp('WARNING - syl not classified');
                asdcascef
            end
        end
        
        % ==== store for this expt
        
        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_STSS=STSS;
        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_STDS=STDS;
        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_DTSS=DTSS;
        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_DTDS=DTDS;
            
        
        
    end
end