%% modify test part of code below - displays name of syls that pass criteria.

for i=1:length(SeqDepPitch_AcrossBirds.birds);
    
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpts
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        
        for j=1:length(SylsUnique);
            
            syl=SylsUnique{j};
            
            
            
            % ======= TEST IF IS THESE TRAITS
            
            
            if                     SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl==0 & ...
                    SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==0;
%                     SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).two_syl_back_same_as_targ==0 & ...
                
                
                disp([birdname '-' exptname '-' syl]);
            end
                        
            
        end
    end
end