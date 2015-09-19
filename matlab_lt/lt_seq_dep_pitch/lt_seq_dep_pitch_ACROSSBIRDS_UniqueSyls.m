function [SeqDepPitch_AcrossBirds]=lt_seq_dep_pitch_ACROSSBIRDS_UniqueSyls(SeqDepPitch_AcrossBirds)
        
NumBirds=length(SeqDepPitch_AcrossBirds.birds);

%% CHECK ALL EXPERIMENTS WHETHER CONTAIN UYNQIUE SYLS FIELD
for j=1:NumBirds;
    NumExperiments=length(SeqDepPitch_AcrossBirds.birds{j}.experiment);
    
    for jj=1:NumExperiments;
        if ~isfield(SeqDepPitch_AcrossBirds.birds{j}.experiment{jj}.Data_PlotLearning.Params.PlotLearning, 'SylFields_Unique');
            
            disp(['Missing unique syls for bird: ' num2str(j), ' expt: ' num2str(jj)]);
            
            % ==== then missing field. create it now
            Params=SeqDepPitch_AcrossBirds.birds{j}.experiment{jj}.Data_PlotLearning.Params;
            SylFieldsSingle=SeqDepPitch_AcrossBirds.birds{j}.experiment{jj}.Data_PlotLearning.Params.PlotLearning.SylFieldsSingle;
            
            % ========== Get relevant syls - i.e. only syls uniquely defined by sequence context
            % (i.e. don't double count b and a[b]);
            SylFields_Unique={};
            for i=1:length(Params.SeqFilter.SylLists.FieldsInOrder)
                SylFields_Unique=[SylFields_Unique Params.SeqFilter.SylLists.FieldsInOrder{i}];
            end
            
            % some syls are unique but not specified in sequence.  extract those as
            % well.
            UpperSyls=regexp(SylFields_Unique, '[A-Z]'); % find all upper case syls alrady in unique syls. if any single syls are not represented, then put the single syl into unique syls
            upper_syls_lowered={};
            for i=1:length(UpperSyls);
                if isempty(UpperSyls{i});
                    upper_syls_lowered{i}=SylFields_Unique{i};
                else
                    upper_syls_lowered{i}=lower(SylFields_Unique{i}(UpperSyls{i}));
                end
            end
            
            % for all single syls, if not represented in lower syls (i.e. syls that
            % are already called unqiue), then add it to the unique group
            for i=1:length(SylFieldsSingle);
                % only continue if the single syl is lowercase (upper is redundant
                % as has to have a lower case of same thing)
                if regexp(SylFieldsSingle{i}, '[A-K]');
                    continue
                end
                
                if ~any(strcmp(SylFieldsSingle{i}, upper_syls_lowered));
                    SylFields_Unique=[SylFields_Unique SylFieldsSingle{i}];
                end
            end
            % ====================================================
            
            % ===================== PUT BACK INTO OUTPUT
            SeqDepPitch_AcrossBirds.birds{j}.experiment{jj}.Data_PlotLearning.Params.PlotLearning.SylFields_Unique=SylFields_Unique;
            SeqDepPitch_AcrossBirds.birds{j}.experiment{jj}.Data_PlotLearning.Params.PlotLearning.SylFields_Unique=SylFields_Unique;
            
        end
        
        % ==== ANOTHER OUTPUT
            SeqDepPitch_AcrossBirds.birds{j}.experiment{jj}.INFORMATION.SylFields_Unique=SeqDepPitch_AcrossBirds.birds{j}.experiment{jj}.Data_PlotLearning.Params.PlotLearning.SylFields_Unique;
            SeqDepPitch_AcrossBirds.birds{j}.experiment{jj}.INFORMATION.SylFields_Unique=SeqDepPitch_AcrossBirds.birds{j}.experiment{jj}.Data_PlotLearning.Params.PlotLearning.SylFields_Unique;
        
    end
end




        
        

