function lt_seq_dep_pitch_ACROSSBIRDS_PlotSpec2(SeqDepPitch_AcrossBirds, PARAMS, plotForIllustrator)
%% LT 1/12/16 - uses v1 multiple times to plot all bird/expts of choice

%% === global params


%% specific bird/expt
if ~isempty(PARAMS.PlotSpec.BirdToPlot);
    if isempty(PARAMS.PlotSpec.ExptToPlot);
        % then plot all expts
        for i=1:length(SeqDepPitch_AcrossBirds.birds);
            birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
            
            if ~strcmp(PARAMS.PlotSpec.BirdToPlot, birdname);
                continue
            end
                
            numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
            for ii=1:numexpts
                
                exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
                
                PARAMS.PlotSpec.ExptToPlot=exptname;
                PARAMS.PlotSpec.MotifsToPlot=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.RegularExpressionsList;
                
                lt_seq_dep_pitch_ACROSSBIRDS_PlotSpec(SeqDepPitch_AcrossBirds, PARAMS, plotForIllustrator);
            end
        end
    else
        for i=1:length(SeqDepPitch_AcrossBirds.birds);
            birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
            
            
            numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
            for ii=1:numexpts
                
                exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
                
                if ~strcmp(PARAMS.PlotSpec.BirdToPlot, birdname) | ~strcmp(PARAMS.PlotSpec.ExptToPlot, exptname)
                    continue
                end
                
                PARAMS.PlotSpec.MotifsToPlot=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.RegularExpressionsList;
                
                lt_seq_dep_pitch_ACROSSBIRDS_PlotSpec(SeqDepPitch_AcrossBirds, PARAMS, plotForIllustrator);
            end
        end
    end
else
    %% go thru each bird (plot all birds)
    
    for i=1:length(SeqDepPitch_AcrossBirds.birds);
        birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
        
        numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
        
        for ii=1:numexpts
            
            exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
            
            
            PARAMS.PlotSpec.BirdToPlot=birdname;
            PARAMS.PlotSpec.ExptToPlot=exptname;
            
            PARAMS.PlotSpec.MotifsToPlot=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.RegularExpressionsList;
            
            lt_seq_dep_pitch_ACROSSBIRDS_PlotSpec(SeqDepPitch_AcrossBirds, PARAMS, plotForIllustrator);
            
            
        end
    end
end