function [SeqDepPitch_AcrossBirds, PARAMS]=lt_seq_dep_pitch_ACROSSBIRDS_LMANrenorm(SeqDepPitch_AcrossBirds, PARAMS, musc_day_window_WNdayinds, musc_day_window_Bidirdayinds_Bidir, pu53_use_later_days, debugON)
% musc_day_window_WNdayinds=[3 10]; % from WN day 3 to 10, collect all for analysis
% pu53_use_later_days=1; then for (birdname, 'pu53wh88') & strcmp(exptname, 'SeqDepPitchLMAN'), uses later days, because inactivation did not work for earlier days

debug=0;
disp('NOTE: only modifies within time window data');

%% 1) SORT OUT ONLY THE THE EXPEIRMENTS THAT HAVE LMAN INACTIVATION DATA
% copy strcuture, save backup.
filter = 'LMAN';
[SeqDepPitch_AcrossBirds, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds, filter);


%% == GO thru each expt, determine pitch for each day relative to PBS baseline (instead of relative to MUSC baseline as before)

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        disp([' == ' birdname ' - ' exptname]);
        for j=1:length(SylsUnique)
            syl=SylsUnique{j};
            NumDays=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).meanFF);
            
            % --- extract syl stuff
                baseline_pbs=SeqDepPitch_AcrossBirds.birds{i}.experiment...
                    {ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).meanFF_WithinTimeWindow;
                baseline_musc=SeqDepPitch_AcrossBirds.birds{i}.experiment...
                    {ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData_MUSC.Baseline.(syl).meanFF_WithinTimeWindow;
                baseline_musc_minus_pbs=baseline_musc-baseline_pbs;
                
                if debug==1;
                disp(['vals below should be: ' num2str(baseline_musc_minus_pbs)]);
                end
            for day=1:NumDays;
                
                if isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_WithinTimeWindow{day});
                    continue
                end
                
                
                
                % ==== METHOD 2
                ffvals_raw=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_WithinTimeWindow{day};
                ffvals_minusBase_new=ffvals_raw-baseline_pbs;
                
             
                % sanity check
                ffvals_minusBase_old=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_DevFromBase_WithinTimeWindow{day};
                if (mean(ffvals_minusBase_new)-mean(ffvals_minusBase_old)) - baseline_musc_minus_pbs>0.1
%                 disp([syl ' - day ' num2str(day) ': ' num2str(mean(ffvals_minusBase_new)-mean(ffvals_minusBase_old))]);
                disp([syl ' - day ' num2str(day) ': ' num2str(mean(ffvals_minusBase_new)-mean(ffvals_minusBase_old)-baseline_musc_minus_pbs)]);
                
                end
 
                % ---- OUTPUT
                SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_WithinTimeWindow{day}=ffvals_minusBase_new;
                SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).meanFF_DevFromBase_WithinTimeWindow(day)=mean(ffvals_minusBase_new);
                
                
                
%                 % === OLD METHOD, not best (derived from difeference in
%                 % baseline means)
%                 % --- change day ffvals to ffvals minus either pbs or musc
%                 % baseline pitch
%                 ffvals_minusBase_old=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_DevFromBase_WithinTimeWindow{day};
%                 
%                 % adjust pitch based on difference of MUSC baseline from
%                 % PBS baseline
%                 
%                 ffvals_minusBase_new = ffvals_minusBase_old + baseline_musc_minus_pbs;
%                 
%                 % --- troubleshooting: compare new ffval (based on
%                 % baseline) to actual ffval
%                 ffvals_raw_new=ffvals_minusBase_new + baseline_pbs;
%                 ffvals_actual=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(syl).FFvals_WithinTimeWindow{day};
%                 
%                 if mean(ffvals_raw_new) ~= mean(ffvals_actual)
%                 disp([syl ' - day ' num2str(day) ': ' num2str(mean(ffvals_raw_new)) ' -- ' num2str(mean(ffvals_actual))]);
%                 end
                
            end
            
        end
    end
end

