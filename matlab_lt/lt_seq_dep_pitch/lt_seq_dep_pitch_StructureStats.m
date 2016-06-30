function [Params, AllDays_StructStatsStruct]=lt_seq_dep_pitch_StructureStats(Params, AllDays_RawDatStruct, DoLMAN, AllDays_PlotLearning)
%% LT 6/7/16 - power instead of amplitude
usePower=0; % for all measures using spectral amplitude, converted to power.
% note: using 1 actually got less separation (for what are clearly diff
% types - i.e. more clustered together). so use 0. not sure why using
% amplitude (0) is better...



%% LT 5/7/16 - Completed DoLMAN
% Saves a completely separate structure that only contains LMAN data. (if
% DoLMAN =0, then is as previous: only saves PBS data (is good: only within
% time window data if is LMAN experiment)



%% LT 12/29/15 - automatically detects whether this is LMAN experiment. If it is, only takes PBS data within time window (i.e. as determined from AllDays_PlotLearning structure)
% NEED TO INPUT BOTH STRUCTURES AllDays_RawDatStruct and
% AllDays_PlotLearning.  

%% LT 8/5/15 - LMAN musc data analysis added
% requires AllDays_PlotLearning to get the actual musc data within time window
% DoLMAN=1; then does both musc and pbs data in one run (default is 0)

            % NOTE: NOT DONE. Currently: if this is LMAN data, then this code DOES only take PBS
            % (but takes all PBS, not just in time window).
            
            


%% LT 4/28/15 - copied to "OLD" and continueing here with changes:
% director method of saving updated - now overwrites, and not saving
% timestamp in name.


%% LT 2/2/15 - takes raw data from lt_seq_dep_pitch and churns out stats on structure of syls.




%% PARAMS

if ~exist('DoLMAN', 'var');
    DoLMAN=0;
else
    % make sure has AllDays_PlotLearning
    
    if ~exist('AllDays_PlotLearning','var');
        disp('Problem, can"t do LMAN without AllDays_PlotLearning. exiting');
        dsafasfavs
    end     
end


SylFieldsAll=fieldnames(AllDays_RawDatStruct{1}.data); % all syls
PlotColors=lt_make_plot_colors(length(SylFieldsAll),0);
NumTargSyls=length(Params.SeqFilter.SylLists.TargetSyls); % number of targ syls
NumDays=datenum(Params.SeqFilter.LastDay)-datenum(Params.SeqFilter.FirstDay)+1;
fs=Params.DayRawDat.fs;

%% For similar syls, compare pitch
clear AllDays_StructStatsStruct

if DoLMAN==0
datafield='data';
elseif DoLMAN==1
    datafield='data_MUSC';
else
    disp('ERROR, what is DOLMAN?');
    dafcaearawcae;
end


%% BASELINE STUFF (OBSOLETE) (PITCH)
% 1) For all syls, get baseline pitch stats, combining all renditions
for i=1:length(SylFieldsAll);
    syl=SylFieldsAll{i};
    
    % 1) Baseline
    X=[];
    for ii=Params.SeqFilter.BaselineDays;
        try % in case day lacks data
            X=[X cell2mat(AllDays_RawDatStruct{ii}.(datafield).(syl)(:,1))']; % collect all data points from baseline days
        catch err
        end
    end
    
    AllDays_StructStatsStruct.IndivSyls.(syl).Baseline.valsFF=X;
    AllDays_StructStatsStruct.IndivSyls.(syl).Baseline.meanFF=mean(X);
    AllDays_StructStatsStruct.IndivSyls.(syl).Baseline.stdFF=std(X);
    AllDays_StructStatsStruct.IndivSyls.(syl).Baseline.n=length(X);
    AllDays_StructStatsStruct.IndivSyls.(syl).Baseline.sem=std(X)/sqrt((length(X)-1));
end


% 2) Plot
figure; hold on;
title('Mean pitch at baseline for all syls');

for i=1:length(SylFieldsAll);
    syl=SylFieldsAll{i};
    
    % color depending on if is target, similar, or diff syl
    if sum(strcmp(syl,Params.SeqFilter.SylLists.TargetSyls))==1;
        errorbar(i,AllDays_StructStatsStruct.IndivSyls.(syl).Baseline.meanFF,...
            AllDays_StructStatsStruct.IndivSyls.(syl).Baseline.sem, 'o','Color','k','MarkerSize',10);
    elseif sum(strcmp(syl,Params.SeqFilter.SylLists.SylsSame))==1;
        errorbar(i,AllDays_StructStatsStruct.IndivSyls.(syl).Baseline.meanFF,...
            AllDays_StructStatsStruct.IndivSyls.(syl).Baseline.sem, 'o','Color','g','MarkerSize',10);
    elseif  isfield(Params.SeqFilter.SylLists,'SylsDifferent');
        if sum(strcmp(syl,Params.SeqFilter.SylLists.SylsDifferent))==1;
            errorbar(i,AllDays_StructStatsStruct.IndivSyls.(syl).Baseline.meanFF,...
                AllDays_StructStatsStruct.IndivSyls.(syl).Baseline.sem, 'o','Color','r','MarkerSize',10);
        else
            errorbar(i,AllDays_StructStatsStruct.IndivSyls.(syl).Baseline.meanFF,...
                AllDays_StructStatsStruct.IndivSyls.(syl).Baseline.sem, 'o','Color','b','MarkerSize',10);
            disp(['Problem - ' syl ' is not targ, similar, or diff. what is it?']);
        end
    else
        errorbar(i,AllDays_StructStatsStruct.IndivSyls.(syl).Baseline.meanFF,...
            AllDays_StructStatsStruct.IndivSyls.(syl).Baseline.sem, 'o','Color','b','MarkerSize',10);
        disp(['Problem - ' syl ' is not targ, similar, or diff. what is it?']);
    end
end



set(gca,'XTick',1:length(SylFieldsAll));
set(gca,'XTickLabel',SylFieldsAll)



%% For all syls, get a vector of various features
% Calculate the vector for all renditions, so that can ask about
% trial-trial correlations as well as correlations across different
% duration bins.

ticID=tic;


for i=1:length(SylFieldsAll); % for each syllable
    syl=SylFieldsAll{i};
    
    c=1; % counter for rends for given syl
    FeatVect=[];
    for ii=1:NumDays; % for each day
        if ~isempty(AllDays_RawDatStruct{ii}); % chekc if day has data
            
            % check if day has syl
            if ~isfield(AllDays_RawDatStruct{ii}.(datafield), syl);
                continue; 
            end
            
            % check if this is LMAN experiment - if so, will get PBS data
            % in time window.
            if isfield(AllDays_RawDatStruct{ii}, 'data_MUSC');
                IsThisLMANExpt=1;
            else
                IsThisLMANExpt=0;
            end
            
            if IsThisLMANExpt==1
                % then use all raw dat, as some windows blurred
                datafield='data_WithOutlier';
            end
            
                NumRends=size(AllDays_RawDatStruct{ii}.(datafield).(syl),1); % how many renditions?
            
            for iii=1:NumRends; % for each rendition on that day
                            
                %  ===== IF THIS IS LMAN DATA THEN FOR EACH REND
                %  MUST CHECK WHETHER IT IS WITHIN TIME WINDOW, USING PLOT
                %  LEARNING STRUCTURE
                if IsThisLMANExpt==1;
                    TimeOfSong=AllDays_RawDatStruct{ii}.(datafield).(syl){iii, 6};
                    
                    if DoLMAN==1
                        % then check this rend vs. DataMatrix_MUSC
                        if ~any(AllDays_PlotLearning.DataMatrix_MUSC.(syl).Tvals_WithinTimeWindow{ii}==TimeOfSong);
                            % then this rend is not within time window;
                            disp(['threw out day ' num2str(ii) '; rend: ' num2str(iii) '/' num2str(NumRends) ', because not in LMAN time window']);
                            continue;
%                         else
%                             disp(['kep day ' num2str(ii) '; rend: ' num2str(iii) '/' num2str(NumRends)])
                        end
                        
                    else
                        % then check against DataMatrix
                        if ~any(AllDays_PlotLearning.DataMatrix.(syl).Tvals_WithinTimeWindow{ii}==TimeOfSong);
                            % then this rend is not within time window;
                            disp(['threw out day ' num2str(ii) '; rend: ' num2str(iii) '/' num2str(NumRends) ', because not in LMAN time window']);
                            continue;
                        end
                        
                    end
                end
                    
                    % =====
                    % 1) calculate feature vector for this rendition - based on
                    % Wohlgemuth et al., 2010, and Sakata 2006
                    
                    % First, get raw sound data and smooth it/bandpass it
                    Dat=AllDays_RawDatStruct{ii}.(datafield).(syl){iii,3}; % raw sound data
                % Process
                Dat=bandpass(Dat,fs,500,10000,'hanningfir'); % bandpass
                Dat_sq=Dat.^2; % squared
                
                % smooth
                sm_win=2; % ms
                len = round(Params.DayRawDat.fs*sm_win/1000); % 2ms smooth wind
                h   = ones(1,len)/len; % uniform kernel
                Dat_sqsm = conv(h, Dat_sq);
                offset = round((length(Dat_sqsm)-length(Dat))/2); % remove offset from conv.
                Dat_sqsm=Dat_sqsm(1+offset:length(Dat)+offset);
                
                
                % FIND ONSET - ampl threshold already defined in params
                % note: Onset is about 13ms from start of data - but I put
                % 5ms preceding data start in findwnote to get the
                % data originally.  why is that?
                % Look inbetween bins 410 422. empriicalyl it is always there.
                % this removes possibility that will hit onset at some point
                % before those bins (i.e. a spike before onset).
                tmp=find(Params.SeqFilter.AmplThr-Dat_sqsm(410:422)<0,1,'first'); % up cross of ampl threshold
                if isempty(tmp); % its possible didnt find crossing in that window- e.g. hand segmented syl
                    tmp=find(Params.SeqFilter.AmplThr-Dat_sqsm(315:510)<0,1,'first'); % give a +/- 3 ms window to find crossing
                    if isempty(tmp);
                        % then throw out this rendition.
                        disp(['Problem for syl ' syl ' on day ' num2str(ii) ', rendition ' num2str(iii) ' - cant find onset, will throw out']);
                        continue
                    else
                        onset=355+tmp-1;
                    end
                else
                    onset=410+tmp-1;
                end
                AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.Info.onset(c)=onset;
                
                
                % EXTRACT SYL WINDOWED sound data (removing pre and post stuff);
                Dur=AllDays_RawDatStruct{ii}.(datafield).(syl){iii,10}; % syl dur
                if onset+floor(fs*Dur/1000)>length(Dat); % problem - not enough data
                    disp(['Problem for syl ' syl ' on day ' num2str(ii) ', rendition ' num2str(iii) ' - raw dat too short, will throw out.']);
                    continue
                end
                Dat_syl=Dat(onset:onset+floor(fs*Dur/1000)); % raw dat
                Dat_sqsm_syl=Dat_sqsm(onset:onset+floor(fs*Dur/1000)); % work with this from now on
                Dat_sqsm_syl_norm=Dat_sqsm_syl/sum(Dat_sqsm_syl); % normalize
                
                
                % EXTRACT FEATURES --------------------------------------
                
                % DURATION
                Dur=AllDays_RawDatStruct{ii}.(datafield).(syl){iii,10}; % in ms
                
                
                % FF
                FF=AllDays_RawDatStruct{ii}.(datafield).(syl){iii,1}; % FF calculated based flat pitch window
                
                
                % SPECTRAL ENTROPY
                % Get Spectrum of entire syllable
                [Spectrum,Fsp, ~]=spectrogram(Dat_syl,hanning(length(Dat_syl)),1,length(Dat_syl),fs);
                % remove unwanted freqs
                Spectrum=Spectrum(Fsp<8000);
                Fsp=Fsp(Fsp<8000);
                Spectrum=Spectrum(Fsp>300);
                Fsp=Fsp(Fsp>300);
                % get norm
                Spectrum=abs(Spectrum);
                % --- square to get power
                if (usePower==1)
                   Spectrum=Spectrum.^2; 
                end
                % Normalize
                Spectrum=Spectrum/sum(Spectrum);
                
                % Spectral Entropy
                SpecEntr=-(Spectrum'*log2(Spectrum))/log2(length(Spectrum)); % 0 to 1, pure to WN
                
                
                % AMPLITUDE ENTROPY (over time)
                AmplEntr=-(Dat_sqsm_syl_norm'*log2(Dat_sqsm_syl_norm))/log2(length(Dat_sqsm_syl_norm)); % 0 to 1, temporal peak to flat
                
                
                % SPECTROTEMPORAL ENTROPY
                nfft=256; % 8ms
                olap=0; % no overlap
                
                [Sp, Fsp, ~]=spectrogram(Dat_syl,nfft,olap,nfft,fs);
                Sp=abs(Sp);
                % filter out unwanted frequencies
                Sp=Sp(Fsp<8000,:);
                Fsp=Fsp(Fsp<8000);
                
                Sp=Sp(Fsp>300,:);
                Fsp=Fsp(Fsp>300);
                
                % --- square to get power
                if (usePower==1)
                Sp=Sp.^2;
                end
                
                Sp=Sp/sum(sum(Sp)); % normalize
                
                % Get entropy
                STentropy=-sum(sum(Sp.*log2(Sp)))/sum(sum(log2(numel(Sp))));
                
                
                % AMPLITUDE AND FREQUENCY "SLOPE" (half2 - half2)/(sum)
                % spectrum of each half.
                [Sp, Fsp, ~]=spectrogram(Dat_syl,ceil(length(Dat_syl)/2),1,ceil(length(Dat_syl)/2),fs);
                
                % filter out high and low freq
                Sp=Sp(Fsp<8000,:);
                Fsp=Fsp(Fsp<8000);
                Sp=Sp(Fsp>300,:);
                Fsp=Fsp(Fsp>300);
                
                % take halves and absolute
                sd1=abs(Sp(:,1));
                sd2=abs(Sp(:,2));
                
                if usePower==1
                    sd1=sd1.^2;
                    sd2=sd2.^2;
                end
                
                % Amplitude modulation
                sum1=sum(sd1);
                sum2=sum(sd2);
                AmplSlope=(sum2-sum1)/(sum2+sum1);
                
                % Frequency slope
                f1=sum(sd1.*Fsp)/sum1;
                f2=sum(sd2.*Fsp)/sum2;
                FreqSlope=(f2-f1)/(f2+f1);
                
                
                % TIME TO HALF PEAK AMPL (as fraction of dur)
                [peak, loc]=max(Dat_sqsm_syl_norm);
                
                halftime=find(Dat_sqsm_syl_norm-peak/2>0,1,'first'); % index where sound reaches half max ampl
                if halftime==1; % i.e. onset is already higher than half max - problem
                    disp(['Problem - onset already higher than half max - syl ' syl ]);
                end
                halftime=halftime/(fs*Dur/1000); % time to half max, as fraction of syl dur
                
                
                % DONE GETTING FEATURES --------------------------------------
                
                % ASSIGN ALL to feature vector
                
                FeatVect(c,1)=Dur;
                FeatVect(c,2)=FF;
                FeatVect(c,3)=SpecEntr;
                FeatVect(c,4)=AmplEntr;
                FeatVect(c,5)=STentropy;
                FeatVect(c,6)=AmplSlope;
                FeatVect(c,7)=FreqSlope;
                FeatVect(c,8)=halftime;
                
                
                % 2) extract desired raw data
                AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.Info.datenum(c)=AllDays_RawDatStruct{ii}.(datafield).(syl){iii,6}; % datenum
                AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.Info.dayInd(c)=ii; % day Index
                AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.Info.filename{c}=AllDays_RawDatStruct{ii}.(datafield).(syl){iii,5}; % filename
                AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.Info.labels{c}=AllDays_RawDatStruct{ii}.(datafield).(syl){iii,7}; % filename
                AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.Info.notepos(c)=AllDays_RawDatStruct{ii}.(datafield).(syl){iii,8}; % filename
                AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.Info.trig(c)=AllDays_RawDatStruct{ii}.(datafield).(syl){iii,9}; % filename
                
                

                c=c+1; % iterate to next rend
            end
            
        end
    end
    AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.FeatVect=FeatVect;
    
end

Params.StructureStats.FeatureLegend={'Dur','FF','SpecEntr','AmplEntr','STentropy','AmplSlope','FreqSlope','halftime'};

RunTime=toc(ticID);




%-----------------------------------------
% troubleshooting plot all waveforms and onsets
% figure; hold on;
% plot(Dat_tmp')
% for j=1:length(onset(i,1,:));
%     line([onset(i,1,j) onset(i,1,j)],ylim);
% end
%
%
%      for i=1:length(SylFieldsAll);
% std(onset(i,ii,:));
% disp('');
% mean(onset(i,ii,:))
% end

% Duration - first get amplitude threshold. I know data starts -5ms before
% the amplitude I chose using evsonganaly, so set threshold to that value

%             figure;
%             plot(Dat);
%             line([x x],ylim)
%             figure; plot(Dat_sq);
%                         line([x x],ylim)
%
%             figure;  plot(Dat_sm);
%             line([x x],ylim)
%             figure;  plot(Dat_sm_log);
%             line([x x],ylim)
%-----------------------------------------        END TROUBLESHOOT




%% PLOT FEATURE VALUES TO VISUALIZE
if (0)
    % NOTE: issue: for MUSC experimesnt (PBS or MUSC data) takes from with
    % outliers, so example might not be actual data
SylFieldsAll=fieldnames(AllDays_StructStatsStruct.IndivSyls);
fs=Params.DayRawDat.fs;

for i=1:length(SylFieldsAll);
    syl=SylFieldsAll{i};
    
    FV=AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.FeatVect;
    
    if isempty(FV);
        continue;
    end
       
    if isempty(AllDays_RawDatStruct{1}.(datafield).(syl));
        continue; 
    end
    
    figure; hold on;
    
    % PLOT example amplitude profile and spectrogram
    % First, amplitude profile and related features
    subplot(1,2,1); hold on;
    title('Example amplitude contour');
    
    Dur=AllDays_RawDatStruct{1}.(datafield).(syl){1,10}; % in ms   
    Dat=AllDays_RawDatStruct{1}.(datafield).(syl){1,3}; % choose the first rendition as example
    Dat=bandpass(Dat,fs,500,10000,'hanningfir'); % bandpass
    Dat_sq=Dat.^2; % squared
    
    sm_win=2; % ms
    len = round(Params.DayRawDat.fs*sm_win/1000); % 2ms smooth wind
    h   = ones(1,len)/len; % uniform kernel
    
    Dat_sqsm = conv(h, Dat_sq);
    offset = round((length(Dat_sqsm)-length(Dat))/2); % remove offset from conv.
    Dat_sqsm=Dat_sqsm(1+offset:length(Dat)+offset);
    onset=find(Params.SeqFilter.AmplThr-Dat_sqsm<0,1,'first'); % up cross of ampl threshold
    
    try % happen to choose a rendition without enough data ...
        Dat_sqsm_syl=Dat_sqsm(onset:onset+floor(fs*Dur/1000)); % work with this from now on
    catch err
                    disp(['syl ' syl ' 1st examples are not long enough, skipping']);
                    disp(' ');
            continue
% 
%         try
%             % choose next syl
%             Dur=AllDays_RawDatStruct{1}.data.(syl){1,10}; % in ms
%             Dat=AllDays_RawDatStruct{1}.data.(syl){1,3}; % choose the first rendition as example
%             Dat=bandpass(Dat,fs,500,10000,'hanningfir'); % bandpass
%             Dat_sq=Dat.^2; % squared
%             
%             Dat_sqsm = conv(h, Dat_sq);
%             offset = round((length(Dat_sqsm)-length(Dat))/2); % remove offset from conv.
%             Dat_sqsm=Dat_sqsm(1+offset:length(Dat)+offset);
%             onset=find(Params.SeqFilter.AmplThr-Dat_sqsm<0,1,'first'); % up cross of ampl threshold
%         catch err
%         end
    end
    
    x=1/fs:1/fs:length(Dat_sqsm_syl)/fs;
    plot(x,Dat_sqsm_syl);
    
    % overlay amplitude related features
    
    meanFV=mean(FV,1);
    stdFV=std(FV,1);
    
    
    % Second, spectrogram and related features
    subplot(1,2,2); hold on;
    title('Example Spectrogram');
    
    nfft=256; % 8ms
    olap=0; % no overlap
    Dat_syl=Dat(onset:onset+floor(fs*Dur/1000)); % raw dat
    [Sp, Fsp, Tsp]=spectrogram(Dat_syl,nfft,floor(0.9*nfft),nfft,fs);
    Sp=abs(Sp);
    Sp=log(Sp);
    
    imagesc(Tsp,Fsp,Sp);
    ylim([300 8000]);
    
    subtitle(['Syllable: ' syl ' (first day)']);
    
    % List calculated features
    
    legend=Params.StructureStats.FeatureLegend;
    disp(['SYLLABLE: ' syl '; ' num2str(size(FV,1)) ' samples.']);
    for ii=1:length(meanFV);
        disp([legend{ii} ' = ' num2str(meanFV(ii)) ' +/- ' num2str(stdFV(ii))]);
    end
    disp(' ');
    
end

end


%% List sample sizes for all syls

disp('Baseline sample size for:');
for i=1:length(SylFieldsAll);
    
    syl=SylFieldsAll{i};
    
    %     disp([syl ': ' num2str(length(AllDays_StructStatsStruct.IndivSyls.(syl).AllRends.Info.datenum))]);
    disp([syl ': ' num2str((AllDays_StructStatsStruct.IndivSyls.(syl).Baseline.n))]);
    
end



%% Figure out which data are for baseline



%% SAVE
tstamp=lt_get_timestamp(0);
Params.StructureStats.savedir=Params.SeqFilter.savedir; % obsolete

cd(Params.SeqFilter.savedir);

if DoLMAN==1
    save('AllDays_StructStatsStruct_MUSC','AllDays_StructStatsStruct');
    save('Params.mat','Params');
    
    %     write a text file that tells you when files were made
    fid1=fopen(['DONE_StructureStatsMUSC_' tstamp '.txt'],'w');
    fclose(fid1);
    
    
    try
        cd FIGURES/StructureStats_MUSC
    catch err
        mkdir FIGURES/StructureStats_MUSC
        cd FIGURES/StructureStats_MUSC
    end
    
    lt_save_all_figs
else
    save('AllDays_StructStatsStruct','AllDays_StructStatsStruct');
    save('Params.mat','Params');
    
    %     write a text file that tells you when files were made
    fid1=fopen(['DONE_StructureStats_' tstamp '.txt'],'w');
    fclose(fid1);
    
    
    try
        cd FIGURES/StructureStats
    catch err
        mkdir FIGURES/StructureStats
        cd FIGURES/StructureStats
    end
    
    lt_save_all_figs
    
    
end


cd ../../






