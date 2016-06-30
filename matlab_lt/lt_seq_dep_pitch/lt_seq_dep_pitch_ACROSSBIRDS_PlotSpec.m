function lt_seq_dep_pitch_ACROSSBIRDS_PlotSpec(SeqDepPitch_AcrossBirds, PARAMS, plotForIllustrator)
%% LT 1/12/16 

% Note: gets specs from last baseline day.
% NOTE: optimal if DayRawDat has data field for syl 'b', as uses that to
% get list of song labels. Otherwise it picks the first on the fieldnames
% list for that DayRawDat.

%%
BirdToPlot=PARAMS.PlotSpec.BirdToPlot;
ExptToPlot=PARAMS.PlotSpec.ExptToPlot;

% --- motifs
NumRendsToPlot=PARAMS.PlotSpec.NumRendsToPlot;
MotifsToPlot=PARAMS.PlotSpec.MotifsToPlot;


% --- syls
NumPreNotes=PARAMS.PlotSpec.NumPreNotes;
NumPostNotes=PARAMS.PlotSpec.NumPostNotes;
PlotAcousticDist=PARAMS.PlotSpec.PlotAcousticDist;
PlotGeneralization=PARAMS.PlotSpec.PlotGeneralization;


% --- both
PreDur=0.048; %
PostDur=0.048;



NumBirds=length(SeqDepPitch_AcrossBirds.birds);

%% MOTIFS
if NumRendsToPlot>0;
for i=1:NumBirds;
    
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpts;
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
    
        if ~strcmp(birdname, BirdToPlot) | ~strcmp(exptname,ExptToPlot);
        continue;
        end
        
        figcount=1;
subplotrows=6;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

if plotForIllustrator==1;
    subplotrows=2;
    subplotcols=1;
end

hsplotsAll=[];

        % ====
            
            WNday1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
            FirstDay_date=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.FirstDay;
            FirstDay_datenum=datenum(FirstDay_date, 'ddmmmyyyy');
            LastBaselineDay_datenum=FirstDay_datenum+WNday1-2;
            LastBaselineDay_date=datestr(LastBaselineDay_datenum, 'ddmmmyyyy');
            
            cd(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.PlotLearningDir);
            cd ..
            
            load(['RawDatStruct_' LastBaselineDay_date '.mat']);
            ParamsDay=load(['Params_' LastBaselineDay_date '.mat']);
            
            % === collect labels
            sylfields=fieldnames(RawDatStruct.data);
            
            if any(strcmp(sylfields, 'b'));
                SylToUse='b';
            else
                SylToUse=sylfields{1};
            end
            LabelsAll=RawDatStruct.data.(SylToUse)(:, 7);
            

        %% ==== plot multiple renditions of the desired motif
        
        
        
        for j=1:length(MotifsToPlot);
            motif=MotifsToPlot{j};

            for jj=1:NumRendsToPlot;
                
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title(motif);
            hsplotsAll=[hsplotsAll hsplot];
            
            NoteInd=[];
            ccc=1;
            while isempty(NoteInd);
                
                ind=randi(length(LabelsAll), 1); % pick a random song
                [NoteInd, MatchString]=regexp(LabelsAll{ind}, motif, 'start','match'); % get matches for this motif
                ccc=ccc+1;
                if ccc==20;
                    disp(['PROBLEM - CANT FIND SONG WITH THIS MOTIF ' motif]);
                    break
                end
            end
            
            indtmp=randi(length(NoteInd));
            NoteInd=NoteInd(indtmp); % pick a random note
            NoteIndRange=NoteInd:NoteInd+length(MatchString{indtmp})-1;  % range of note inds for entire motif    
            
            % === extract sound data for that note range and plot
            fnameNotMat=RawDatStruct.data.(SylToUse){ind, 5};
            fnameCbin=fnameNotMat(1:strfind(fnameNotMat, '.not.mat')-1);
            
            % find directory containing this filename
            cd ..
            cmd=['!find . -maxdepth 2 -name ' fnameNotMat ' > temptemp.txt'];
            eval(cmd)
            
            fid=fopen('temptemp.txt');
            filedir=fgetl(fid);
            while filedir(9)~='_';
                % then this dir is not formatted as I normalyl do. keep
                % looking
                filedir=fgetl(fid);
                if ~ischar(filedir)
                    % then give up. take first directory, reopen the file.
                    fid2=fopen('temptemp.txt');
                    filedir=fgetl(fid2);
                    break
                end
            end
            
            slashes=findstr(filedir, '/');
            filedir=filedir(slashes(1)+1:slashes(2)-1); % directory holding song files
            disp(filedir);
            
            
            % extract timing of labeled notes
            cd(filedir);
            NotMatDat=load(fnameNotMat);
            if ~strcmp(NotMatDat.labels, LabelsAll{ind})
                disp([birdname '-' exptname '-' fnameNotMat '; ERROR: saved labels not same as just loaded labels']); % confirm that labels are as from before
            end
            
            % read cbin file
            [dat, Fs, DOFILT, ext]=ReadDataFile(fnameCbin,'0'); 
            
            % get segment of interest
            TimeOnset=NotMatDat.onsets(NoteIndRange(1)); % ms
            TimeOnset=TimeOnset-1000*PreDur;           
            TimeOffset=NotMatDat.offsets(NoteIndRange(end));
            TimeOffset=TimeOffset+1000*PostDur;
            
            TimeOnset=ceil(Fs*TimeOnset/1000); % convert to samples
            TimeOffset=ceil(Fs*TimeOffset/1000); % convert to samples
            
            dattmp=dat(TimeOnset:TimeOffset);
            
            [sm,sp,t,f]=SmoothData(dattmp,Fs,DOFILT,'‘hanningfirff’'); % ends up doing buttor, with filtfilt. (because there is extra quote marks in filter type)
            sp=sp(10:140,:);
            f=f(10:140);
            
            % -- plot
            sptemp=sp;
            pp=find(sptemp>0);
            mntmp = min(min(sptemp(pp)));
            pp=find(sptemp==0);
            sptemp(pp) = mntmp;
            sptemp=log(sptemp);
            
            % subtract min from all
            sptemp=sptemp-min(min(sptemp));
            
            % -- determine clim and cmax (thresholds)
            sptemp_vector=reshape(sptemp, numel(sptemp), 1);
            Xcenters=linspace(min(sptemp_vector), max(sptemp_vector), length(sptemp_vector)/200);
            [Ybinned, ~, ~]=lt_plot_histogram(sptemp_vector, Xcenters, 0, '','',1,'k');
            plot(Xcenters, Ybinned);
            
            [~, maxind]=max(Ybinned); % first peak in distrubtion of magnitudes
            Cmin=Xcenters(maxind);
            Cmax=max(sptemp_vector);
            
            % Plot
            colormap('hot')
            imagesc(t, f, sptemp, [Cmin+Cmin/10 Cmax]);
            axis([t(1) t(end) f(1) f(end)]);
            

            
            % labels
            xlabel('time (s)');
            ylabel('frequency (hz)');
            end
        end

          linkaxes(hsplotsAll, 'xy');
lt_subtitle([birdname '-' exptname]);
lt_plot_annotation(1, 'note: each plot zaxis norm to itself', 'r');
    
    end
end
end




%% ====== FOR EACH EXPT PLOT EXAMPLES OF ALL SYL TYPES (EVERYTHING IN SYLS UNIQUE) [sort by type and sequence]
% figcount=1;
% subplotrows=6;
% subplotcols=3;
% fignums_alreadyused=[];
% hfigs=[];



for i=1:NumBirds;
    
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpts;
        
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        
        if ~strcmp(birdname, BirdToPlot) | ~strcmp(exptname,ExptToPlot);
            continue;
        end
        
        lt_figure; hold on;
        hsplotsAll=[];
        Climits=[];

        
        % ==== GO TO FOLDER WITH SONGS
        WNday1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
        FirstDay_date=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.FirstDay;
        FirstDay_datenum=datenum(FirstDay_date, 'ddmmmyyyy');
        LastBaselineDay_datenum=FirstDay_datenum+WNday1-2;
        LastBaselineDay_date=datestr(LastBaselineDay_datenum, 'ddmmmyyyy');
        
        cd(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.PlotLearningDir);
        cd ..
        
        try
        load(['RawDatStruct_' LastBaselineDay_date '.mat']);
        ParamsDay=load(['Params_' LastBaselineDay_date '.mat']);
        catch err
            % go back one day, day has no data
            LastBaselineDay_date=datestr(LastBaselineDay_datenum-1, 'ddmmmyyyy');
                    load(['RawDatStruct_' LastBaselineDay_date '.mat']);
        ParamsDay=load(['Params_' LastBaselineDay_date '.mat']);

        end
            
        % === collect labels
        sylfields=fieldnames(RawDatStruct.data);
        
        if any(strcmp(sylfields, 'b'));
            SylToUse='b';
        else
            SylToUse=sylfields{1};
        end
        LabelsAll=RawDatStruct.data.(SylToUse)(:, 7);
        
        
        % ==== plot multiple renditions of each syl ======================
        
        SylListsToPlot={'targsyl', 'SylFields_Unique_STSS', 'SylFields_Unique_STDS', 'SylFields_Unique_DTSS', 'SylFields_Unique_DTDS'};
        
        % figure out num columns of subplots based on number of syls
        numcols=6;
        for k=1:length(SylListsToPlot);
            numcols=max([numcols length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.(SylListsToPlot{k}))]);
        end
        
        for k=1:length(SylListsToPlot);
            SylList=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.(SylListsToPlot{k});
            
            if strcmp(SylListsToPlot{k}, 'targsyl');
                SylList={SylList};
            end
            
            for j=1:length(SylList);
                syl=SylList{j};
                
%                 numcols=6;
%                 if length(SylList)>6;
%                     numcols=length(SylList);
%                 end
                
                hsplot=lt_subplot(length(SylListsToPlot), numcols, numcols*(k-1)+j); hold on;
                if strcmp(SylListsToPlot{k}, 'targsyl');
                     title([syl '-targ']);
                else
                    
                    title([syl '-' SylListsToPlot{k}(end-3:end)]);
                end
                %                     [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                %                     title(syl);
                hsplotsAll=[hsplotsAll hsplot];
                
                NoteInd=[];
                ccc=1;
                while isempty(NoteInd);
                    ind=randi(length(LabelsAll), 1); % pick a random song
                    NoteInd=regexp(LabelsAll{ind}, lower(syl)); % get matches for this motif
                    ccc=ccc+1;
                    if ccc==20;
                        disp(['PROBLEM - CANT FIND SONG WITH THIS MOTIF ' syl]);
                        break
                    end
                end
                
                NoteInd=NoteInd(randi(length(NoteInd),1)); % pick a random note
                % get the exact note location (in the inputed string)
                if length(syl)>1;
                    CapsLocation=regexp(syl, '[A-Z]');
                    assert(~isempty(CapsLocation), 'PROBLEM - length(syl)>1,but cant find caps...?!');
                    NoteInd=NoteInd+CapsLocation-1;
                end
                
                
                NoteIndRange=NoteInd-NumPreNotes:NoteInd+NumPostNotes;  % range of note inds for entire motif
                % make sure range is not out of song bounds
                if NoteIndRange(1)<1;
                    NoteIndRange=NoteIndRange(NoteIndRange>0);
                elseif NoteIndRange(end)>length(LabelsAll{ind});
                    NoteIndRange=NoteIndRange(NoteIndRange<=length(LabelsAll{ind}));
                end
                
                
                % === extract sound data for that note range and plot
                fnameNotMat=RawDatStruct.data.(SylToUse){ind, 5};
                fnameCbin=fnameNotMat(1:strfind(fnameNotMat, '.not.mat')-1);
                
                % find directory containing this filename
                cd ..
                cmd=['!find . -maxdepth 2 -name ' fnameNotMat ' > temptemp.txt'];
                eval(cmd)
                
             fid=fopen('temptemp.txt');
            filedir=fgetl(fid);
            while filedir(9)~='_';
                % then this dir is not formatted as I normalyl do. keep
                % looking
                filedir=fgetl(fid);
                if ~ischar(filedir)
                    % then give up. take first directory, reopen the file.
                    fid2=fopen('temptemp.txt');
                    filedir=fgetl(fid2);
                    break
                end
            end
               
                slashes=findstr(filedir, '/');
                filedir=filedir(slashes(1)+1:slashes(2)-1); % directory holding song files
                disp(filedir);
                
                % extract timing of labeled notes
                cd(filedir);
                NotMatDat=load(fnameNotMat);
            if ~strcmp(NotMatDat.labels, LabelsAll{ind})
                disp([birdname '-' exptname '-' fnameNotMat '; ERROR: saved labels not same as just loaded labels']); % confirm that labels are as from before
            end
                
                % read cbin file
                [dat, Fs, DOFILT, ext]=ReadDataFile(fnameCbin,'0');
                
                % get segment of interest
                TimeOnset=NotMatDat.onsets(NoteIndRange(1)); % ms
                TimeOnset=TimeOnset-1000*PreDur;
                TimeOffset=NotMatDat.offsets(NoteIndRange(end));
                TimeOffset=TimeOffset+1000*PostDur;
                
                TimeOnset=ceil(Fs*TimeOnset/1000); % convert to samples
                TimeOffset=ceil(Fs*TimeOffset/1000); % convert to samples
                
                dattmp=dat(TimeOnset:TimeOffset);
                
                [sm,sp,t,f]=SmoothData(dattmp,Fs,DOFILT,'‘hanningfirff’'); % ends up doing buttor, with filtfilt. (because there is extra quote marks in filter type)
                sp=sp(10:140,:);
                f=f(10:140);
                
                
                % -- plot
                sptemp=sp;
                pp=find(sptemp>0);
                mntmp = min(min(sptemp(pp)));
                pp=find(sptemp==0);
                sptemp(pp) = mntmp;
                sptemp=log(sptemp);
                
                
                % Plot
                if isempty(Climits)
                    % -- determine clim and cmax (thresholds)
                    sptemp_vector=reshape(sptemp, numel(sptemp), 1);
                    Xcenters=linspace(min(sptemp_vector), max(sptemp_vector), length(sptemp_vector)/200);
                    [Ybinned, ~, ~]=lt_plot_histogram(sptemp_vector, Xcenters, 0, '','',1,'k');
                    plot(Xcenters, Ybinned);
                    
                    [~, maxind]=max(Ybinned); % first peak in distrubtion of magnitudes
                    Cmin=Xcenters(maxind);
                    Cmax=max(sptemp_vector);
                    Climits=[Cmin+Cmin/10 Cmax];
                end
                colormap('hot')
                imagesc(t, f, sptemp, Climits);
                axis([t(1) t(end) f(1) f(end)]);
                
                % plot line demarcating syl
                NoteTimeOnset=NotMatDat.onsets(NoteInd)-20; % 16 is emprical (due to smoothing)
                NoteTimeOnset_RelSylOnset=NoteTimeOnset/1000-(TimeOnset/Fs); % sec
                NoteTimeOffset=NotMatDat.offsets(NoteInd)+10;
                NoteTimeOffset_RelSylOnset=NoteTimeOffset/1000-(TimeOnset/Fs); % sec
                
                line([NoteTimeOnset_RelSylOnset NoteTimeOnset_RelSylOnset], ylim, 'Color','m', 'LineStyle','--');
                line([NoteTimeOffset_RelSylOnset NoteTimeOffset_RelSylOnset], ylim, 'Color','m', 'LineStyle','--');
                
                % labels
                if j==1;
                    ylabel('frequency (hz)');
                    if k==length(SylListsToPlot);
                        xlabel('time (s)');
                    end
                end
                
                % === annotate with acoustic dist and such if desired
                Xlim=xlim;
                Ylim=ylim;
                
                if PlotGeneralization==1;
                generali=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).LEARNING.learning_metric.mean_rel_targ;
                lt_plot_text(Xlim(1), Ylim(end)-0.1*(Ylim(2)-Ylim(1)), ['Gen: ' num2str(generali, '%3.2g')], 'w');
                end
                if PlotAcousticDist==1;
                acousti=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).eucldist_from_targ_PCA;
                lt_plot_text(Xlim(1), Ylim(end)-0.2*(Ylim(2)-Ylim(1)), ['AcDist: ' num2str(acousti, '%3.2g')], 'w');
                end
                
            end
        end
        
        

        
        linkaxes(hsplotsAll, 'xy');
        lt_subtitle([birdname '-' exptname]);
        lt_plot_annotation(1, 'note: within expt all using same z scale (based on first plot)', 'r');
        
        
    end
end





