function PitchContour_AutoCorr_4MimisData

%if nargin<1
    BatchFile = input('Enter the name of batch file name (e.g. batch.undir.20060802): ','s');
    NoteType = input('Enter note type that you want to analyze: ','s')
%end

Q2 = input('Do you want to calculate pitch contours? (y/n): ','s');
SavedDataName = ['PitchContour_',NoteType,'.mat'];
if Q2(1) == 'y'

    Batch = readtextfile(BatchFile)

    clear PitchContour; K=1;
    for n=1:length(Batch(:,1))
        wavfilename = deblank(Batch(n,:));
        
        if ~isempty(wavfilename)
        
        disp(['analyzing ',wavfilename,'...'])
%         labelfile = [wavfilename(1:end-5),'.not.mat'];
        labelfile = [wavfilename(1:end-5),'.cbin.not.mat'];
        NotMatFile = load(labelfile);
        labels = NotMatFile.labels;
        onsets = NotMatFile.onsets;
        offsets = NotMatFile.offsets;
        NotePosition = find(labels==NoteType);

        if length(NotePosition)>0
            NoteOnset = onsets(NotePosition);
            NoteOffset = offsets(NotePosition);

            % load .wav file and fileter it.
            %[wavfile,Freq,Nbits] = wavread(wavfilename);
            [wavfile, Freq] = read_filt(wavfilename);
            flow = 250;
            nfir = 512;
            ndelay = fix(nfir/2) ;
            bfir = fir1(nfir, flow*2/Freq, 'high');
            filt_wavfile = filtfilt(bfir,1,wavfile);

             % segment the wav file, normalize its amplitude
            for m=1:length(NotePosition)
                Note = filt_wavfile(round(NoteOnset(m)/1000*Freq) :...
                    round(NoteOffset(m)/1000*Freq));
                Note = Note./max(abs(Note)).*0.9;
                Silence = zeros(round(Freq./1000*5),1); 
                Note = [Silence;Note;Silence]; % add 5 ms of silent periods before and after the note
                FileName = [wavfilename,'.',NoteType,num2str(m),'.PichContour.fig'];
                %wavwrite(Note,Freq,Nbits,NoteFileName);


                if K==1;
                    figure
                    specgram(Note,1024,44100,512,500); hold on
                    colormap(flipud(gray)); 
                    caxis([-30 30])
                    ylim([0 3000])
                    xlabel('time (sec)'); ylabel('frequency (Hz)');

                    F_low = input('Enter lower limit of fundamental frequency:');
                    F_high = input('Enter higher limit of fundamental frequency:');
                    close
                end

                % %Default values
                N = 1024;
                OVERLAP = 1020;
                sigma = 2;
                harms = 1;
                filetype = 'w';
                returnAV = 0;
                PlotSpec = 1;

                [pitch_data,TimeAxis] =...
                    jc_pitchmat1024sk(Note,N,OVERLAP,sigma,F_low,F_high,harms,filetype,returnAV,PlotSpec);
                title(FileName)
                saveas(gcf,FileName)
                disp([FileName,' was saved!'])
                %close all

                PitchContour.FileName{K} = FileName;
                PitchContour.PitchData{K} = pitch_data;
                PitchContour.TimeAxis{K} = TimeAxis;
                K=K+1;
                if K>20
                    close all
                end
            end
        end
        end
    end
    PitchContour.F_low = F_low;
    PitchContour.F_high = F_high;
    
    save(SavedDataName, 'PitchContour')    
    close all

end

load(SavedDataName)

% plot all notes
figure('Position',[0 1200 700 900]);
subplot(2,1,1)
for n=1:length(PitchContour.PitchData)
        plot(1000.*PitchContour.TimeAxis{n}, PitchContour.PitchData{n},'r-'); hold on
        text(1000*PitchContour.TimeAxis{n}(end),PitchContour.PitchData{n}(end),num2str(n),'FontSize',12)
        %text(1000*PitchContour.TimeAxis{n}(1)-2,PitchContour.PitchData{n}(1),num2str(n),'FontSize',12)
end
ylim([PitchContour.F_low,PitchContour.F_high])
xlabel('time (msec)'); ylabel('frequency (Hz)');
title(['Pitch contours of syllable "',NoteType,'"'])
saveas(gcf,['PitchContour_Overlay_',NoteType,'.fig'])


Q = input('Do you want to do auto-corerlation analysis? (y/n): ','s');
if Q(1)=='n'
    return
end


% Specify the time window to analyze
OnsetTime = input('Enter the onset fime of the segment that you measure (msec): ');
OffsetTime = input('Enter the offset fime of the segment that you measure (msec): ');

PitchData = PitchContour.PitchData;
TimeAxis = PitchContour.TimeAxis;

m=1;
for n=1:length(TimeAxis)
    Time = TimeAxis{n};
    Contour = PitchData{n};
    if Time(1)<=OnsetTime/1000 && Time(end)>=OffsetTime/1000
        G = find(Time>=OnsetTime/1000 & Time<=OffsetTime/1000);
        ContourSeg(m,:) = Contour(G);
        TimeSeg = Time(G);
        m=m+1;
    end
end


subplot(2,1,2)
for n=1:length(ContourSeg(:,1))
    plot(1000.*TimeSeg,ContourSeg(n,:),'r'); hold on
end    
xlabel('time (msec)'); ylabel('frequency (Hz)');
ylim([PitchContour.F_low PitchContour.F_high])
title('Segment to analyze auto-correlation')


% % Discard outliers.
% DiscardData = input('Enter the number of traces to be discarded (e.g. 2,5,10): ','s')
% if length(DiscardData)<1; DiscardData2 = [];
% else
%     com = find(DiscardData==',');
%     if length(com)<1
%         DiscardData2 = str2num(DiscardData);
%     else
%         com2 = [0,com];
%         for n=1:length(com2)
%             if n<length(com2)
%                 DiscardData2(n) = str2num(DiscardData(com2(n)+1:com2(n+1)-1))
%             else
%                 DiscardData2(n) = str2num(DiscardData(com2(n)+1:end))
%             end
%         end
%     end
% end

% Get threshold for excluding outliers
Q3 = input('Do you want to discard outliers? (y/n) ','s');
if Q3(1) == 'y';
    Ex_low = input('Enter the lower limit of frequency to exclude outliers: ');
    Ex_hi = input('Enter the higher limit of frequency to exclude outliers: ');
    plot([OnsetTime OffsetTime],[Ex_low Ex_low],'g--'); hold on
    plot([OnsetTime OffsetTime],[Ex_hi Ex_hi],'g--'); hold on

    % Discard outliers
    DiscardData = []; i=1;
    for n=1:length(ContourSeg(:,1))
        if min(ContourSeg(n,:))<Ex_low || max(ContourSeg(n,:))>Ex_hi
            DiscardData(i) = n;
            i=i+1;
        end
    end   
    ContourSeg(DiscardData,:) = [];
    for n=1:length(ContourSeg(:,1))
        plot(1000.*TimeSeg,ContourSeg(n,:),'b'); hold on
    end 
else
    Ex_low = []; Ex_hi = []; DiscardData = [];
end

saveas(gcf,['PitchContour_Overlay_',NoteType,'.fig'])


figure
subplot(2,1,1)
ContourSegMean = mean(ContourSeg);
for n=1:length(ContourSeg(:,1))
    plot(1000.*TimeSeg,ContourSeg(n,:)); hold on
    plot(1000.*TimeSeg,ContourSegMean,'r'); hold on
end
title('Raw pitch contours and their mean')
ylabel('Frequency')

% Calculate % change from mean and plot them
subplot(2,1,2)
for n=1:length(ContourSeg(:,1))
    ContourSegChange(n,:) = (ContourSeg(n,:)-ContourSegMean)./ContourSegMean.*100;
    plot(1000.*TimeSeg,ContourSegChange(n,:)); hold on
end
plot([min(1000.*TimeSeg) max(1000.*TimeSeg)],[0 0],'r'); hold on
title('% different from mean')
ylabel('Percent'); xlabel('Time(msec)')
saveas(gcf,['PitchContour_Overlay_',NoteType,'_Norm.fig'])

% Calculate xcov
figure
subplot(2,1,1)
xc = []; R2 = []; R2_time=[];
for n=1:length(ContourSegChange(:,1))
    Seg = ContourSegChange(n,:);
    xc(n,:) = xcorr(Seg,'coeff');
    R2(n,:) = xc(n,:).^2;
    R2_time = ((1:length(R2(n,:)))-length(Seg)).*(Time(2)-Time(1)).*1000; 
    plot(R2_time,R2(n,:),'b:'); hold on
end
R2_mean = mean(R2);
plot(R2_time,R2_mean,'r'); hold on
xlim([-50 50])
xlabel('Time lag (ms)')
ylabel('R^2')
title('Coefficient of determination (R^2) of pitch contours')

R2_mean_minus = fliplr(R2_mean(1:length(Seg)));
R2_mean_plus = R2_mean(length(Seg):end);
R2_mean_abs = (R2_mean_minus+R2_mean_plus)./2;
R2_mean_abs_time = R2_time(R2_time>=0);

% Measure width
% [a,b] = min(abs(R2_mean_abs-0.5)); % half width
% [c,d] = min(abs(R2_mean_abs-(1/exp(1)))); % time constant
% HalfWidth = b*(Time(2)-Time(1)).*1000;
% TimeConstant = d*(Time(2)-Time(1)).*1000;
HalfWidth = interp1(R2_mean_abs,R2_mean_abs_time,0.5);
TimeConstant = interp1(R2_mean_abs,R2_mean_abs_time,1/exp(1));

% plot results
subplot(2,1,2)
plot(R2_mean_abs_time,R2_mean_abs,'r'); hold on
plot([0 HalfWidth],[.5 .5],'g');hold on
plot([HalfWidth HalfWidth],[0 .5],'g');hold on
plot([0 TimeConstant],[(1/exp(1)) (1/exp(1))],'m');hold on
plot([TimeConstant TimeConstant],[0 (1/exp(1))],'m');hold on
text(HalfWidth-.2,-.03,num2str(HalfWidth),'Color','g'); hold on
text(TimeConstant-.2,-.06,num2str(TimeConstant),'Color','m'); hold on
xlabel('Time lag (ms)')
ylabel('R^2')
title('Half width (green) and time constant (magenta)')
saveas(gcf,['PitchContour_Overlay_',NoteType,'_AutoCorr.fig'])

% Combine variables and save them
AutoCorrAnal = [];
AutoCorrAnal.DiscardThreshold_low = Ex_low;
AutoCorrAnal.DiscardThreshold_hi = Ex_hi;
AutoCorrAnal.DiscardData = DiscardData;
AutoCorrAnal.OnsetTime = OnsetTime ;
AutoCorrAnal.OffsetTime = OffsetTime;
AutoCorrAnal.Contour_All = ContourSeg;
AutoCorrAnal.Contour_Mean = ContourSegMean;
AutoCorrAnal.Contour_Deviation = ContourSegChange;
AutoCorrAnal.Contour_Time = TimeSeg;
AutoCorrAnal.R2_All = R2;
AutoCorrAnal.R2_Time = R2_time;
AutoCorrAnal.R2_Mean = R2_mean;
AutoCorrAnal.R2_MeanAbs = R2_mean_abs;
AutoCorrAnal.R2_MeanAbs_Time = R2_mean_abs_time;
AutoCorrAnal.R2_HalfWidth = HalfWidth;
AutoCorrAnal.R2_TimeConstant = TimeConstant;

SavedDataName2 = ['PitchContour_',NoteType,'_AutoCorr.mat'];
save(SavedDataName2, 'AutoCorrAnal')  




    