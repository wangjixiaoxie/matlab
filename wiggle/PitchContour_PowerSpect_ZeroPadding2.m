function PitchContour_PowerSpect_ZeroPadding2(note)

% This is siimlar to PitchContour_PowerSpect.m but calculate power spectrum
% with zero padding to have greater frequency resolution.
% Also, tapler signal edges with 50 and 100 Hz cosine functions, instead of
% 10 & 20 Hz in the previous script.

% This is siimlar to PitchContour_PowerSpect_ZeroPadding.m but using input
% signal with mean subtracted to remove DC component.

NewNsample = 2^12; % Nunmber of samples after zero-padding.

load(['PitchContour_',note,'_AutoCorr.mat'])
%open(['PitchContour_Overlay_',note,'_Norm.fig'])

data = AutoCorrAnal.Contour_Deviation;
time = AutoCorrAnal.Contour_Time;
TimeOffset = time(1);
Ndata = length(data(:,1));
Nsample = length(data(1,:));
Fs = 1/((time(end)-time(1))/Nsample);
ContourLength = time(end)-time(1);


% Make taper function at the ends of signal
% 100Hz sine function
taperFs = 100;
Ntaper = floor(1/taperFs/2/(1/Fs));
taperT= 0:1/Fs:(1/Fs)*Ntaper;
taper1 = fliplr(cos(2*pi*taperFs*taperT));
taper2 = cos(2*pi*taperFs*taperT);
if Nsample>2*length(taperT)
    TP100 = 1;
    taper100Hz = [taper1 ones(1,Nsample-(Ntaper*2)-2) taper2];
    taper100Hz = taper100Hz/2+0.5;
else
    TP100=0;
end

% 200Hz sine function
taperFs=200;
Ntaper = floor(1/taperFs/2/(1/Fs));
taperT= 0:1/Fs:(1/Fs)*Ntaper;
taper1 = fliplr(cos(2*pi*taperFs*taperT));
taper2 = cos(2*pi*taperFs*taperT);
if Nsample>2*length(taperT)
    TP200 = 1;
    taper200Hz = [taper1 ones(1,Nsample-(Ntaper*2)-2) taper2];
    taper200Hz = taper200Hz/2+0.5;
else
    TP200=0;
end


figure
for n=1:Ndata
    signal = data(n,:);
    
    subplot(3,1,1)
    RawCont(n,:) = signal-mean(signal);
    %RawCont(n,:) = signal; % Do not remove DC component
    plot(time,RawCont(n,:)); hold on
    title('Raw data')
    xlabel('Time (sec)');ylabel('Pitch (%)')
    
    subplot(3,1,2)
    if TP100>0
    TaperCont100Hz(n,:) = RawCont(n,:).*taper100Hz;
    plot(time,TaperCont100Hz(n,:)); hold on
    title('Tapered data (100Hz cosine function)')
    xlabel('Time (sec)');ylabel('Pitch (%)')
    end
    
    subplot(3,1,3)
    if TP200>0
    TaperCont200Hz(n,:) = RawCont(n,:).*taper200Hz;
    plot(time,TaperCont200Hz(n,:)); hold on
    title('Tapered data (200Hz cosine function)')
    xlabel('Time (sec)'); ylabel('Pitch (%)')
    end
    
end
sp{1} = ['Pitch contour of note "',note,'"']; 
sp{2} = ['(% change from x-trial mean, mean subtracted in each trial)'];
suptitle(sp)
saveas(gcf,['PitchContour_PowerSpect1_',note,'_MeanSubtr.fig'])


%%%%%% Compute power spectrum %%%%%%
% There seem to be two methods to calculate power specral density
% (1) http://www.mathworks.com/help/signal/ug/psd-estimate-using-fft.html
% (2) http://www.mathworks.com/help/matlab/examples/fft-for-spectral-analysis.html
% I use the 1st one in the following script but the results are similar btw
% the two mothods although acrutal values are different...
% Goldberg & Fee used the 1st method too.

figure

% Raw data

f = 0:Fs/NewNsample:Fs/2;
for n=1:Ndata    
    Y = fft(RawCont(n,:),NewNsample);
    Y = Y(1:length(f));
    Y2 = (1/(Fs*NewNsample)).*abs(Y).^2;
    %Y2(2:end-1) = 2*Y2(2:end-1);
    Power(n,:)= Y2;
end

for t=1:length(Power(1,:))
    MeanPower(t) = mean(Power(:,t));
    MedianPower(t) = median(Power(:,t));
    PowerSD(t) = std(Power(:,t));
end


subplot(3,1,1)
for n=1:Ndata
    plot(f,Power(n,:)); hold on
end
plot(f,MeanPower,'r'); hold on
%plot(f,MeanPower+PowerSD,':r'); hold on
%plot(f,MeanPower-PowerSD,':r'); hold on
plot(f,MedianPower,'g'); hold on
xlim([-1 500])
ylabel('Power'); % xlabel('Frequency(Hz)')
title('Raw data')

RawData.PowerAll = Power;
RawData.PowerMean = MeanPower;
RawData.PowerMedian = MedianPower;
RawData.PowerSD = PowerSD;
RawData.Contours = RawCont;


% Taper 100Hz data
if TP100>0
    
    f = 0:Fs/NewNsample:Fs/2;
    for n=1:Ndata
        Sig = TaperCont100Hz(n,:)-mean(TaperCont100Hz(n,:));
        Y = fft(Sig,NewNsample);
        Y = Y(1:length(f));
        Y2 = (1/(Fs*NewNsample)).*abs(Y).^2;
        %Y2(2:end-1) = 2*Y2(2:end-1);
        Power(n,:)= Y2;
    end

    MeanPower = mean(Power);
    MedianPower = median(Power);
    PowerSD = std(Power);

    subplot(3,1,2)
    for n=1:Ndata
        plot(f,Power(n,:)); hold on
    end
    plot(f,MeanPower,'r'); hold on
    %plot(f,MeanPower+PowerSD,':r'); hold on
    %plot(f,MeanPower-PowerSD,':r'); hold on
    plot(f,MedianPower,'g'); hold on
    xlim([-1 500])
    ylabel('Power'); % xlabel('Frequency(Hz)')
    title('Tapered data (100Hz)')

    Taper100HzData.PowerAll = Power;
    Taper100HzData.PowerMean = MeanPower;
    Taper100HzData.PowerMedian = MedianPower;
    Taper100HzData.PowerSD = PowerSD;
    Taper100HzData.Contours = TaperCont100Hz;

else

    Taper100HzData.PowerAll = [];
    Taper100HzData.PowerMean = [];
    Taper100HzData.PowerMedian = [];
    Taper100HzData.PowerSD = [];
    Taper100HzData.Contours = [];
    
end
    

% Taper 200Hz data
if TP200>0
    
    f = 0:Fs/NewNsample:Fs/2;
    for n=1:Ndata    
        Sig = TaperCont200Hz(n,:)-mean(TaperCont200Hz(n,:));
        Y = fft(Sig,NewNsample);
        Y = Y(1:length(f));
        Y2 = (1/(Fs*NewNsample)).*abs(Y).^2;
        %Y2(2:end-1) = 2*Y2(2:end-1);
        Power(n,:)= Y2;
    end

    MeanPower = mean(Power);
    MedianPower = median(Power);
    PowerSD = std(Power);

    subplot(3,1,3)
    for n=1:Ndata
        plot(f,Power(n,:)); hold on
    end
    plot(f,MeanPower,'r'); hold on
    %plot(f,MeanPower+PowerSD,':r'); hold on
    %plot(f,MeanPower-PowerSD,':r'); hold on
    plot(f,MedianPower,'g'); hold on
    xlim([-1 500])
    ylabel('Power'); % xlabel('Frequency(Hz)')
    title('Tapered data (200Hz)')

    Taper200HzData.PowerAll = Power;
    Taper200HzData.PowerMean = MeanPower;
    Taper200HzData.PowerMedian = MedianPower;
    Taper200HzData.PowerSD = PowerSD;
    Taper200HzData.Contours = TaperCont200Hz;
    
else
    
    Taper200HzData.PowerAll = [];
    Taper200HzData.PowerMean = [];
    Taper200HzData.PowerMedian = [];
    Taper200HzData.PowerSD = [];
    Taper200HzData.Contours = [];
    
end

suptitle(['Power spectral density of note "',note,'"'])
saveas(gcf,['PitchContour_PowerSpect2_',note,'_MeanSubtr.fig'])


FreqAxis = f;
ContourTime = time;

save(['PitchContour_PowerSpect_',note,'_MeanSubtr.mat'],'RawData','Taper100HzData','Taper200HzData',...
    'FreqAxis','ContourLength','ContourTime','note')




