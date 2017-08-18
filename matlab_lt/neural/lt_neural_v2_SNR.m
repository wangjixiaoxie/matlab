function [SNR, SignalPower, NoisePower] = lt_neural_v2_SNR(FRmat, plotOn)

% FRmat: timebins (assumes ms resolution) x trials 
% plotOn : then plots on a new figure;

%% LT 7/9/17 -
% takes smoothed FR traces, assumed to be aligned, and extracts signal and
% noise power as in Sahani & Linden, NIPS, 2002

%% calculate


N = size(FRmat, 2);

% 1st take average, then take power
Rmean = mean(FRmat, 2); % mean response
PowerOfMean = sum((Rmean - mean(Rmean)).^2);


% 1st take power, then take average
tmp = mean(FRmat,1);
tmp = repmat(tmp, size(FRmat,1), 1);
Power_each_trial = sum((FRmat - tmp).^2, 1);
MeanPowerOverTrials = mean(Power_each_trial);

% === get estimated power of signal
SignalPower = (1/(N-1)) * (N*PowerOfMean - MeanPowerOverTrials);
NoisePower = MeanPowerOverTrials - SignalPower;

SNR = SignalPower/NoisePower;


%% === plot trials and mean.
if plotOn ==1
    
lt_figure; hold on;
plot(1:size(FRmat, 1), FRmat, '-');

plot(1:size(FRmat,1), Rmean, '-k', 'LineWidth', 3);

lt_plot_annotation(1, ['N = ' num2str(N) ', S=' num2str(SignalPower, '%3.2g') ', N=' ....
    num2str(NoisePower, '%3.2g') ', SNR='  num2str(SNR, '%3.2g')], 'r')
end




