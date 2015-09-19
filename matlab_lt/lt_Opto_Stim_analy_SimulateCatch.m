function [DatStructCompiled, Params]=lt_Opto_Stim_analy_SimulateCatch(DatStructCompiled,Params);
%% LT 1/18/15 - SIMULATE CATCH BY TAKING RANDOM SAMPLES
% Is useful if using catch songs - there all trials are called catch (in .rec file), even
% if you have catch trials as well. Therefore, if want to get stim catch
% and notcatch (trial-wise) take subsample.

NumTrials=length(DatStructCompiled.All);

%% Subsample

SimCatchInd=[];
SimNotCatchInd=[];

for i=1:NumTrials;
if rand<Params.SimCatchFrac;
 SimCatchInd=[SimCatchInd i];
else 
 SimNotCatchInd=[SimNotCatchInd i];
end
end

%% Divide all trials into SimCatch or SimNotCatch

DatStructCompiled.SIMULATE_Catch=DatStructCompiled.All(SimCatchInd);
DatStructCompiled.SIMULATE_NotCatch=DatStructCompiled.All(SimNotCatchInd);

%% Save

cd(Params.savefolder);

save('DatStructCompiled_SIMULATEDCatch.mat','DatStructCompiled');
save('Params_SIMULATEDCatch.mat','Params');




