function [max_power, best_per]=costrans(spectrum);

win_flag=0;
plot_flag=1;

%allowable frequencies for fundamental, in bins of input spectrum
%for spctrogram with 16 ms window bin width = 62 ms
f_min=5;
f_max=20;


%length of spect
lspect=size(spectrum,1);
nbins=size(spectrum,2);


if win_flag==1
 %window spectrum
 window=hanning(lspect);
 win_mat=window*ones(1,size(spectrum,2));
 spectrum=spectrum.*win_mat;
end

%calculate dct of each spectrum: the values used for zero padding and window size
%  will have affect on interopolation and smoothing;


cost_matrix=dct(spectrum);
%cost_matrix=cost_matrix.^2;

%find maximum in acceptable range
coef_freq=[0:lspect-1]./2;  %cycles per length of spectrum
periods=lspect./coef_freq;  %periods in bins
acc_ind=find(f_min <= periods & periods <= f_max);
periods=periods';

  
cut_cost_matrix=cost_matrix(acc_ind,:);
cut_periods=periods(acc_ind,:);
[max_power,ind]=max(cut_cost_matrix);
best_per=cut_periods(ind);


%dissallow values at extremes of range (require real peaks);

%min_acc=min(acc_ind);
%max_acc=max(acc_ind);

%for i= 1:size(cepmatrix,2)
  %for each column, get all maxima
%  [maxima,g,max_ind]=findmax1(cepmatrix(:,i))
  %lookfor maxima in acceptable range
%  [max_ind_acc]=find(min_acc <= max_ind <= max_acc); 
  


if plot_flag==1
  fit=[best_per:best_per:length(spectrum)];
  figure
  subplot(2,1,1);
  coef_freq=coef_freq./lspect;
  plot(coef_freq,cost_matrix);
  axis([0,.5,0,max(cost_matrix)]);
  subplot(2,1,2)
  plot(spectrum,'m');
  hold on
    for i=1:length(fit)
     line([fit(i),fit(i)],[0,max(spectrum)])
   end
  hold off
end



