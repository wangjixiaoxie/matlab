function [max_power, best_per]=cepstrum(spectrum);

plot_flag=1;

%allowable frequencies for fundamental, in bins of input spectrum
%for spctrogram with 16 ms window bin width = 62 ms
f_min=6;
f_max=18;


%length of spect
lspect=size(spectrum,1);
nbins=size(spectrum,2);

%calculate psd of each spectrum: the values used for zero padding and window size
%  will have affect on interopolation and smoothing;

cepmatrix=[];
for i=1:nbins
  [cep,f]=psd(spectrum(:,i),256,1,length(spectrum(:,i)));      
  cepmatrix=[cepmatrix,cep];
end


%find maximum in acceptable range
periods=1./f;  %periods in bins
acc_ind=find(f_min <= periods & periods <= f_max);

  
cut_cepmatrix=cepmatrix(acc_ind,:);
cut_periods=periods(acc_ind,:);
[max_power,ind]=max(cut_cepmatrix);
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
  plot(f,cepmatrix);
  axis([0,.5,0,max(cepmatrix)]);
  subplot(2,1,2)
  plot(spectrum,'m');
  hold on
    for i=1:length(fit)
     line([fit(i),fit(i)],[0,max(spectrum)])
   end
  hold off
end



