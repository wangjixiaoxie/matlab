function [best_freq, best_idx, ratio_prods, ratio_freq_scale]=deltaprods(spectrum,f_step);

%Note: function deltaprods assumes that first bin of spectrum represents a frequency of f_step, NOT 0.
% ratio_prods has 2 columns, the first contains diff_prods
%  the second contains percent_harm,  see below for how these are calculated

plot_flag=0; %set to 1 or 2 for debugging & display

%set range of ff over which to search
ff_min=300;
ff_max=1200;

%calculate min & max periods in bins

min_per=round(ff_min/f_step);
step_per=.05;
max_per=round(ff_max/f_step);

%length of spect
lspect=size(spectrum,1);
base=zeros(1,lspect);
cos_matrix=[];
sin_matrix=[];
n=0;
% factors of appropriate periodicities
indexes=[0:lspect-1];
for i=min_per:step_per:max_per
  n=n+1;
  hidx=round([i:i:lspect-.5*i]);
  nhidx=round([i+.5*i:i:lspect]);
  if length(hidx) > length(nhidx)
     hidx=hidx(1:length(nhidx));
  elseif length(nhidx) > length(hidx)
     nhidx=nhidx(1:length(hidx));
  end
  factor1=base;
  factor1(hidx)=ones(size(hidx))/length(hidx);
  factor2=base;
  factor2(nhidx)=ones(size(nhidx))/length(nhidx);;
  cos_matrix=[cos_matrix;factor1];
  sin_matrix=[sin_matrix;factor2];
%  plot(cos_matrix(n,:));
%  hold
%  plot(sin_matrix(n,:),'m');
%  hold off
%  pause
end

periods=[min_per:step_per:max_per];

%calculate cos products
cos_prods=cos_matrix*spectrum;
sin_prods=sin_matrix*spectrum;

%diff prods is the average modulation in dB between a best fit peak and the adjacent trough.
diff_prods=(cos_prods-sin_prods);

%calculate quality of fit for power spectrum
power_spect=exp(spectrum);
p_cos_prods=cos_matrix*power_spect;
p_sin_prods=sin_matrix*power_spect;
percent_harm=(p_cos_prods-p_sin_prods)./p_cos_prods;


ratio_prods=[diff_prods,percent_harm];
[best_ratio, best_idx]=max(ratio_prods(:,1));


%calculate energy at peaks versus average energy
%for i=1:length(best_idx)
% best_cos_prod(i)=cos_prods(best_idx(i),i);
%end
%mean_spect=mean(spectrum);
%peak_vs_avg=best_cos_prod./mean_spect;


best_per=periods(best_idx);
[bad_per_ind]=find(best_per == max_per | best_per == min_per);
best_per(bad_per_ind)=zeros(size(bad_per_ind));

spect_freq_scale=f_step*[1:lspect];
best_freq=best_per.*f_step;
max_freq=max(spect_freq_scale);
ratio_freq_scale=periods*f_step;



if plot_flag ==1
 figure
 subplot(2,1,1);
 plot(ratio_prods);
 subplot(2,1,2)
 plot(spectrum);
 hold on;
 plot(length(hidx)*max(spectrum)*cos_matrix(best_idx,:),'m');
 hold off
end

if plot_flag ==2;
  figure
  subplot(2,1,1);
  plot(ratio_freq_scale,ratio_prods);
  subplot(2,1,2)
  plot(spect_freq_scale,spectrum);
  hold on;
  for i = best_freq:best_freq:max(spect_freq_scale)
    plot([i,i],[0,max(spectrum)],'m')
  end
end

