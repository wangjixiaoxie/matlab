function [cos_matrix, cos_prods, best_val, best_freq]=deltaprods2(spectrum);

%Note: function assumes that first bin of spectrum represents a frequency of f_step, NOT 0.

min_per=3;
step_per=.1;
max_per=15;


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
  factor1(hidx)=ones(size(hidx));
  factor2=base;
  factor2(nhidx)=ones(size(nhidx));
  cos_matrix=[cos_matrix;factor1];
  sin_matrix=[sin_matrix;factor2];
  %plot(cos_matrix(n,:));
  %hold
  %plot(sin_matrix(n,:),'m');
  %hold off
  %pause
end

periods=[min_per:step_per:max_per];

%calculate cos products
cos_prods=cos_matrix*spectrum;
sin_prods=sin_matrix*spectrum;
ratio_prods=cos_prods./sin_prods;
diff_prods=cos_prods-sin_prods;

[best_val, best_idx]=max(ratio_prods);

best_freq=periods(best_idx);
figure
subplot(3,1,1);
plot(ratio_prods);
subplot(3,1,2);
plot(diff_prods);
subplot(3,1,3)
plot(spectrum);
hold on
plot(max(spectrum)*cos_matrix(best_idx,:),'m');
hold off


