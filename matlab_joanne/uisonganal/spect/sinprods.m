function [cos_matrix, cos_prods, best_val, best_freq]=sinprods(spectrum);

min_per=5;
step_per=.1;
max_per=20;


%length of spect
lspect=size(spectrum,1);
cos_matrix=[];
sin_matrix=[];
n=0;
% factors of appropriate periodicities
indexes=[0:lspect-1];
for i=min_per:step_per:max_per
  n=n+1;
  factor=cos(2*pi/i*indexes);
  %truncate at integral number of cycles
  n_points=min(lspect,ceil(i*floor((lspect)/i)));
  factor1=[factor(1,1:n_points),-1*zeros(1,lspect-n_points)];
  factor2=[factor(1,1:n_points),ones(1,lspect-n_points)];
  %factor1=(factor1+1)/2;
  factor2=1-(factor2+1)/2;
  cos_matrix=[cos_matrix;factor1];
  sin_matrix=[sin_matrix;factor2];
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

[best_val, best_idx]=max(ratio_prods);

best_freq=periods(best_idx);

plot(cos_prods);
figure
plot(spectrum);
hold on
plot(max(spectrum)*cos_matrix(best_cos,:),'m');
hold off


