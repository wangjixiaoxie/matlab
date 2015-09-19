

function [W,M,V,L,E] = jc_EM_GM(X,k,ltol,maxiter,pflag,Init)
% Inputs:
% X(n,d) - input data, n=number of observations, d=dimension of variable
% k - maximum number of Gaussian components allowed
% ltol - percentage of the log likelihood difference between 2 iterations ([] for none)
% maxiter - maximum number of iteration allowed ([] for none)
% pflag - 1 for plotting GM for 1D or 2D cases only, 0 otherwise ([] for none)
% Init - structure of initial W, M, V: Init.W, Init.M, Init.V ([] for none)
%
% Ouputs:
% W(1,k) - estimated weights of GM
% M(d,k) - estimated mean vectors of GM
% V(d,d,k) - estimated covariance matrices of GM
% L - log likelihood of estimates
%

%%%% Validate inputs %%%%
if nargin <= 1,
disp('EM_GM must have at least 2 inputs: X,k!/n')
return
elseif nargin == 2,
ltol = 0.1; maxiter = 1000; pflag = 0; Init = [];
err_X = Verify_X(X);
err_k = Verify_k(k);
if err_X | err_k, return; end
elseif nargin == 3,
maxiter = 1000; pflag = 0; Init = [];
err_X = Verify_X(X);
err_k = Verify_k(k);
[ltol,err_ltol] = Verify_ltol(ltol);
if err_X | err_k | err_ltol, return; end
elseif nargin == 4,
pflag = 0; Init = [];
err_X = Verify_X(X);
err_k = Verify_k(k);
[ltol,err_ltol] = Verify_ltol(ltol);
[maxiter,err_maxiter] = Verify_maxiter(maxiter);
if err_X | err_k | err_ltol | err_maxiter, return; end
elseif nargin == 5,
Init = [];
err_X = Verify_X(X);
err_k = Verify_k(k);
[ltol,err_ltol] = Verify_ltol(ltol);
[maxiter,err_maxiter] = Verify_maxiter(maxiter);
[pflag,err_pflag] = Verify_pflag(pflag);
if err_X | err_k | err_ltol | err_maxiter | err_pflag, return; end
elseif nargin == 6,
err_X = Verify_X(X);
err_k = Verify_k(k);
[ltol,err_ltol] = Verify_ltol(ltol);
[maxiter,err_maxiter] = Verify_maxiter(maxiter);
[pflag,err_pflag] = Verify_pflag(pflag);
[Init,err_Init]=Verify_Init(Init);
if err_X | err_k | err_ltol | err_maxiter | err_pflag | err_Init, return; end
else
disp('EM_GM must have 2 to 6 inputs!');
return
end

%%%% Initialize W, M, V,L %%%%
t = cputime;
if isempty(Init),
[W,M,V] = Init_EM(X,k); L = 0;
else
W = Init.W;
M = Init.M;
V = Init.V;
end
[E,Ln] = Expectation(X,k,W,M,V); % Initialize log likelihood and compute expectation
Lo = 2*Ln;

%%%% EM algorithm %%%%
niter = 0;
while (abs(100*(Ln-Lo)/Lo)>ltol) & (niter<=maxiter),
%**************************************************************************
% For the first loop, expectation is computed above (cf. Line 69). For the
% next loops, it is computed in their preceding loops. Therefore, the
% following line [Line 80] is omitted.
%**************************************************************************
% E = Expectation(X,k,W,M,V); % E-step

[W,M,V] = Maximization(X,k,E); % M-step
Lo = Ln;
[E,Ln] = Expectation(X,k,W,M,V); % Log-likelihood and expectation computation
niter = niter + 1;
LnValues(niter) = Ln;
disp(['Iteration # = ' num2str(niter) ', Ln = ' num2str(Ln) ', Stopping Criterion = ' num2str(abs(100*(Ln-Lo)/Lo)) ', Tolerance = ' num2str(ltol)])
end
L = Ln;
theTime = num2str(fix(clock));
theTime(theTime==' ') = '_';
save(['EM_CostFunction_' theTime '_MoreLoopsRemoved.mat'],'LnValues')

%%%% Plot 1D or 2D %%%%
if pflag==1,
[n,d] = size(X);
if d>2,
disp('Can only plot 1 or 2 dimensional applications!/n');
else
Plot_GM(X,k,W,M,V);
end
elapsed_time = sprintf('CPU time used for EM_GM: %5.2fs',cputime-t);
disp(elapsed_time);
disp(sprintf('Number of iterations: %d',niter-1));
end
%%%%%%%%%%%%%%%%%%%%%%
%%%% End of EM_GM %%%%
%%%%%%%%%%%%%%%%%%%%%%





%**************************************************************************
%************************ LIKELIHOOD REMOVED COMPLETELY *******************
%**************************************************************************
% function L = Likelihood(X,k,W,M,V)
% Compute L based on K. V. Mardia, "Multivariate Analysis", Academic Press, 1979, PP. 96-97
% % to enchance computational speed
% [n,d] = size(X);
% U = mean(X)';
% S = cov(X);
% L = 0;
% for i=1:k,
% iV = inv(V(:,:,i));
% L = L + W(i)*(-0.5*n*log(det(2*pi*V(:,:,i))) ...
% -0.5*(n-1)*(trace(iV*S)+(U-M(:,i))'*iV*(U-M(:,i))));
% end
%**************************************************************************
%*************************** MODIFICATION ENDS HERE ***********************
%**************************************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% End of Likelihood %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%












