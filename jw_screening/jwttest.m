function [h,p,ci,stats] = jwttest(xm,ym,xstd,ystd,xn,yn,alpha,tail)
%TTEST2 Hypothesis test: Compares the averages of two samples.
%   [H,P,CI,STATS] = JWTTEST(XM,YM,XSTD,YSTD,YN,XN,ALPHA,TAIL)
%same as matlab6 buit-in ttest2, but user supplies means and std instead of distribution samples.
%   performs a t-test to
%   determine whether two samples from a normal distribution (with
%   unknown but equal variances) could have the same mean.
%
%   Xm and Ym are sample means.Xstd and Ystd are sample std devs.
%   Xn and Yn are the sample sizes.
%   The null hypothesis is: "means are equal".
%   For TAIL =  0  the alternative hypothesis is: "means are not equal."
%   For TAIL =  1, alternative: "mean of X is greater than mean of Y."
%   For TAIL = -1, alternative: "mean of X is less than mean of Y."
%   TAIL = 0 by default.
%
%   ALPHA is desired significance level (ALPHA = 0.05 by default). 
%   P is the p-value, or the probability of observing the given result
%     by chance given that the null hypothesis is true. Small values
%     of P cast doubt on the validity of the null hypothesis.
%   CI is a confidence interval for the true difference in means.
%   STATS is a structure with two elements named 'tstat' (the value
%     of the t statistic) and 'df' (its degrees of freedom).
%
%   H=0 => "Do not reject null hypothesis at significance level of alpha."
%   H=1 => "Reject null hypothesis at significance level of alpha."

%   References:
%      [1] E. Kreyszig, "Introductory Mathematical Statistics",
%      John Wiley, 1970, section 13.4. (Table 13.4.1 on page 210)

%   Copyright 1993-2000 The MathWorks, Inc. 
%   $Revision: 2.12 $  $Date: 2000/05/26 18:53:51 $

if nargin < 2, 
    error('Requires at least two input arguments'); 
end

%[m1 n1] = size(x);
%[m2 n2] = size(y);
%if (m1 ~= 1 & n1 ~= 1) | (m2 ~= 1 & n2 ~= 1)
%    error('Requires vector first and second inputs.');
%end
%x = x(~isnan(x));
%y = y(~isnan(y));
 
if nargin < 8, 
    tail = 0; 
end 

if nargin < 7, 
    alpha = 0.05; 
end 

if (prod(size(alpha))>1), error('ALPHA must be a scalar.'); end
if (alpha<=0 | alpha>=1), error('ALPHA must be between 0 and 1.'); end

dfx = xn - 1;
dfy = yn - 1;
dfe  = dfx + dfy;
msx = dfx * xstd^2;
msy = dfy * ystd^2;

difference = (xm) - (ym);
pooleds    = sqrt((msx + msy) * (1/(dfx + 1) + 1/(dfy + 1)) / dfe);

ratio = difference / pooleds;
if (nargout>3), stats = struct('tstat', ratio, 'df', dfe); end

% Find the p-value for the tail = 1 test.
p  = 1 - tcdf(ratio,dfe);

% Adjust the p-value for other null hypotheses.
if (tail == 0)
    p = 2 * min(p, 1-p);
    spread = tinv(1 - alpha / 2,dfe) * pooleds;
    if (nargout>2), ci = [(difference - spread) (difference + spread)]; end
else
    spread = tinv(1 - alpha,dfe) * pooleds;
    if (tail == 1)
       if (nargout>2), ci = [(difference - spread), Inf]; end
    else
       p = 1 - p;
       if (nargout>2), ci = [-Inf, (difference + spread)]; end
    end
end

% Determine if the actual significance exceeds the desired significance
h = 0;
if p <= alpha, 
    h = 1; 
end 

if isnan(p), 
    h = NaN; 
end
