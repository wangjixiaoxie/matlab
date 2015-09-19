%function gh=plotcov(cent,cov,style,alpha)
function [ghline,ghmarker]=plotcov(cent,cov,style,alpha)

% FUNCTION gh=plotcov(cent,cov,st,alpha)
%
% cent = [x,y] or [x;y]
% cov  is symmetric 2x2
%
% Plot (1-alpha) Conf Interval for Covariance Matrix cov
%%
%% alpha defaults to .05

N=40;
if(nargin<4), alpha=.05; end

if(length(cent)==1),
    cent=[cent,cent];
end

sc=sqrtm(cov);

a=2*pi*[0:N]'/N;
x=sqrt(chi2inv(1-alpha,2))*[cos(a),sin(a)]*sc';

x=[cent(1)+x(:,1),cent(2)+x(:,2)];
if(nargin<3),
    ghmarker=plot(cent(1),cent(2),'o');hold on
    ghline=plot(x(:,1),x(:,2));
else
    if isnumeric(style)
        ghmarker=plot(cent(1),cent(2),'o','color',style);hold on
        ghline=plot(x(:,1),x(:,2),'color',style);
    else
        ghmarker=plot(cent(1),cent(2),[style 'o']);hold on
        ghline=plot(x(:,1),x(:,2),style);
    end
end


% x=[cent(1)+x(:,1),cent(2)+x(:,2)];
% if(nargin<3), ggh=plot(x(:,1),x(:,2));
% else ggh=plot(x(:,1),x(:,2),style);
% end
%if(nargout>0), gh=ggh; end


%%%% THE MAL. DIST IS DIST CHI2 WITH D.O.F = DIMENSIONS (i.e. 2)
%%>
%%> A=rand(2,2);
%%> x=A*randn(2,10000);
%%>
%%> x=[x(1,:)-mean(x(1,:));x(2,:)-mean(x(2,:))];
%%> clf; plot(x(1,:),x(2,:),'.'); hold on; axis equal
%%> C = cov(x'); %% == x*x'/10000
%%>
%%> a = 2*pi*[0:40]/40;
%%> u = [cos(a);sin(a)];
%%>
%%> KK = sqrt(chi2inv(.95,2))
%%> c=KK*sqrtm(C)*u;
%%> set(plot(c(1,:),c(2,:),'r-'),'linew',3)
%%>
%%> mal_distE = sqrt( sum( (inv(sqrtm(C))*x).^2 ) );
%%> mal_distT = sqrt( sum( (inv(sqrtm(A))*x).^2 ) );
%%> mean( mal_distE < KK )   %% ~= .95
%%> mean( mal_distT < KK )   %% ~= .95
%%
%% sqrt(chi2inv(.95,2))  = 2.4477
