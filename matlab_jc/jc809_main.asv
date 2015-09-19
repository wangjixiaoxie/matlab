%------------------------------------------------------------
% The main program
function jc809_main
ma = 9;  % number of parameters
a(ma)=zeros(1,ma); % create the parameters array
alpha=zeros(ma); % create the curvature matrix
beta=zeros(1,ma); % create the beta vector
chisq = 0;
lam = 0;

% start the process
x = 0; y = 0; sig = 0; % are defined as arrays by 'getdata'
jc809_getdata(x,y,sig);
yfit = y; % initialize array of fitted values

% initialize a, alpha, beta, chisq, and lam
init(x,y,sig,a,alpha,beta,chisq,lam)
ochisq = chisq;

% The following loop calls 'mrqmin' and checks for a termination condition
while 1=1
     if mrqmin(x,y,sig,a,alpha,beta,chisq,lam) % true = improvement
   		 % define array of fitted values with the current set of parameters
         % and display the current yfit
         modelf(x,a,yfit)

        % check termination condition
     if chisq == 0 
         break
     end
     if chisq <= ochisq if (ochisq-chisq)/chisq < 1e-40 break
     ochisq = chisq
     
         end
     end
     end
     end
end

% end of process: calculate the covariance matrix
covar = inv(alpha)
j = 1 to ma % index of diagonal elements
vari = covar[j,j] % the variance elements
meanv = mean(vari) % the mean variance
stdv = stdev(vari) % standard deviation of variance