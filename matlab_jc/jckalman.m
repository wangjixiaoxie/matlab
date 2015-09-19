function [mn,se]=jckalman(dataIN)
% This takes a 1xn vector of FF values (or anything) and uses linear
% filtering to estimate the mean and standard error of the mean.
% Uses code from Zoubin Gharahmani
% http://www.gatsby.ucl.ac.uk/~zoubin/software.html
% under "EM for linear dynamical systems". 
% Which Anne Smith pointed me to.
% Also acknowledge Loren Frank.


% Pre-processing
    % Subtract baseline median from data
    % Append baseline data to beginning of data
    % Append zeros(1,50) to beginning of data

% make a 1xn vector
data(1,1:length(dataIN))=dataIN(1:end);

% Data order of 1 - since x is a 1xn vector of FF values
dataorder=1;

% Generate a state-space model of the data
net=lds(data',dataorder,length(data),1000,0.0001);

% Use state-space model to smooth the data
y(1,1,1:length(data))=data;
[lik,Xfin,Pfin]=kalmansmooth(net.A,net.C,net.Q,net.R,net.x0,net.P0,y);

% Xfin is the state
% Pfin is the variance of the state estimate
    % thus sqrt(Pfin) is the standard error of the state estimate
for i=1:length(data)
    mn(i)=Xfin(:,:,i)*net.C;
    varmn=Pfin(:,:,i);
    se(i)=sqrt(varmn)*net.C;
end

