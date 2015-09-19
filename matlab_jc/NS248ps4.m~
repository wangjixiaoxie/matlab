
% Problem 1a
N=10;
[X,Y]=GetSample(N);
mean(X)-mean(Y)
% -0.7003

B=1000;
Dboot=zeros(B,1);
for b=1:B
    Xboot=X(ceil(N*rand(N,1)));
    Yboot=Y(ceil(N*rand(N,1)));
    Dboot(b)=mean(Xboot)-mean(Yboot);
end
std(Dboot) = 0.6248 
% -0.7003 is between 2*0.6248 of zero, thus not significant

prctile(Dboot,[2.5 97.5])=-1.9402 & 0.5554 
%-- spans zero thus not significant
% No evidence of adaptation.

% Problem 1b
B=1000;
Dboot=zeros(B,1);
for b=1:B
    a=ceil(N*rand(N,1));
    Xboot=X(ceil(a));
    Yboot=Y(ceil(a));
    Dboot(b)=mean(Xboot)-mean(Yboot);
end
std(Dboot) = 0.1525 
% -0.7003 is between 2*0.1525 of zero, thus significant

prctile(Dboot,[2.5 97.5])=-1.0167 & -0.3986 
%-- does not span zero thus significant
% CHANGES THE INTERPRETATION

% Problem 2a
D=mean(X)-mean(Y);
Dperm=zeros(B,1);
Z=[X;Y];
for r=1:B
    [tmp,i]=sort(rand(N+N,1));
    Zperm=Z(i,:);
    Dperm(r)=mean(Zperm(1:N,:))-mean(Zperm(N+1:2*N,:));
end
p=mean(abs(Dperm)>abs(D))=0.2960
% not significant

% Problem 2b
Dperm=zeros(B,1);
Z=[X Y];
for r=1:B
    xind=ceil(rand(1,N)*2); % x is either 1 or 2
    yind=3-xind; % y is either 2 or 1 (whichever x is not)
    for k=1:N
        xvals(k)=Z(k,xind(k));
        yvals(k)=Z(k,yind(k));
    end
    Dperm(r)=mean(xvals)-mean(yvals);
end
p=mean(abs(Dperm)>abs(D))=0.001 % significant
% YES it does change the results of the test

% Problem 3
% UNPAIRED
Dvals=[0.12 0.25 0.5 1 2];
B=500;
for N=[10 20]
    for D=1:5
        D
        guess=zeros(1,200);
        for j=1:200
            [X,Y]=GetSample(N,Dvals(D));
            Dperm=zeros(B,1);
            Z=[X;Y];
            for r=1:B
                [tmp,i]=sort(rand(N+N,1));
                Zperm=Z(i,:);
                Dperm(r)=mean(Zperm(1:N,:))-mean(Zperm(N+1:2*N,:));
            end
            guess(j)=mean(abs(Dperm)>abs(Dvals(D)))<0.05; % correctly reject H0
        end
        PowerEst(N/10,D)=(sum(guess))/200; % proportion of time that you reject H0
    end
end
figure;plot(Dvals,PowerEst')
% PAIRED
clear all
Dvals=[0.12 0.25 0.5 1 2];
B=500;
for N=[10 20]
    for D=1:5
        D
        guess=zeros(1,200);
        for j=1:200
            [X,Y]=GetSample(N,Dvals(D));
            Dperm=zeros(B,1);
            Z=[X Y];
            for r=1:B
                xind=ceil(rand(1,N)*2); % x is either 1 or 2
                yind=3-xind; % y is either 2 or 1 (whichever x is not)
                for k=1:N
                    xvals(k)=Z(k,xind(k));
                    yvals(k)=Z(k,yind(k));
                end
                Dperm(r)=mean(xvals)-mean(yvals);
            end
            guess(j)=mean(abs(Dperm)>abs(Dvals(D)))<0.05; % correctly reject H0
        end
        PowerEst(N/10,D)=(sum(guess))/200; % proportion of time that you reject H0
    end
end


hold on;plot(Dvals,PowerEst','Linewidth',2)


