
% Problem 1
% A
    % Data are PAIRED
    [h,p]=ttest(y2-y1)
    % p=0.0079
    
    [p,h]=signtest(y2-y1)
    % p=0.0354
    
    % Significant even if you don't assume normality
% B
r=corrcoef(x,y1)
% r(2)^2=0.6179
r=corrcoef(x,y2)
% r(2)^2=0.5309
mean([0.6179 0.5309])=0.5744 % a little more than half
% Other variables ---> arousal/attention state of the monkey, biophysical variation
% in visual processing, uncontrolled properties of the stimulus...

% Problem 2
% A - figure1
figure;plot([1:40],g1,'*','Color','b')
hold on;plot([1:40],g2,'*','Color','r')


% B - figure2
figure;errorbar(mean(g1),std(g1)/sqrt(10-1))
hold on;errorbar(mean(g2),std(g2)/sqrt(10-1),'r')
% Significant at higher values
% It is questionable whether normality is appropriate here - data are
% skewed towards higher values.

% C - figure 3
latency = [g1;g2];
animalspergroup=10;
[p,table,stats]=anova2(latency,animalspergroup);
% Both columns (trial number) and rows (group) have a highly significant
% effect on latency.

%% multcompare is strange --- doesn't do what I want?

% D 
x=1:1:40;
figure;hold on;plot(mean(g1)','.');plot(mean(g2)','.','Color','r')
[b1,bint1,r1,rint1,stats1]=regress([mean(g1)'],[ones(40,1),x'],0.05);
[b2,bint2,r2,rint2,stats2]=regress([mean(g2)'],[ones(40,1),x'],0.05);
plot(x,b1(1)+b1(2)*x','Linewidth',3);
plot(x,b2(1)+b2(2)*x','r','Linewidth',3);
plot(x,(bint1(1)+bint1(2)*x));plot(x,(bint1(3)+bint1(4)*x))
plot(x,(bint2(1)+bint2(2)*x),'r');plot(x,(bint2(3)+bint2(4)*x),'r')

for i=1:400
x(i)=ceil(i/10);
y1(i)=g1(i);
y2(i)=g2(i);
end
[b1,bint1,r1,rint1,stats]=regress(y1',[ones(400,1) x'],0.05);
[b2,bint2,r2,rint2,stats]=regress(y2',[ones(400,1) x'],0.05);
figure;plot(x,b1(1)+b1(2)*x','b','Linewidth',3);hold on;plot(x,y1,'.','Color','b')
plot(x,b2(1)+b2(2)*x','r','Linewidth',3);plot(x,y2,'.','Color','r')
plot(x,(bint1(1)+bint1(2)*x));plot(x,(bint1(3)+bint1(4)*x))
plot(x,(bint2(1)+bint2(2)*x),'r');plot(x,(bint2(3)+bint2(4)*x),'r')

% Analysis does not indicate a group difference

% E
% Fit with an exponential because that is appropriate for a learning curve
% and because data are U-shape curved around the linear fit.

% Use nlinfit - requires a function and a guess
% function (exponential)
    function y=expfit(x,a,b,c)
    y=a*exp(-b*x)+c;
% guess - use same guess for both to avoid bias
    y=45*exp(-0.04*[1:1:40])+10;
% nlinfit
[beta1,r1,j1]=nlinfit(x,y1,@expfit,[45 0.04 10]);
[beta2,r2,j2]=nlinfit(x,y2,@expfit,[45 0.04 10]);
sum(r1.^2)=7.4994e+04
sum(r2.^2)=1.0550e+05
ci1=nlparci(beta1,r1,j1,0.05);
ci2=nlparci(beta2,r2,j2,0.05);
figure;hold on;
plot(x,y1,'.','Color','b')
plot(x,expfit(beta1,x),'Linewidth',3)
plot(x,expfit(ci1(:,1),x),'Linewidth',1)
plot(x,expfit(ci1(:,2),x),'Linewidth',1)
plot(x,y2,'.','Color','r')
plot(x,expfit(beta2,x),'Linewidth',3,'Color','r')
plot(x,expfit(ci2(:,1),x),'Linewidth',1,'Color','r')
plot(x,expfit(ci2(:,2),x),'Linewidth',1,'Color','r')


%%%% to do...
% Fit with an exponential because that is appropriate for a learning curve
% and because data are U-shape curved around the linear fit.

% Use nlinfit - requires a function and a guess
% function (exponential)
    function y=expfit(x,a,b,c)
    y=a*exp(-b*x)+c;
% guess - use same guess for both to avoid bias
    y=45*exp(-0.04*[1:1:40])+10;
% nlinfit
y1=mean(g1)';
y2=mean(g2)';
[beta1,r1,j1]=nlinfit(x,y1',@expfit,[45 0.04 10]);
[beta2,r2,j2]=nlinfit(x,y2',@expfit,[45 0.04 10]);
sum(r1.^2)=7.4994e+04
sum(r2.^2)=1.0550e+05
ci1=nlparci(beta1,r1,j1,0.05);
ci2=nlparci(beta2,r2,j2,0.05);
figure;hold on;
plot(x,y1,'.','Color','b')
plot(x,expfit(beta1,x),'Linewidth',3)
plot(x,expfit(ci1(:,1),x),'Linewidth',1)
plot(x,expfit(ci1(:,2),x),'Linewidth',1)
plot(x,y2,'.','Color','r')
plot(x,expfit(beta2,x),'Linewidth',3,'Color','r')
plot(x,expfit(ci2(:,1),x),'Linewidth',1,'Color','r')
plot(x,expfit(ci2(:,2),x),'Linewidth',1,'Color','r')



% Problem 3
% Can the spiketimes predict the stimulus? i.e. act as a filter
% A
for i=1:50
    for j=1:30
        spikes1(i,j)=sum(spiketimes1{i}<j*0.01);
        spikes2(i,j)=sum(spiketimes2{i}<j*0.01);
    end
end
spikes1(:,2:30)=spikes1(:,2:30)-spikes1(:,1:29);
spikes2(:,2:30)=spikes2(:,2:30)-spikes2(:,1:29);
mnspikes1=mean(spikes1);
serrspikes1=std(spikes1)/sqrt(50-1);
mnspikes2=mean(spikes2);
serrspikes2=std(spikes2)/sqrt(50-1);
figure;hold on;
plot(mnspikes1,'Linewidth',3);plot(mnspikes2,'Linewidth',3,'Color','r')
plot(mnspikes1+serrspikes1);plot(mnspikes1-serrspikes1)
plot(mnspikes2+serrspikes2,'r');plot(mnspikes2-serrspikes2,'r')
% Might be a difference in the magnitude of the response to song onset.
% Response to blue (non-tutor song) is greater

% Something is wrong with errorbar.m - not sure why - it was okay saturday

% B
% x is firing rate (Hz)
% y is stimulus --- 1 or 2
y=[ones(1,50) 2*ones(1,50)];
x=[mean(spikes1')/0.3 mean(spikes2')/0.3];
[b,dev,stats]=glmfit(x,y,'poisson');
stats.p(2)=0.6526 % - nowhere near significant
% C - 10ms bins
y=[ones(1,50) 2*ones(1,50)]';
x=[spikes1;spikes2];
[b,dev,stats]=glmfit(x,y,'poisson');
stats.p(2)=0.9038 % - nowhere near significant
% D ????