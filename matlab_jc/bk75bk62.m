dirf('*.cbin','batch')
findcatch('batch')
evsonganaly

load /cardinal/bk75bk62/seqData050510.mat
tvall914wn=[tv914wnB tv914wnC];
numsall2=[zeros(1,length([tv914wnB])) ones(1,length([tv914wnC]))];
[b2,ind2]=sort(tvall914wn);
numsorted914wn=numsall2(ind2);
tvsorted914wn=tvall914wn(ind2);

tvall914pre=[tv914preB tv914preC];
numsall2=[zeros(1,length([tv914preB])) ones(1,length([tv914preC]))];
[b2,ind2]=sort(tvall914pre);
numsorted914pre=numsall2(ind2);
tvsorted914pre=tvall914pre(ind2);
% probability of B as a function of time
figure;hold on;
subplot(211);hold on;
plot(runningaverage(tvsorted914pre,40),1-runningaverage(numsorted914pre,40),'*','Color','k')
plot(runningaverage(tvsorted914wn(1:450),40),1-runningaverage(numsorted914wn(1:450),40),'*','Color','r')
subplot(212);hold on;
plot(runningaverage(tvsorted1001pre,40),1-runningaverage(numsorted1001pre,40),'*','Color','k')
plot(runningaverage(tvsorted1001apv,40),1-runningaverage(numsorted1001apv,40),'*','Color','r')



% dominant A
fvalsDOM=findwnoteJC('batchJCnotes','a','b','',0,[3000 4200],8500,1,'obs0',1);
% recessive A
fvalsREC=findwnoteJC('batchJCnotes','a','cc','',0,[3000 4200],8500,1,'obs0',1);

size(fvalsDOM,2)/(size(fvalsDOM,2)+size(fvalsREC,2))
% 63.3percent dominance

% sequence2.m

%%% Pre-surgery experiment
% 9/14/09 at 11:12am chick time - WN on
% 9/14/09 at 5:14pm chick time - WN off


% Surgery on 9/17/09
% mass=13.2g
% 2200um to RA on right
% placement ~1250um

% probes in on 9/23/09
% singing began on 9/25/09 around noon
% 9/27/09 at 1:44pm chick time - 2mM apv on at 1mL/min
% new template at 3:20pm chick time - hit above 2640Hz
% 6:40pm - 2630Hz
% 7:47pm - acsf on at 1.5uL/min
% 8pm - acsf on at 0.6uL/min
% 8:21pm - door closed /wn off/ acq on

% 9/30/09 at 1:30pm - switched to syntax program

% 10/01/09 at 9:00am - AP5 on / WN on - okay b/c didn't sing until 9:15am
    % and didn't sing through for another 10minutes
% 10/01/09 at 3:02pm - ACSF on at 1.5uL/min and WN off and door open and
    % acquisition off
   % at 3:17pm - ACSF on at 0.6uL/min
   % sang with door open
% 10/01/09 at 3:37pm - door closed, acquisition on, new folder
% 10/01/09 at 3:50-4pm - door open but singing
% 10/01/09 at 4pm - door wide open - no singing
% 10/01/09 at 5:40pm - door closed

% 10/02/09 at 11:00am chick time - APV on at 1.0uL/min
% 10/02/09 at 12:38pm - WN on - hit below 2610Hz
%              1:35pm - hit below 2580Hz
%              2:38pm - hit below 2600Hz
% 10/02/09 at 5:35pm - WN off - door open - ACSF at 1.5uL/min
% 10/03/09 at ~7:40am - clog detected - pressure off

% 10/04/09 at ~8pm - bird deceased - put into fix within 2hrs by TW

figure;plot(tvals927,mean(pitch927(320:380,:)),'.','MarkerSize',5)
hold on;plot(tvals927apv,mean(pitch927apv(320:380,:)),'.','MarkerSize',10,'Color','r')
hold on;plot(tvals927apvwn,mean(pitch927apvwn(320:380,:)),'.','MarkerSize',10,'Color','g')
hold on;plot(tvals928acsf,mean(pitch928acsf(320:380,:)),'.','MarkerSize',10,'Color','b')
hold on;plot(tvals929acsf,mean(pitch929acsf(320:380,:)),'.','MarkerSize',10,'Color','b')
hold on;plot(tvals930acsf,mean(pitch930acsf(320:380,:)),'.','MarkerSize',10,'Color','b')
hold on;plot(tvals930Bacsf,mean(pitch930Bacsf(320:380,:)),'.','MarkerSize',10,'Color','b')
hold on;plot(tvals1001apv,mean(pitch1001apv(320:380,:)),'.','MarkerSize',10,'Color','r')
hold on;plot(tvals1001acsf,mean(pitch1001acsf(320:380,:)),'.','MarkerSize',10,'Color','b')
hold on;plot(tvals1002apv,mean(pitch1002apv(320:380,:)),'.','MarkerSize',10,'Color','r')
hold on;plot(tvals1002apvwn,mean(pitch1002apvwn(320:380,:)),'.','MarkerSize',10,'Color','g')
hold on;plot(tvals1003acsf,mean(pitch1003acsf(320:380,:)),'.','MarkerSize',10,'Color','b')

% relative
fvals928acsfCCA=findwnoteJC('batchJCnotes','a','c-','',0,[2000 2700],8500,1,'obs0',1);
fvals928acsf=findwnoteJC('batchJCnotes','c','','-a',0,[2000 2700],8500,1,'obs0',1);
fvals930acsfCCA=findwnoteJC('batchJCnotes','a','c-','',0,[2000 2700],8500,1,'obs0',1);
fvals930acsf=findwnoteJC('batchJCnotes','c','','-a',0,[2000 2700],8500,1,'obs0',1);
for i=1:length(fvals928acsf)
    shifted928acsf(i,:)=fvals928acsf(i).datt;
    shifted928acsfCCA(i,:)=fvals928acsfCCA(i).datt;
end
for i=1:length(fvals930acsf)
    shifted930acsf(i,:)=fvals929acsf(i).datt;
    shifted930acsfCCA(i,:)=fvals929acsfCCA(i).datt;
end
pitch1002acsfCCA=jc_pitchmat1024(shifted1002acsfCCA,1024,1020,1,2800,3800,[1],'obs0',1);
pitch1002acsf=jc_pitchmat1024(shifted1002acsf,1024,1020,1,2100,2800,[1],'obs0',1);

figure;plot(tvals927(31:end),mean(pitch927(320:380,31:end))-mean(pitch927CCA(250:350,31:end)),'.','MarkerSize',5)
hold on;plot(tvals927apv,mean(pitch927apv(320:380,:))-mean(pitch927apvCCA(250:350,:)),'.','MarkerSize',10,'Color','r')
hold on;plot(tvals927apvwn,mean(pitch927apvwn(320:380,:))-mean(pitch927apvwnCCA(250:350,:)),'.','MarkerSize',10,'Color','g')
hold on;plot(tvals928acsf,mean(pitch928acsf(320:380,:))-mean(pitch928acsfCCA(250:350,:)),'.','MarkerSize',10,'Color','b')
hold on;plot(tvals929acsf,mean(pitch929acsf(320:380,:))-mean(pitch929acsfCCA(250:350,:)),'.','MarkerSize',10,'Color','b')
hold on;plot(tvals930acsf,mean(pitch930acsf(320:380,:))-mean(pitch930acsfCCA(250:350,:)),'.','MarkerSize',10,'Color','b')
hold on;plot(tvals930Bacsf,mean(pitch930Bacsf(320:380,:))-mean(pitch930BacsfCCA(250:350,:)),'.','MarkerSize',10,'Color','b')
hold on;plot(tvals1001apv,mean(pitch1001apv(320:380,:))-mean(pitch1001apvCCA(250:350,:)),'.','MarkerSize',10,'Color','r')
hold on;plot(tvals1001acsf,mean(pitch1001acsf(320:380,:))-mean(pitch1001acsfCCA(250:350,:)),'.','MarkerSize',10,'Color','b')
hold on;plot(tvals1002apv,mean(pitch1002apv(320:380,:))-mean(pitch1002apvCCA(250:350,:)),'.','MarkerSize',10,'Color','r')
hold on;plot(tvals1002apvwn,mean(pitch1002apvwn(320:380,:))-mean(pitch1002apvwnCCA(250:350,:)),'.','MarkerSize',10,'Color','g')
hold on;plot(tvals1003acsf,mean(pitch1003acsf(320:380,:))-mean(pitch1003acsfCCA(250:350,:)),'.','MarkerSize',10,'Color','b')

figure;plot(tvals927(31:end),mean(pitch927(320:380,31:end))-mean(pitch927(1200:1450,31:end)),'.','MarkerSize',5)
hold on;plot(tvals927apv,mean(pitch927apv(320:380,:))-mean(pitch927apv(1200:1450,:)),'.','MarkerSize',10,'Color','r')
hold on;plot(tvals927apvwn,mean(pitch927apvwn(320:380,:))-mean(pitch927apvwn(1200:1450,:)),'.','MarkerSize',10,'Color','g')
hold on;plot(tvals928acsf,mean(pitch928acsf(320:380,:))-mean(pitch928acsf(1200:1450,:)),'.','MarkerSize',10,'Color','b')
hold on;plot(tvals929acsf,mean(pitch929acsf(320:380,:))-mean(pitch929acsf(1200:1450,:)),'.','MarkerSize',10,'Color','b')
hold on;plot(tvals930acsf,mean(pitch930acsf(320:380,:))-mean(pitch930acsf(1200:1450,:)),'.','MarkerSize',10,'Color','b')
hold on;plot(tvals930Bacsf,mean(pitch930Bacsf(320:380,:))-mean(pitch930Bacsf(1200:1450,:)),'.','MarkerSize',10,'Color','b')
hold on;plot(tvals1001apv,mean(pitch1001apv(320:380,:))-mean(pitch1001apv(1200:1450,:)),'.','MarkerSize',10,'Color','r')
hold on;plot(tvals1001acsf,mean(pitch1001acsf(320:380,:))-mean(pitch1001acsf(1200:1450,:)),'.','MarkerSize',10,'Color','b')
hold on;plot(tvals1002apv,mean(pitch1002apv(320:380,:))-mean(pitch1002apv(1200:1450,:)),'.','MarkerSize',10,'Color','r')
hold on;plot(tvals1002apvwn,mean(pitch1002apvwn(320:380,:))-mean(pitch1002apvwn(1200:1450,:)),'.','MarkerSize',10,'Color','g')
hold on;plot(tvals1003acsf,mean(pitch1003acsf(320:380,:))-mean(pitch1003acsf(1200:1450,:)),'.','MarkerSize',10,'Color','b')


%%% justification
clear instACT
clear instREL
width=50;
for i=1:size(pitch913,2)-width
    instACT(i)=median(mean(pitch913(320:380,i:i+width)));
    instREL(i)=median(mean(pitch913(320:380,i:i+width)))-median(mean(pitch913(1200:1450,i:i+width)));
end
std(instACT)
std(instREL)

