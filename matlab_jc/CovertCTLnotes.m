a=150;
% pu57
% note after instead? YES YES YEST
figure;subplot(231);hold on;plot(Experiment(1).timeACpre,mean(Experiment(1).pitchACpre(154:216,:))-mean(Experiment(1).pitchACpre(1000:1050,:)),'*')
%plot(Experiment(1).timeAPV,mean(Experiment(1).pitchAPV(154:216,:))-mean(Experiment(1).pitchAPV(1000:1050,:)),'*','Color','r')
plot(Experiment(1).timeAPVwn,mean(Experiment(1).pitchAPVwn(154:216,:))-mean(Experiment(1).pitchAPVwn(1000:1050,:)),'*','Color','r')
plot(Experiment(1).timeACpost,mean(Experiment(1).pitchACpost(154:216,:))-mean(Experiment(1).pitchACpost(1000:1050,:)),'*')
hold on;plot(Experiment(2).timeACpre,mean(Experiment(2).pitchACpre(170:234,:))-mean(Experiment(2).pitchACpre(1000:1050,:)),'*')
%plot(Experiment(2).timeAPV,mean(Experiment(2).pitchAPV(170:234,:))-mean(Experiment(2).pitchAPV(1000:1050,:)),'*','Color','r')
plot(Experiment(2).timeAPVwn,mean(Experiment(2).pitchAPVwn(170:234,:))-mean(Experiment(2).pitchAPVwn(1000:1050,:)),'*','Color','r')
plot(Experiment(2).timeACpost,mean(Experiment(2).pitchACpost(170:234,:))-mean(Experiment(2).pitchACpost(1000:1050,:)),'*')
clear instACT
clear instREL
width=50;
for i=1:size(Experiment(1).pitchACpre,2)-width
    instACT(i)=median(mean(Experiment(1).pitchACpre(154:216,i:i+width)));
    instREL(i)=median(mean(Experiment(1).pitchACpre(154:216,i:i+width)))-median(mean(Experiment(1).pitchACpre(1000:1050,i:i+width)));
end
std(instACT)
std(instREL)



% b10o7
subplot(234);hold on;plot(E7timepre,mean(pitchE7pretargeted(150:170,:))-mean(pitchE7pre(300:350,:)),'*')
%hold on;plot(Experiment(7).timeAPV(1:51),mean(Experiment(7).pitchAPV(150:170,1:51))-mean(pitchE7apv(300:350,:)),'*','Color','r')
plot(Experiment(7).timeAPVwn(1:110),mean(Experiment(7).pitchAPVwn(150:170,1:110))-mean(pitchE7apvwn(300:350,:)),'*','Color','r')
plot(Experiment(7).timeACpost,mean(Experiment(7).pitchACpost(150:170,:))-mean(pitchE7post(300:350,:)),'*')

hold on;plot(Experiment(7).timeACpost,mean(Experiment(7).pitchACpost(150:170,:))-mean(pitchE8pre(300:350,:)),'*')
hold on;plot(Experiment(8).timeAPV,mean(Experiment(8).pitchAPV(150:170,:))-mean(pitchE8apv(300:350,:)),'*','Color','r')
%%% justification
clear instACT
clear instREL
width=50;
for i=1:size(pitchE7pretargeted,2)-width
    instACT(i)=median(mean(pitchE7pretargeted(150:170,i:i+width)));
    instREL(i)=median(mean(pitchE7pretargeted(150:170,i:i+width)))-median(mean(pitchE7pre(300:350,i:i+width)));
end
std(instACT)
std(instREL)

% o30bk18
subplot(235);hold on;plot(Experiment(10).timeACpre,mean(Experiment(10).pitchACpre(150:204,:))-mean(Experiment(10).pitchACpre(850:1100,:)),'*')
%plot(Experiment(10).timeAPV,mean(Experiment(10).pitchAPV(150:204,:))-mean(Experiment(10).pitchAPV(850:1100,:)),'*','Color','r')
plot(Experiment(10).timeAPVwn,mean(Experiment(10).pitchAPVwn(150:204,:))-mean(Experiment(10).pitchAPVwn(850:1100,:)),'*','Color','r')
plot(Experiment(10).timeACpost,mean(Experiment(10).pitchACpost(150:204,:))-mean(Experiment(10).pitchACpost(850:1100,:)),'*')
hold on;plot(Experiment(11).timeACpre,mean(Experiment(11).pitchACpre(150:204,:))-mean(Experiment(11).pitchACpre(850:1100,:)),'*')
plot(Experiment(11).timeAPV,mean(Experiment(11).pitchAPV(150:204,:))-mean(Experiment(11).pitchAPV(850:1100,:)),'*','Color','r')
clear instACT
clear instREL
width=50;
for i=1:size(Experiment(10).pitchACpre,2)-width
    instACT(i)=median(mean(Experiment(10).pitchACpre(150:204,i:i+width)));
    instREL(i)=median(mean(Experiment(10).pitchACpre(150:204,i:i+width)))-median(mean(Experiment(10).pitchACpre(850:1100,i:i+width)));
end
std(instACT)
std(instREL)


% bk75bk62
subplot(236);plot(tvals927(31:end),mean(pitch927(320:380,31:end))-pitch927CCA(350,31:end),'*')
hold on;plot(tvals927apv,mean(pitch927apv(320:380,:))-pitch927apvCCA(350,:),'*','Color','r')
hold on;plot(tvals927apvwn,mean(pitch927apvwn(320:380,:))-pitch927apvwnCCA(350,:),'*','Color','g')
hold on;plot(tvals928acsf,mean(pitch928acsf(320:380,:))-pitch928acsfCCA(350,:),'*','Color','b')
hold on;plot(tvals929acsf,mean(pitch929acsf(320:380,:))-pitch929acsfCCA(350,:),'*','Color','b')


% pu57
% GOOD
figure;hold on;plot(Experiment(1).timeACpre,mean(Experiment(1).pitchACpre(154:216,:))-pitchE1pre(a,:),'*')
plot(Experiment(1).timeAPV,mean(Experiment(1).pitchAPV(154:216,:))-pitchE1apv(a,:),'*','Color','r')
plot(Experiment(1).timeAPVwn,mean(Experiment(1).pitchAPVwn(154:216,:))-pitchE1apvwn(a,1:449),'*','Color','g')
plot(Experiment(1).timeACpost(1:890),mean(Experiment(1).pitchACpost(154:216,1:890))-pitchE1post(a,1:890),'*')
% slight shift at the end excludes it from consideration
figure;hold on;plot(Experiment(2).timeACpre,mean(Experiment(2).pitchACpre(170:234,:))-mean(pitchE2pre(150:200,:)),'*')
plot(Experiment(2).timeAPV,mean(Experiment(2).pitchAPV(170:234,:))-mean(pitchE2apv(150:200,:)),'*','Color','r')
plot(Experiment(2).timeAPVwn,mean(Experiment(2).pitchAPVwn(170:234,:))-mean(pitchE2apvwn(150:200,:)),'*','Color','g')
plot(Experiment(2).timeACpost(1:1120),mean(Experiment(2).pitchACpost(170:234,1:1120))-mean(pitchE2post(150:200,1:1120)),'*')
% note after instead? YES YES YEST
figure;hold on;plot(Experiment(1).timeACpre,mean(Experiment(1).pitchACpre(154:216,:))-mean(Experiment(1).pitchACpre(1000:1050,:)),'*')
plot(Experiment(1).timeAPV,mean(Experiment(1).pitchAPV(154:216,:))-mean(Experiment(1).pitchAPV(1000:1050,:)),'*','Color','r')
plot(Experiment(1).timeAPVwn,mean(Experiment(1).pitchAPVwn(154:216,:))-mean(Experiment(1).pitchAPVwn(1000:1050,:)),'*','Color','g')
plot(Experiment(1).timeACpost,mean(Experiment(1).pitchACpost(154:216,:))-mean(Experiment(1).pitchACpost(1000:1050,:)),'*')
hold on;plot(Experiment(2).timeACpre,mean(Experiment(2).pitchACpre(170:234,:))-mean(Experiment(2).pitchACpre(1000:1050,:)),'*')
plot(Experiment(2).timeAPV,mean(Experiment(2).pitchAPV(170:234,:))-mean(Experiment(2).pitchAPV(1000:1050,:)),'*','Color','r')
plot(Experiment(2).timeAPVwn,mean(Experiment(2).pitchAPVwn(170:234,:))-mean(Experiment(2).pitchAPVwn(1000:1050,:)),'*','Color','g')
plot(Experiment(2).timeACpost,mean(Experiment(2).pitchACpost(170:234,:))-mean(Experiment(2).pitchACpost(1000:1050,:)),'*')

% pu56 data is just plain weird - damage likely
figure;hold on;plot(Experiment(3).timeACpre,mean(Experiment(3).pitchACpre(190:200,:))-mean(Experiment(3).pitchACpre(1000:1050,:)),'*')
plot(Experiment(3).timeAPV,mean(Experiment(3).pitchAPV(190:200,:))-mean(Experiment(3).pitchAPV(1000:1050,:)),'*','Color','r')
plot(Experiment(3).timeAPVwn,mean(Experiment(3).pitchAPVwn(190:200,:))-mean(Experiment(3).pitchAPVwn(1000:1050,:)),'*','Color','g')
plot(Experiment(3).timeACpost,mean(Experiment(3).pitchACpost(190:200,:))-mean(Experiment(3).pitchACpost(1000:1050,:)),'*')
figure;hold on;plot(Experiment(4).timeACpre,mean(Experiment(4).pitchACpre(259:323,:))-mean(Experiment(4).pitchACpre(1000:1050,:)),'*')
plot(Experiment(4).timeAPV,mean(Experiment(4).pitchAPV(259:323,:))-mean(Experiment(4).pitchAPV(1000:1050,:)),'*','Color','r')
plot(Experiment(4).timeAPVwn,mean(Experiment(4).pitchAPVwn(259:323,:))-mean(Experiment(4).pitchAPVwn(1000:1050,:)),'*','Color','g')
plot(Experiment(4).timeACpost,mean(Experiment(4).pitchACpost(259:323,1:221))-mean(Experiment(4).pitchACpost(1000:1050,1:221)),'*')
figure;hold on;plot(Experiment(5).timeACpre,mean(Experiment(5).pitchACpre(231:295,1:221))-mean(Experiment(5).pitchACpre(950:1000,1:221)),'*')
plot(Experiment(5).timeAPV,mean(Experiment(5).pitchAPV(231:295,:))-mean(Experiment(5).pitchAPV(950:1000,:)),'*','Color','r')
plot(Experiment(5).timeAPVwn,mean(Experiment(5).pitchAPVwn(231:295,:))-mean(Experiment(5).pitchAPVwn(950:1000,:)),'*','Color','g')
plot(Experiment(5).timeACpost,mean(Experiment(5).pitchACpost(231:295,:))-mean(Experiment(5).pitchACpost(950:1000,:)),'*')
figure;hold on;plot(Experiment(6).timeACpre,mean(Experiment(6).pitchACpre(222:286,:))-mean(Experiment(6).pitchACpre(950:1000,:)),'*')
plot(Experiment(6).timeAPV,mean(Experiment(6).pitchAPV(222:286,:))-mean(Experiment(6).pitchAPV(950:1000,:)),'*','Color','r')
plot(Experiment(6).timeAPVwn,mean(Experiment(6).pitchAPVwn(222:286,:))-mean(Experiment(6).pitchAPVwn(950:1000,:)),'*','Color','g')
plot(Experiment(6).timeACpost,mean(Experiment(6).pitchACpost(222:286,:))-mean(Experiment(6).pitchACpost(950:1000,:)),'*')

% b10o7
figure;hold on;plot(Experiment(7).timeACpre,mean(Experiment(7).pitchACpre(150:170,:))-mean(Experiment(7).pitchACpre(1580,:)),'*')
plot(Experiment(7).timeAPV,mean(Experiment(7).pitchAPV(150:170,:))-mean(Experiment(7).pitchAPV(1580,:)),'*','Color','r')
plot(Experiment(7).timeAPVwn,mean(Experiment(7).pitchAPVwn(150:170,:))-mean(Experiment(7).pitchAPVwn(1580,:)),'*','Color','g')
plot(Experiment(7).timeACpost,mean(Experiment(7).pitchACpost(150:170,:))-mean(Experiment(7).pitchACpost(1580,:)),'*')
figure;hold on;plot(Experiment(8).timeACpre,mean(Experiment(8).pitchACpre(150:170,:))-mean(Experiment(8).pitchACpre(1580,:)),'*')
plot(Experiment(8).timeAPV,mean(Experiment(8).pitchAPV(150:170,:))-mean(Experiment(8).pitchAPV(1580,:)),'*','Color','r')
plot(Experiment(8).timeAPVwn,mean(Experiment(8).pitchAPVwn(150:170,:))-mean(Experiment(8).pitchAPVwn(1580,:)),'*','Color','g')
plot(Experiment(8).timeACpost,mean(Experiment(8).pitchACpost(150:170,:))-mean(Experiment(8).pitchACpost(1580,:)),'*')

% two notes later 
figure;hold on;plot(E7timepre,mean(pitchE7pretargeted(150:170,:))-mean(pitchE7pre(300:350,:)),'*')
hold on;plot(Experiment(7).timeAPV(1:51),mean(Experiment(7).pitchAPV(150:170,1:51))-mean(pitchE7apv(300:350,:)),'*','Color','r')
plot(Experiment(7).timeAPVwn(1:110),mean(Experiment(7).pitchAPVwn(150:170,1:110))-mean(pitchE7apvwn(300:350,:)),'*','Color','g')
plot(Experiment(7).timeACpost,mean(Experiment(7).pitchACpost(150:170,:))-mean(pitchE7post(300:350,:)),'*')

hold on;plot(Experiment(7).timeACpost,mean(Experiment(7).pitchACpost(150:170,:))-mean(pitchE8pre(300:350,:)),'*')
hold on;plot(Experiment(8).timeAPV,mean(Experiment(8).pitchAPV(150:170,:))-mean(pitchE8apv(300:350,:)),'*','Color','r')
plot(Experiment(8).timeAPVwn,mean(Experiment(8).pitchAPVwn(150:170,:))-mean(pitchE8apvwn(300:350,:)),'*','Color','g')
plot(Experiment(8).timeACpost,mean(Experiment(8).pitchACpost(150:170,:))-mean(pitchE8post(300:350,:)),'*')


% next note
figure;hold on;plot(Experiment(10).timeACpre,mean(Experiment(10).pitchACpre(150:204,:))-mean(Experiment(10).pitchACpre(850:1100,:)),'*')
plot(Experiment(10).timeAPV,mean(Experiment(10).pitchAPV(150:204,:))-mean(Experiment(10).pitchAPV(850:1100,:)),'*','Color','r')
plot(Experiment(10).timeAPVwn,mean(Experiment(10).pitchAPVwn(150:204,:))-mean(Experiment(10).pitchAPVwn(850:1100,:)),'*','Color','g')
plot(Experiment(10).timeACpost,mean(Experiment(10).pitchACpost(150:204,:))-mean(Experiment(10).pitchACpost(850:1100,:)),'*')
hold on;plot(Experiment(11).timeACpre,mean(Experiment(11).pitchACpre(150:204,:))-mean(Experiment(11).pitchACpre(850:1100,:)),'*')
plot(Experiment(11).timeAPV,mean(Experiment(11).pitchAPV(150:204,:))-mean(Experiment(11).pitchAPV(850:1100,:)),'*','Color','r')
plot(Experiment(11).timeAPVwn,mean(Experiment(11).pitchAPVwn(150:204,:))-mean(Experiment(11).pitchAPVwn(850:1100,:)),'*','Color','g')
plot(Experiment(11).timeACpost,mean(Experiment(11).pitchACpost(150:204,:))-mean(Experiment(11).pitchACpost(850:1100,:)),'*')


% bk75
figure;plot(tvals927(31:end),mean(pitch927(320:380,31:end))-pitch927CCA(350,31:end),'*')
hold on;plot(tvals927apv,mean(pitch927apv(320:380,:))-pitch927apvCCA(350,:),'*','Color','r')
hold on;plot(tvals927apvwn,mean(pitch927apvwn(320:380,:))-pitch927apvwnCCA(350,:),'*','Color','g')
hold on;plot(tvals928acsf,mean(pitch928acsf(320:380,:))-pitch928acsfCCA(350,:),'*','Color','b')
