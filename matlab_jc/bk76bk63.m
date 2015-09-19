% bk76bk63
% song is eligible for both pitch and syntax shifting
% low and high stack notes
% 10.7.09 - started screening

% 10.14.09 - baseline pitch shift
% WN on at 2:05pm launchpad time - WN hits first long stack if below 2505Hz
% WN off at 5:50pm

% 11.5.09 - bilateral RA cannulae implantation 
%          - targeting uncertain because RA recording was inconclusive
% 11.17.09 (tuesday) - probes in
% 11.19.09 - began singing in the morning
% 11.20.09 - 4mMapv on at 12:31pm chick time at 1.0uL/min
%          - acsf on at 1.5uL/min at around 3:20 chick time
%       oddly enough, apv increases the length of bouts, this increasing
%       the relative proportion of B's to A's, since A's are always sung
%       first in a bout
% 11.23.09 - 4mM apv on at 2:10pm chick time at 1.0uL/min
% clog at 11.24.09 morning - probes out around noon

% 11.25.09 - new probes in at 1:30pm
% 11.26.09 - singing - but filter was high pass so lost data
% 11.26.09 at 2:44pm - corrected the filter and began recording song
% 11.27.09 at 10:15am - 4mM apv on at 1.0uL/min
%          at 11:42am - WN on hit below
% 11.30.09 at 12:15pm - hit below 2510Hz, acsf on
% 11.30.09 at 5:23pm - amp off
% 12.01.09 at 11:08pm - 4mM apv on at 1.0uL/min
%          at 1:20pm - amp on hit below 2525Hz
%          at 9pm (lights on) - amp off
% 12.03.09 at 11:38pm - 4mM apv on at 1.0uL/min
%          at 1:48pm - hit below 3400Hz (0.8uL/min)
%          at 8:07pm - ampoff/lights off/ ACSF on at 1.5uL/min (switched
%                           to 1.0uL/min at 8:15pm)
% 12.06.09 at 12:17pm - 4mM apv on at 1.0uL/min
%          at 2:10pm - wn on - hit above 3350Hz
%          at 5:09pm - probes out
%          at 6pm - lights off (no hits before lights off)
% 12.07.09 morning - singing again
% 12.14.09 - new probes in at 4pm

fvals1119B2=findwnoteJC('batchJCnotes','b','','',0,[2000 2700],8500,1,'obs0',1);
for i=1:length(fvals1119B2)
    shifted1119B2(i,:)=fvals1119B2(i).datt;
end
pitch1119B2=jc_pitchmat1024(shifted1119B2,1024,1020,1,2000,3000,[1],'obs0',1);
fvals1119B=findwnoteJC('batchJCnotes','b','','',0,[2000 2700],8500,1,'obs0',1);
for i=1:length(fvals1119B)
    shifted1119B(i,:)=fvals1119B(i).datt;
end
pitch1119B=jc_pitchmat1024(shifted1119B,1024,1020,1,2000,3000,[1],'obs0',1);

fvals1120apvA=findwnoteJC('batchJCnotes','a','','',0,[2000 2700],8500,1,'obs0',1);
for i=1:length(fvals1120apvA)
    shifted1120apvA(i,:)=fvals1120apvA(i).datt;
end
pitch1120apvA=jc_pitchmat1024(shifted1120apvA,1024,1020,1,2800,3800,[1],'obs0',1);

figure;hold on;
plot(runningaverage(tvals1123B,15),runningaverage(pitch1123B(250,:),15))
plot(runningaverage(tvals1123acpost,15),runningaverage(pitch1123acpost(250,:),15))
plot(runningaverage(tvals1123apvB,15),runningaverage(pitch1123apvB(250,:),15),'r')
plot(runningaverage(tvals1123acpost,15),runningaverage(pitch1123acpost(250,:),15))
plot(runningaverage(tvals1124B,15),runningaverage(pitch1124B(250,:),15))
plot(runningaverage(tvals1124amB,15),runningaverage(pitch1124amB(250,:),15))
plot(runningaverage(tvals1127Pre,15),runningaverage(pitch1127Pre(250,:),15))
figure;hold on;
plot(runningaverage(tvals1130pre,15),runningaverage(pitch1130pre(240,:),15))
plot(runningaverage(tvals1130wnon,15),runningaverage(pitch1130wnon(240,:),15))
plot(runningaverage(tvals1130wnoff,15),runningaverage(pitch1130wnoff(240,:),15),'Color','k')
plot(runningaverage(tvals1201pre,15),runningaverage(pitch1201pre(240,:),15))
plot(runningaverage(tvals1201apv,15),runningaverage(pitch1201apv(240,:),15),'Color','g')
plot(runningaverage(tvals1201apvwn,15),runningaverage(pitch1201apvwn(240,:),15),'Color','r')
plot(runningaverage(tvals1201post,15),runningaverage(pitch1201post(240,:),15))
plot(runningaverage(tvals1203preB,15),runningaverage(pitch1203preB(240,:),15))
plot(runningaverage(tvals1204apvB,15),runningaverage(pitch1204apvB(240,:),15),'g')
plot(runningaverage(tvals1204acsfB,15),runningaverage(pitch1204acsfB(240,:),15))
plot(runningaverage(tvals1206apvB,15),runningaverage(pitch1206apvB(240,:),15),'g')
plot(runningaverage(tvals1206postB,15),runningaverage(pitch1206postB(240,:),15))
figure;hold on;
plot(runningaverage(tvals1130wnoffC,15),runningaverage(pitch1130wnoffC(260,:),15))
plot(runningaverage(tvals1201postC(150:end),15),runningaverage(pitch1201postC(260,150:end),15))
plot(runningaverage(tvals1201apvC,15),runningaverage(pitch1201apvC(260,:),15),'g')
plot(runningaverage(tvals1203preC,15),runningaverage(pitch1203preC(260,:),15))
plot(runningaverage(tvals1203apvC,15),runningaverage(pitch1203apvC(260,:),15),'g')
plot(runningaverage(tvals1203apvwnC,15),runningaverage(pitch1203apvwnC(260,:),15),'r')
plot(runningaverage(tvals1204acsfC,30),runningaverage(pitch1204acsfC(260,:),30),'k')
plot(runningaverage(tvals1206apvC,15),runningaverage(pitch1206apvC(260,:),15),'g')
plot(runningaverage(tvals1206apvwnC,15),runningaverage(pitch1206apvwnC(260,:),15),'r')
plot(runningaverage(tvals1206postC,30),runningaverage(pitch1206postC(260,:),30))

% removing the early morning stuff - which distracts
figure;hold on;
plot(runningaverage(tvals1130wnoffC,15),runningaverage(pitch1130wnoffC(260,:),15))
plot(runningaverage(tvals1201postC(150:end),15),runningaverage(pitch1201postC(260,150:end),15))
plot(runningaverage(tvals1201apvC,15),runningaverage(pitch1201apvC(260,:),15),'g')
plot(runningaverage(tvals1203preC,15),runningaverage(pitch1203preC(260,:),15))
plot(runningaverage(tvals1203apvC,15),runningaverage(pitch1203apvC(260,:),15),'g')
plot(runningaverage(tvals1203apvwnC,15),runningaverage(pitch1203apvwnC(260,:),15),'r')
plot(runningaverage(tvals1204acsfC(1:616),30),runningaverage(pitch1204acsfC(260,1:616),30),'k')
plot(runningaverage(tvals1204acsfC(662:end),30),runningaverage(pitch1204acsfC(260,662:end),30),'k')
plot(runningaverage(tvals1206apvC,15),runningaverage(pitch1206apvC(260,:),15),'g')
plot(runningaverage(tvals1206apvwnC,15),runningaverage(pitch1206apvwnC(260,:),15),'r')
plot(runningaverage(tvals1206postC(49:393),30),runningaverage(pitch1206postC(260,49:393),30))
plot(runningaverage(tvals1206postC(481:end),30),runningaverage(pitch1206postC(260,481:end),30))