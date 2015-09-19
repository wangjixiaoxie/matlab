% go into a data directory
% all commands are from matlab


%STEP 1: FILE
dirf('*.cbin','batch')
findcatch('batch')
%STEP 2: rename batch file and label notes -- 'a' is default
evsonganaly
%STEP 2: summary_stats2.m
[fvals1205A,FFstats1205A,avg1205A,toffs1205A,pitch1205A]=summary_statsZF(templa1,'batch1205Afiles','batch1205Anotes',cntrng1);
[fvals1205B,FFstats1205B,avg1205B,toffs1205B,pitch1205B]=summary_statsZF(templa2,'batch1205Bfiles','batch1205Bnotes',cntrng2);
[fvals1205C,FFstats1205C,avg1205C,toffs1205C,pitch1205C]=summary_statsZF(templa2,'batch1205Cfiles','batch1205Cnotes',cntrng2);

[fvals1120,FFstats1120,avg1120,toffs1120,pitch1120]=summary_statsJC(templa1,'batchJCfiles','batchJCnotes',cntrng1);

%bk50
[fvalsT,FFstatsT,avgT,toffsT,pitchT]=summary_statsJC(templa1,'batchJCfiles','batchJCnotes',cntrng1);
%STEP 3: summary_plot.m
[xax1]=summary_plot(FFstatsT3,'batchT',avgT3);
%STEP 4: plot a baseline curve
hold on;plot(xax1,a913,'k')


%BK50w18-before - just taking "catch" trials
[fvals1,FFstats1,avg1,toffset1,pitch1]=summary_stats5(templa,'batch927Afiles','batch927Anotes',cntrng);
[xax1]=summary_plot2(FFstats1,'batchtest',avg1);

%STEP 5: look at the actual pitch curves
figure;hold on;
for i=1:size(pitch1,2)
    plot(pitch1(:,i)+250*i)
end

%%%TSAnaly-- note templa1/cntrng1 vs. templa/cntrng
fvals917=findwnote4('batchTSAnotes','a','','',0,[2000 2700],8500,1,'obs0',0);
mk_tempf('batchTSAfiles',templa1,2,'obs0');
get_trigt2('batchTSAfiles',cntrng1,0.3,128,0,1);
for i=1:length(fvals)
    shifted(i,:)=fvals(i).datt;
end
[pitchCTLtsa]=jc_pitchmat1204(shifted,1204,1020,1,1950,2600,[1 ],'obs0',1);
for i=0:20
[mm919,postesc919(i+1),posthit919(i+1)]=jcTSA919('batch919TSAfiles',pitch919tsa,fvals919,i);
end

%Simulation
[vals,trigs]=triglabel('batch918Bfiles','a',1,1,0,1);
toff=[];
for ii=1:length(trigs)
toff=[toff;trigs(ii).toffset];
end
toffset=((toff/1000)*(32000)-512)/4;
aavfin=ContingSim(fvals,toffset,2500); %where 2500 is the contingency to go above


%%%%WAV FILES from other room%%%%%%
fvals1001=findwnote4('batchtestnotes','a','','',0,[2000 2700],8500,1,'w');
for i=1:length(fvals1001)
    shifted1001(i,:)=fvals1001(i).datt;
end
[pitch1001,avg1001]=jc_pitchmat1204(shifted1001,1204,1020,1,1950,2600,[1 ],'w',1);
avg1001res=resample(avg1001,320,441);
figure;hold on;
plot(avg1001res),plot(avgBaseline,'r')


%%%%Plot the stuff as standard deviations%%%%
%1. Look at get_avn to determine note onset and offset times
%Actually the best way is to look at where the pitch curve levels off
B.notonset=150;
B.notoffset=800;
B.notlength=B.notoffset-B.notonset+1;
%2. Choose the shifted data and the baseline data
B.avgBaseline=avgbk1+avgbk2;
B.avgBaseline=B.avgBaseline./2;
B.pitchBaseline=[pitchbk1,pitchbk2];
B.stdBaseline=std(B.pitchBaseline');
B.avgShifted=avgbk7;
%avgShifted=avgShifted./6;
%Calculate and plot the trace
B.SDtrace=(B.avgShifted-B.avgBaseline)./B.stdBaseline;
B.SDTplot=B.SDtrace(B.notonset:B.notoffset);
for i=1:B.notlength
    B.xax(i)=i/8;
end
figure;plot(B.xax,B.SDTplot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(fvals913a)
    shifted(i,:)=fvals913a(i).datt;
end
[pitch913a,avg913a]=jc_pitchmat1204(shifted,1204,1020,1,1950,2600,[1 2 3],'obs0',1);
for i=1:length(fvals913b)
    shifted(i,:)=fvals913b(i).datt;
end
[pitch913b,avg913b]=jc_pitchmat1204(shifted,1204,1020,1,1950,2600,[1 2 3],'obs0',1);
for i=1:length(fvals913c)
    shifted(i,:)=fvals913c(i).datt;
end
[pitch913c,avg913c]=jc_pitchmat1204(shifted,1204,1020,1,1950,2600,[1 2 3],'obs0',1);
for i=1:length(fvals914a)
    shifted(i,:)=fvals914a(i).datt;
end
[pitch914a,avg914a]=jc_pitchmat1204(shifted,1204,1020,1,1950,2600,[1 2 3],'obs0',1);
for i=1:length(fvals914c)
    shifted(i,:)=fvals914c(i).datt;
end
[pitch914c,avg914c]=jc_pitchmat1204(shifted,1204,1020,1,1950,2600,[1 2 3],'obs0',1);
for i=1:length(fvals915a)
    shifted(i,:)=fvals915a(i).datt;
end
[pitch915a,avg915a]=jc_pitchmat1204(shifted,1204,1020,1,1950,2600,[1 2 3],'obs0',1);
for i=1:length(fvals915b)
    shifted(i,:)=fvals915b(i).datt;
end
[pitch915b,avg915b]=jc_pitchmat1204(shifted,1204,1020,1,1950,2600,[1 2 3],'obs0',1);
for i=1:length(fvals915c)
    shifted(i,:)=fvals915c(i).datt;
end
[pitch915c,avg915c]=jc_pitchmat1204(shifted,1204,1020,1,1950,2600,[1 2 3],'obs0',1);

