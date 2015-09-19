% bk80w28
% Did something (hit/miss) with him in August 2009 - see /cardinal6
% September 2010
% 09.03 - put in launchpad soundbox - for NatNeuro ctls
%   wntest from 12-2pm - never sang through but good initiation of song

% 09.07 - wn test for notch - too loud
%   template test - first threshold crossing
%   hit below 2285Hz at 8pm
% 09.08 
%   hit below 2300Hz at 10:06am
%   hit below 2320Hz at 12:55pm
%   hit below 2335Hz at 3:25pm
%   hit below 2350Hz at 4:20pm
%   hit below 2365Hz at 5pm
%   hit below 2380Hz at 8:15pm
%  looks like he learned a lot in the morning and then it has leveled off

% 09.09
%   hit below 2420Hz
% 09.10
%   WN off at 12:25pm - allow recovery and then repeat the same experiment

% 09.12 - median for entire day is 2300Hz
%   10:30pm - wn on, hit below 2300Hz
% 09.13
%   10:38am - hit below 2330Hz
%   10:46am - hit below 2380Hz
%   11:05am - adjusted template slightly to recognize BTAF better
%   12:23pm - hit below 2400Hz


fvBaseline=findwnoteJC('batchnotes','b','','',0,[6000 8100],8500,1,'obs0',0);
for i=1:length(fvBaseline)
    shiftedBaseline(i,:)=fvBaseline(i).datt;
end
pitchBaseline=jc_pitchmat1024(shiftedBaseline,1024,1020,1,6000,8100,[1],'obs0',1);
fvWNon=findwnoteNN('batchnotes','b','','',0,[6000 8100],8500,1,'obs0',0);
        recognized=[];
        hitrange=[];
        catchtrial=[];
        for i=1:length(fvWNon)
            shiftedWNon(i,:)=fvWNon(i).datt;
            recognized(i)=fvWNon(i).CATCH>-1;
            hitrange(i)=fvWNon(i).TEMPLATE==1;
            catchtrial(i)=fvWNon(i).CATCH==1;
        end
        indEscapeAbove=recognized & hitrange & catchtrial;
        indHitAbove=recognized & hitrange & ~catchtrial;
pitchWNon=jc_pitchmat1024(shiftedWNon,1024,1020,1,6000,8100,[1],'obs0',1);    
figure;plot(median(pitchWNon(:,indEscapeAbove)'))
hold on;plot(median(pitchWNon(:,indHitAbove)'),'r')

for i=1:size(pitchWNon,2)
    residWNon(:,i)=pitchWNon(:,i)-mean(pitchBaseline')';
end
b=zeros(1363,1363);
for i=1:1363
    [b(i,1:i),bint(i).data]=regress(residWNon(:,i),residWNon(:,1:i));
    i
end    
regcoefs=b;
clear mb
for i=1:length(regcoefs)
    mb(i)=mean(regcoefs(i+1:length(regcoefs),i));
end



% load /cardinal3/bk80w28/Analyzed910.mat
    % Beautiful!!!!!
        figure;plot(runningaverage(t(indEscapeAbove),50),runningaverage(mb(indEscapeAbove),50),'r')
        hold on;plot(runningaverage(t(indHitAbove),50),runningaverage(mb(indHitAbove),50),'b')












fvals80w7Ahit=findwnoteJC('batch807Anotes','a','','',0,[2000 2700],8500,1,'obs0');
fvals807Amiss=findwnoteJC('batch807Anotes','b','','',0,[2000 2700],8500,1,'obs0');
for i=1:length(fvals807Ahit)
shifted807Ahit(i,:)=fvals807Ahit(i).datt;
end
for i=1:length(fvals807Amiss)
shifted807Amiss(i,:)=fvals807Amiss(i).datt;
end
pitch807Ahit=jc_pitchmat1024(shifted807Ahit,1024,1020,1,2000,2700,[1],'obs0',1);
pitch807Amiss=jc_pitchmat1024(shifted807Amiss,1024,1020,1,2000,2700,[1],'obs0',1);
pitch807A=[pitch807Ahit,pitch807Amiss];
figure;plot(mean(pitchBase'))
hold on;plot(mean(pitch807A'),'r')
hold on;plot(mean(pitch807Amiss'),'g')
hold on;plot(mean(pitch807Ahit'),'k')
