%D:\jcharles\Circadian
%As of 10/22/08 (Wednesday before trip to Princeton)
%BK20summary_final.mat
%BK50summary_final.mat

BK20
'b' - 150:250
%The general idea is to see if the change in pitch overnight is
%anticorrelated with the learning on the subsequent day.
%To measure the change in pitch overnight, we take the last ~50 notes sung 
%at dusk and the first ~50 notes sung at dawn.  We measure pitch using the
%triggers (method 1) and also by taking the median pitch curve over the
%interval after the note has begun but before the wn blast -- this
%interval is chosen based on low standard deviation (method 2).

%To measure learning, we measure the difference between song sung between
%4pm and 9pm on one day vs. song sung between 4pm and 9pm on the next day.
%We choose late in the day to eliminate any effect of morning reversion
%that is our independent variable.  We use all notes sung during catch
%trials during this period (catch trials are 10% of song).

%Measuring change in pitch overnight
    %select only last 20 songs of evening/first 20 of morning
    cleandir4('batch',1000,500,6,10);
    %label all
    evsonganaly
    %These are called batch13dusk/batch14dawn
    %   1. Triggers method --- change 2000 and 3000 depending on the note
    dusk13B=evtaf_freq('batch13dusk',[2000 3000],'b',128,'obs0',0,0);
    dusk13A=evtaf_freq('batch713dusk',[3000 4500],'a',128,'obs0',0,0);
    dusk13B=dusk13B(:,2);
    DuskBff(1)=median(dusk13b);
    %   2. Pitch contour method ---2200 to 2800 for bk20 'b'; 
            %3100 to 3900 for bk20 'a'
    fvals13PMb=findwnoteJC('batch17dusknotes','b','','',0,[2000 3000],6000,1,'obs0')
    for i=1:50
        shifted13PMb(i,:)=fvals13PMb(i).datt;
    end
    [pitch13PMb,avg13PMb]=jc_pitchmat1024(shifted13PMb,1024,1020,1,2200,2800,[1 2 3],'obs0',1);
            %Determine where to draw limits 
                %130:230 for bk20 'a'; 150:250 for bk20 'b'; 
                %600:700 for bk50 'a'
    figure;plot(std(pitch13PMb'))
            %
    medP13PMb=median(pitch13PMb');
    DuskB(1)=median(medP13PMb(150:250));
%Measuring learning during the day
    %choose song between 4pm and 9pm (*Cfiles)
    %make folder of catch trials in this time period
    fv17a=findwnoteJC2('batch17Cnotes','a','-','a',0.01,[3000 4500],3000,1,'obs0');
    vvals17a=getvals(fv17a,1,'TRIG');
    FF17pmA=vvals17a(:,2);
    Night2NightLearning(4)=median(FF17pmA)-median(FF16pmA);
%Compare
    OvernightA_contours=DawnA-DuskA;
    OvernightA_triggers=DawnAff-DuskAff;
    %plot these with Night2NightLearning
    %OvernightA_contours and OvernightA_triggers should generally agree
    %If there were a relationship, you would expect to see the largest positive values of 
    %Night2NightLearning at the same position as the largest negative values
    %of OvernightA in a bird that you are driving up.  This would mean that
    %learning from evening 1 to evening 2 correlates with a deterioration
    %towards baseline during night 1.  %%%Note it might be a good idea to
    %show baseline to indicate that it is indeed a deterioration towards
    %baseline.
    
    
    
    
    
    
    
    
    
    

%make batch of all files
BK50
%'a' - 3100:3900
%'b' - 2200:2800
%select only last 20 songs of evening/first 20 of morning
cleandir4('batch',1000,500,6,10);
%label all
evsonganaly
%
%
fvals06PMb=findwnoteJC('batch16dusknotes','b','','',0,[3000 4500],6000,1,'obs0')
fvals17AMb=findwnoteJC('batch17dawnnotes','b','','',0,[3000 4500],6000,1,'obs0')

[pitch20AMa,avg20AMa]=jc_pitchmat1024(shifted20AMa,1024,1020,1,3100,3900,[1 2 3],'obs0',1);
figure;plot(median(pitch16PMa(100:260,:)'))
DawnB(3)=mean(median(pitch16AMb(160:240)));

vals1008dawnA=evtaf_freq('batch.catch',[3100 3900],'a',128,'obs0',0,0);
FF1008dawnA=vals1008dawnA(:,2);
DawnAff(1)=median(FF1008dawnA)

%Other notes: 
%findwnoteJC2 allows you to search for '-a' when you normally can't
fv=findwnoteJC2('batch19Cnotes','a','-a','a',0.01,[3000 4500],3000,1,'obs0');
%Get all the frequency values for the entire experiment in Tim's data
    %Use this to see how far bird has gone down and compare to
    %inactivations
ffs=avls.vls{1,1}(:,2); 