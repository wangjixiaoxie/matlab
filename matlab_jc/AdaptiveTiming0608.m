% AdaptiveTiming

% Consider the following analyses

    % Compare curves when ACSF(non-targeted time) = MU(non-targeted time)
         % Localized timing is the last thing to get consolidated in the DMP
         % AFP is necessary to encode precise timing information
         % pk16r47 is a great example of this
    % Compare curves when ACSF(entire note) = MU(entire note)
        
% Analysis of Muscimol in LMAN
% bk61w42       /cardinal
        figure;hold on;
        plot(tvACpre,pitchACpre(420,:),'*','Color','k')
        plot(tvACwn1,pitchACwn1(420,:),'*')
        plot(tvACwn2,pitchACwn2(420,:),'*')
        plot(tvACwn3,pitchACwn3(420,:),'*')
        plot(tvACwn4,pitchACwn4(420,:),'*')
        plot(tvMUwn1,pitchMUwn1(420,:),'*','Color','r') % 1207
        plot(tvMUwn2,pitchMUwn2(420,:),'*','Color','r') % 1209
        % plot(tvMUwn3,pitchACwn3(420,:),'*','Color','r') % 1211 -didn't sing
        plot(tvMUwn4,pitchMUwn4(420,:),'*','Color','r') % 1212
        plot(tvMUpre3,pitchMUpre3(420,:),'*','Color','r') % 1204
        plot(tvMUpre1,pitchMUpre1(420,:),'*','Color','r') % 1203
        plot(tvMUpre2,pitchMUpre2(420,:),'*','Color','r') % 1202

        pitchMUpre=[pitchMUpre2 pitchMUpre1 pitchMUpre3];


        figure;hold on;
        plot(mean(pitchMUwn1')-mean(pitchMUpre'),'r')
        plot(mean(pitchMUwn4')-mean(pitchMUpre'),'r')
        plot(mean(pitchMUwn2')-mean(pitchMUpre'),'r')
        plot(mean(pitchACwn1(:,140:200)')-mean(pitchACpre'),'b')
        plot(mean(pitchACwn1(:,200:300)')-mean(pitchACpre'),'b')
        plot(mean(pitchACwn2(:,1:100)')-mean(pitchACpre'),'b')
        xlim([400 800])
figure;subplot(121);hold on;
yvalsAC=(runningaverage([mean(pitchACwn1(400:450,:)) mean(pitchACwn2(400:450,:)) mean(pitchACwn3(400:450,:))],20));
xvalsAC=(runningaverage([mean(pitchACwn1(400:800,:)) mean(pitchACwn2(400:800,:)) mean(pitchACwn3(400:800,:))],20));
plot(xvalsAC,yvalsAC,'*')    
yvalsMU=(runningaverage([mean(pitchMUwn1(400:450,:)) mean(pitchMUwn2(400:450,:))  mean(pitchMUwn4(400:450,:))],20));
xvalsMU=(runningaverage([mean(pitchMUwn1(400:800,:)) mean(pitchMUwn2(400:800,:))  mean(pitchMUwn4(400:800,:))],20));
plot(xvalsMU,yvalsMU,'*','Color','r')        
subplot(122);hold on;
mnbaseACy=mean(mean(pitchACpre(400:450,:)));
mnbaseACx=mean(mean(pitchACpre(400:800,:)));
mnbaseMUy=mean(mean(pitchMUpre(400:450,:)));
mnbaseMUx=mean(mean(pitchMUpre(400:800,:)));
yvalsAC=(runningaverage([mean(pitchACwn1(400:450,:)) mean(pitchACwn2(400:450,:)) mean(pitchACwn3(400:450,:))],20));
xvalsAC=(runningaverage([mean(pitchACwn1(400:800,:)) mean(pitchACwn2(400:800,:)) mean(pitchACwn3(400:800,:))],20));
plot(xvalsAC-mnbaseACx,yvalsAC-mnbaseACy,'*')    
yvalsMU=(runningaverage([mean(pitchMUwn1(400:450,:)) mean(pitchMUwn2(400:450,:))  mean(pitchMUwn4(400:450,:))],20));
xvalsMU=(runningaverage([mean(pitchMUwn1(400:800,:)) mean(pitchMUwn2(400:800,:))  mean(pitchMUwn4(400:800,:))],20));
plot(xvalsMU-mnbaseMUx,yvalsMU-mnbaseMUy,'*','Color','r')        
plot([0 500],[0 500],'k')


% bk63w43       /cardinal
 % first experiment - December 10 

 j=3
 window=[round(median(Data(j).targs-20)):round(median(Data(j).targs))];
 figure;hold on;
 for i=1:length(Data(j).exp)
     if (Data(j).baseline(i)*Data(j).acsf(i))
         plot(Data(j).exp(i).time,mean(Data(j).exp(i).pitch(window,:)),'*','Color','k')
     else
     if 1-Data(j).acsf(i)
         plot(Data(j).exp(i).time,mean(Data(j).exp(i).pitch(window,:)),'*','Color','r')
     else
         plot(Data(j).exp(i).time,mean(Data(j).exp(i).pitch(window,:)),'*','Color','b')  
     end
     end
 end
 
    figure;hold on;
    plot(tvACpre,pitchACpre(480,:),'*','Color','k')
    plot(tvACwn1,pitchACwn1(480,:),'*')
    plot(tvMUpre1,pitchMUpre1(480,:),'*','Color','r')
    plot(tvMUpre2,pitchMUpre2(480,:),'*','Color','r')
    plot(tvMUwn1,pitchMUwn1(480,:)/3,'*','Color','r')
    plot(tvMUwn2,pitchMUwn2(480,:)/3,'*','Color','r')    
    plot(tvMUwn3,pitchMUwn3(480,:)/3,'*','Color','r')    
 figure;hold on;
 plot(median(pitchMUwn3')/3-median([pitchMUpre1 pitchMUpre2]'),'r')
 plot(median(pitchACwn1(:,240:300)')-median(pitchACpre'))
 xlim([370 760])
 
 
 
figure;subplot(121);hold on;
yvalsAC=(runningaverage([mean(pitchACwn1(440:490,:))],20));
xvalsAC=(runningaverage([mean(pitchACwn1(370:750,:))],20));
plot(xvalsAC,yvalsAC,'*')    
yvalsMU=(runningaverage([mean(pitchMUwn1(440:490,:)) mean(pitchMUwn2(440:490,:)) mean(pitchMUwn3(440:490,:)) mean(pitchMUwn4(440:490,:)) mean(pitchMUwn5(440:490,:))]/3,20));
xvalsMU=(runningaverage([mean(pitchMUwn1(370:750,:)) mean(pitchMUwn2(370:750,:)) mean(pitchMUwn3(370:750,:)) mean(pitchMUwn4(370:750,:)) mean(pitchMUwn5(370:750,:))]/3,20));
plot(xvalsMU,yvalsMU,'*','Color','r')        
subplot(122);hold on;
mnbaseACy=mean(mean(pitchACpre(440:490,:)));
mnbaseACx=mean(mean(pitchACpre(370:750,:)));
mnbaseMUy=mean(mean(pitchMUpre(440:490,:)));
mnbaseMUx=mean(mean(pitchMUpre(370:750,:)));
plot(-1*(xvalsAC-mnbaseACx),-1*(yvalsAC-mnbaseACy),'*')    
plot(-1*(xvalsMU-mnbaseMUx),-1*(yvalsMU-mnbaseMUy),'*','Color','r')        
plot([0 200],[0 200],'k')



 % second experiment May 28- DMP never learned much
    figure;hold on;
    plot(tvEXP2_ACpre1,pitchEXP2_ACpre1(500,:),'*','Color','k')
    plot(tvEXP2_ACwn1,pitchEXP2_ACwn1(500,:),'*')
    plot(tvEXP2_MUwn1,pitchEXP2_MUwn1(500,:),'*','Color','r')
    plot(tvEXP2_MUpre1,pitchEXP2_MUpre1(500,:),'*','Color','r')
    plot(tvEXP2_MUpre2,pitchEXP2_MUpre2(500,:),'*','Color','r')


    figure;hold on;
    plot(mean(pitchEXP2_ACpre1'))
    plot(mean(pitchEXP2_ACwn1(:,end-70:end)'))
    plot(mean(pitchEXP2_MUwn1'),'r')
    plot(mean(pitchEXP2_MUpre2'),'r')
figure;subplot(121);hold on;
yvalsAC=(runningaverage([mean(pitchEXP2_ACwn1(460:500,:))],20));
xvalsAC=(runningaverage([mean(pitchEXP2_ACwn1(370:750,:))],20));
plot(xvalsAC,yvalsAC,'*')    
yvalsMU=(runningaverage([mean(pitchEXP2_MUwn1(460:500,:))],20));
xvalsMU=(runningaverage([mean(pitchEXP2_MUwn1(370:750,:))],20));
plot(xvalsMU,yvalsMU,'*','Color','r')        
subplot(122);hold on;
mnbaseACy=mean(mean(pitchEXP2_ACpre1(460:500,:)));
mnbaseACx=mean(mean(pitchEXP2_ACpre1(370:750,:)));
mnbaseMUy=mean(mean(pitchEXP2_MUpre2(460:500,:)));
mnbaseMUx=mean(mean(pitchEXP2_MUpre2(370:750,:)));
plot((xvalsAC-mnbaseACx),(yvalsAC-mnbaseACy),'*')    
plot((xvalsMU-mnbaseMUx),(yvalsMU-mnbaseMUy),'*','Color','r')        
plot([0 200],[0 200],'k')

% pu56w26       /cardinal5
% baseline - 21709-2, 21709_ACSF, 219_ACSF
% wn on - 222tmp2, 224tmp3, 225ACSF, 227_ACSF
        figure;hold on;
        xlim([110 330])
        plot(mean(pitchMUwn3catch(1:1001,:)')-mean(pitchMUpre1(65:1065,:)'),'r') % 0302
        plot(mean([pitchACwn2catch(65:1065,end-25:end) pitchACwn3(1:1001,1:25)]')-mean(pitchACpre3(65:1065,:)'))
        plot(mean(pitchMUwn2catch(1:1001,:)')-mean(pitchMUpre1(65:1065,:)'),'r') % 0227
        plot(mean([pitchACwn2catch(65:1065,41:90)]')-mean(pitchACpre3(65:1065,:)'))
        
        
           figure;subplot(121);hold on;
        yvalsAC=(runningaverage([mean(pitchACwn1catch(64+250:64+290,:)) mean(pitchACwn2catch(64+250:64+290,:)) ],20));
        xvalsAC=(runningaverage([mean(pitchACwn1catch(64+110:64+330,:)) mean(pitchACwn2catch(64+110:64+330,:)) ],20));
        plot(xvalsAC,yvalsAC,'*')    
        yvalsMU=(runningaverage([mean(pitchMUwn1(250+64:290+64,:))/2 mean(pitchMUwn2catch(250:290,:)) mean(pitchMUwn3catch(250:290,:)) ],20));
        xvalsMU=(runningaverage([mean(pitchMUwn1(110+64:330+64,:))/2 mean(pitchMUwn2catch(110:330,:)) mean(pitchMUwn3catch(110:330,:)) ],20));
        plot(xvalsMU,yvalsMU,'*','Color','r')        
        subplot(122);hold on;
        mnbaseACy=mean(mean(pitchACpre3(250+64:290+64,:)));
        mnbaseACx=mean(mean(pitchACpre3(110+64:330+64,:)));
        mnbaseMUy=mean(mean(pitchMUpre1(250+64:290+64,:)));
        mnbaseMUx=mean(mean(pitchMUpre1(110+64:330+64,:)));
        plot(-1*(xvalsAC-mnbaseACx),-1*(yvalsAC-mnbaseACy),'*')    
        plot(-1*(xvalsMU-mnbaseMUx),-1*(yvalsMU-mnbaseMUy),'*','Color','r')        
        plot([0 200],[0 200],'k')     
% pu57w52
        figure;plot(mean(pitchMUwn1catch(1:1001,:)')-mean(pitchMUpre1catch(65:1065,:)'),'r') % 20909_APV_5mM and 21209_APV_5mM
        hold on;plot(mean(pitchACwn1catch(65:1065,400:500)')-mean(pitchACpre1catch(65:1065,:)'),'b') % baseline is morning of 209 - notes are from wnon2909eve
        plot(mean(pitchMUwn2catch(65:1065,:)')-mean(pitchMUpre1catch(65:1065,:)'),'r') % 20909_APV_5mM and 21409_APV_5mM
        %hold on;plot(mean(pitchACwn2catch(65:1065,1:20)')-mean(pitchACpre1catch(65:1065,:)'),'b') % baseline is morning of 209 - notes are from wnon2909eve
        hold on;plot(mean([pitchACwn1catch(65:1065,end-11:end) pitchAC212catch(1:1001,:)]')-mean(pitchACpre1catch(65:1065,:)'),'b') % baseline is morning of 209 - notes are from wnon2909eve
        plot([median(targs)-64 median(targs)-64],[0 100])
        
        
           figure;subplot(121);hold on;
        yvalsAC=(runningaverage([mean(pitchACwn1catch(64+250:64+290,:)) mean(pitchACwn2catch(64+250:64+290,:)) ],20));
        xvalsAC=(runningaverage([mean(pitchACwn1catch(64+110:64+340,:)) mean(pitchACwn2catch(64+110:64+340,:)) ],20));
        plot(xvalsAC,yvalsAC,'*')    
        yvalsMU=(runningaverage([mean(pitchMUwn1catch(250:290,:)) mean(pitchMUwn2catch(64+250:64+290,:))],20));
        xvalsMU=(runningaverage([mean(pitchMUwn1catch(110:330,:)) mean(pitchMUwn2catch(64+110:64+330,:))],20));
        plot(xvalsMU,yvalsMU,'*','Color','r')        
        subplot(122);hold on;
        mnbaseACy=mean(mean(pitchACpre1catch(250+64:290+64,:)));
        mnbaseACx=mean(mean(pitchACpre1catch(110+64:330+64,:)));
        mnbaseMUy=mean(mean(pitchMUpre1catch(250+64:290+64,:)));
        mnbaseMUx=mean(mean(pitchMUpre1catch(110+64:330+64,:)));
        plot((xvalsAC-mnbaseACx),(yvalsAC-mnbaseACy),'*')    
        plot((xvalsMU-mnbaseMUx),(yvalsMU-mnbaseMUy),'*','Color','r')        
        plot([0 300],[0 300],'k')     
% pk32bk28
        figure;hold on;
        plot(mean(AD12.exp(8).selectedpitchcurves')-mean(AD12.baselineINA'),'r')
        plot(mean(AD12.exp(11).selectedpitchcurves')-mean(AD12.baselineINA'),'r')
        plot(mean(AD12.exp(10).selectedpitchcurves')-mean(AD12.baselineAC'),'b')
        plot(mean(AD12.exp(6).selectedpitchcurves')-mean(AD12.baselineAC'),'b')
        % 420 is targ
        median(AD12.exp(5).toffset)
        
    shiftingAC=[];
    for i=[4 5 6 7 9 10 13]
        shiftingAC=[shiftingAC AD12.exp(i).selectedpitchcurves];
    end
    shiftingMU=[];
    for i=[8 11 14]
        shiftingMU=[shiftingMU AD12.exp(i).selectedpitchcurves];
    end
        
    figure;subplot(121);hold on; 
        yvalsAC=(runningaverage([mean(shiftingAC(380:420,:)) ],20));
        xvalsAC=(runningaverage([mean(shiftingAC(300:800,:)) ],20));
        plot(xvalsAC,yvalsAC,'*')    
        yvalsMU=(runningaverage([mean(shiftingMU(380:420,:)) ],20));
        xvalsMU=(runningaverage([mean(shiftingMU(300:800,:))],20));
        plot(xvalsMU,yvalsMU,'*','Color','r')        
        subplot(122);hold on;
        mnbaseACy=mean(mean(AD12.baselineAC(380:420,:)));
        mnbaseACx=mean(mean(AD12.baselineAC(300:800,:)));
        mnbaseMUy=mean(mean(AD12.baselineINA(380:420,:)));
        mnbaseMUx=mean(mean(AD12.baselineINA(300:800,:)));
        plot((xvalsAC-mnbaseACx),(yvalsAC-mnbaseACy),'*')    
        plot((xvalsMU-mnbaseMUx),(yvalsMU-mnbaseMUy),'*','Color','r')        
        plot([0 150],[0 150],'k')     

% pk16r47
% HERE,I COMPARE BOTH FROM AC BASELINE
% I'm not cherrypicking with the INAs - takes a while for reversion to
% occur and var reduction to occur
    figure;hold on;
    % second ina
    plot(mean(AD10.exp(11).selectedpitchcurves(:,end-50:end)')-mean(AD10.baselineAC'),'b')
    plot(mean(AD10.exp(8).selectedpitchcurves(:,end-50:end)')-mean(AD10.baselineAC'),'b')
    plot(mean(AD10.exp(12).selectedpitchcurves(:,10:end)')-mean(AD10.baselineAC'),'r') 
        median(targs)
    % first ina
    figure;hold on;
    hold on;plot(mean(AD10.exp(9).selectedpitchcurves(:,35:end)')-mean(AD10.baselineAC'),'r')
    plot(mean(AD10.exp(7).selectedpitchcurves')-mean(AD10.baselineAC'),'b')
        plot(mean(AD10.exp(8).selectedpitchcurves(:,end-20:end)')-mean(AD10.baselineAC'),'k')
            plot(mean(AD10.exp(8).selectedpitchcurves(:,end-40:end-20)')-mean(AD10.baselineAC'),'k')
  
    shiftingAC=[];
    for i=[7 8 10 11 13]
        shiftingAC=[shiftingAC AD10.exp(i).selectedpitchcurves];
    end
    shiftingMU=[];
        shiftingMU=[ AD10.exp(9).selectedpitchcurves(:,35:end) AD10.exp(12).selectedpitchcurves(:,10:end)];

      figure;subplot(121);hold on; 
        yvalsAC=(runningaverage([mean(shiftingAC(250:280,:)) ],20));
        xvalsAC=(runningaverage([mean(shiftingAC(250:450,:)) ],20));
        plot(xvalsAC,yvalsAC,'*')    
        yvalsMU=(runningaverage([mean(shiftingMU(250:280,:)) ],20));
        xvalsMU=(runningaverage([mean(shiftingMU(250:450,:))],20));
        plot(xvalsMU,yvalsMU,'*','Color','r')        
        subplot(122);hold on;
        mnbaseACy=mean(mean(AD10.baselineAC(250:280,:)));
        mnbaseACx=mean(mean(AD10.baselineAC(250:450,:)));
        mnbaseMUy=mean(mean(AD10.baselineINA(250:280,:)));
        mnbaseMUx=mean(mean(AD10.baselineINA(250:450,:)));
        plot((xvalsAC-mnbaseACx),(yvalsAC-mnbaseACy),'*')    
        plot((xvalsMU-mnbaseMUx),(yvalsMU-mnbaseMUy),'*','Color','r')        
        plot([0 150],[0 150],'k')               
% pk20r49
figure;hold on;
plot(mean(AD11.exp(16).selectedpitchcurves(:,100:end)')-mean(AD11.baselineAC'),'k')
plot(mean(AD11.exp(17).selectedpitchcurves')-mean(AD11.baselineINA'),'r')
hold on;plot(mean(AD11.exp(15).selectedpitchcurves(:,1:50)')-mean(AD11.baselineAC'),'k')
hold on;plot(mean(AD11.exp(12).selectedpitchcurves(:,end-20:end)')-mean(AD11.baselineAC'),'k')
plot(mean(AD11.exp(14).selectedpitchcurves')-mean(AD11.baselineINA'),'r')
hold on;plot(mean(AD11.exp(9).selectedpitchcurves')-mean(AD11.baselineAC'),'k')
hold on;plot(mean(AD11.exp(11).selectedpitchcurves')-mean(AD11.baselineINA'),'r')

  
    shiftingAC=[];
    for i=[7 9 10 12 13 15 16]
        shiftingAC=[shiftingAC AD11.exp(i).selectedpitchcurves];
    end
    shiftingMU=[];
    for i=[8 11 14 17]
        shiftingMU=[shiftingMU AD11.exp(i).selectedpitchcurves];
    end
      figure;subplot(121);hold on; 
        yvalsAC=(runningaverage([mean(shiftingAC(270:330,:)) ],20));
        xvalsAC=(runningaverage([mean(shiftingAC(200:400,:)) ],20));
        plot(xvalsAC,yvalsAC,'*')    
        yvalsMU=(runningaverage([mean(shiftingMU(270:330,:)) ],20));
        xvalsMU=(runningaverage([mean(shiftingMU(200:400,:))],20));
        plot(xvalsMU,yvalsMU,'*','Color','r')        
        subplot(122);hold on;
        mnbaseACy=mean(mean(AD11.baselineAC(270:330,:)));
        mnbaseACx=mean(mean(AD11.baselineAC(200:400,:)));
        mnbaseMUy=mean(mean(AD11.baselineINA(270:330,:)));
        mnbaseMUx=mean(mean(AD11.baselineINA(200:400,:)));
        plot((xvalsAC-mnbaseACx),(yvalsAC-mnbaseACy),'*')    
        plot((xvalsMU-mnbaseMUx),(yvalsMU-mnbaseMUy),'*','Color','r')        
        plot([0 150],[0 150],'k')   
        
        
 % r37g7
       figure;subplot(121);hold on; 
        yvalsAC=(runningaverage([mean(ShiftingAC(200:250,:)) ],20));
        xvalsAC=(runningaverage([mean(ShiftingAC(190:380,:)) ],20));
        plot(xvalsAC,yvalsAC,'*')    
        yvalsMU=(runningaverage([mean(ShiftingMU(200:250,:)) ],20));
        xvalsMU=(runningaverage([mean(ShiftingMU(190:380,:))],20));
        plot(xvalsMU,yvalsMU,'*','Color','r')        
        subplot(122);hold on;
        mnbaseACy=mean(mean(baselineAC(200:250,:)));
        mnbaseACx=mean(mean(baselineAC(190:380,:)));
        mnbaseMUy=mean(mean(baselineMU(200:250,:)));
        mnbaseMUx=mean(mean(baselineMU(190:380,:)));
        plot((xvalsAC-mnbaseACx),(yvalsAC-mnbaseACy),'*')    
        plot((xvalsMU-mnbaseMUx),(yvalsMU-mnbaseMUy),'*','Color','r')        
        plot([0 150],[0 150],'k')   