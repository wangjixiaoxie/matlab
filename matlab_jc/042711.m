% LMAN analysis
% 

% 04.26.11
% for each folder, enter basic information
    k=10;
    Alldata(k).folder='0428_wnon_L1_1hzto10khz';
    Alldata(k).RightLMAN=[2 5];
    Alldata(k).LeftLMAN=[3 4];
    Alldata(k).REF='left';
    Alldata(k).wnonrand=0;
    Alldata(k).wninstruct=0;
% pitch info, hit/escape info
    preHITms=0;
    Alldata(k).fv=findwnoteSPK('batchcatchnotes','a','','',0,[2000 2700],8000,1,'obs0',0,preHITms);
    fv=Alldata(k).fv;
    Alldata(k).pitch=pitchcontour(fv,2000,3000);    
    clear isTRIG isCATCH
    for i=1:length(fv);isTRIG(i)=fv(i).TRIG;end
    for i=1:length(fv);isCATCH(i)=fv(i).CATCH(1);end
    Alldata(k).wntarg=~isCATCH(find(isTRIG==1)); % for dividing targtimesWN
    Alldata(k).isHIT=find(isTRIG & ~isCATCH); % for identifying fvals
    Alldata(k).isESC=find(~isTRIG | isCATCH);
% neural info
    pretimems=1000;
    Alldata(k).Datt=getchandata('batchcatchnotes','a',[0:1:4],64000,pretimems);
 for i=1:length(Alldata(10).fv)
i
[y,f,t,p]=spectrogram(Alldata(10).Datt(2).data(i,:),1024,1020,512,32000);
Dattspect2(i).data=p(1:80,:);
end
for i=1:length(Alldata(10).fv)
i
[y,f,t,p]=spectrogram(Alldata(10).Datt(3).data(i,:),1024,1020,512,32000);
Dattspect3(i).data=p(1:80,:);
end
for i=1:length(Alldata(10).fv)
i
[y,f,t,p]=spectrogram(Alldata(10).Datt(4).data(i,:),1024,1020,512,32000);
Dattspect4(i).data=p(1:80,:);
end
for i=1:length(Alldata(10).fv)
i
[y,f,t,p]=spectrogram(Alldata(10).Datt(5).data(i,:),1024,1020,512,32000);
Dattspect5(i).data=p(1:80,:);
end
for i=1:length(Alldata(10).fv)
gamma2(i,:)=mean((Dattspect2(i).data(1:3,:)))-mean(mean(Dattspect2(i).data(1:3,1:3000)'));
gamma3(i,:)=mean((Dattspect3(i).data(1:3,:)))-mean(mean(Dattspect3(i).data(1:3,1:3000)'));
gamma4(i,:)=mean((Dattspect4(i).data(1:3,:)))-mean(mean(Dattspect4(i).data(1:3,1:3000)'));
gamma5(i,:)=mean((Dattspect5(i).data(1:3,:)))-mean(mean(Dattspect5(i).data(1:3,1:3000)'));
highpow2(i,:)=mean((Dattspect2(i).data(4:20,:)))-mean(mean(Dattspect2(i).data(4:20,1:3000)'));
highpow3(i,:)=mean((Dattspect3(i).data(4:20,:)))-mean(mean(Dattspect3(i).data(4:20,1:3000)'));
highpow4(i,:)=mean((Dattspect4(i).data(4:20,:)))-mean(mean(Dattspect4(i).data(4:20,1:3000)'));
highpow5(i,:)=mean((Dattspect5(i).data(4:20,:)))-mean(mean(Dattspect5(i).data(4:20,1:3000)'));
end   
    for i=1:length(Alldata)
        figure;hold on;
        fs=32000;
        [b,a]=butter(4,[100/fs]); % corresponds to smoothing of 100Hz - 10ms
        if Alldata(i).wnonrand==0
            clear  Mndata
            mdpt=round(size(Alldata(i).Datt(1).data,1)/2);
            for k=1:length(Alldata(i).Datt) % for each channel
                subplot(length(Alldata(i).Datt),1,k);hold on;
                num1=filtfilt(b,a,median(abs(Alldata(i).Datt(k).data(1:mdpt,:))));
                num2=filtfilt(b,a,median(abs(Alldata(i).Datt(k).data(mdpt+1:end,:))));
                plot(num1,'b');
                plot(num2,'r');
                ylim([min([num1(1000:end-1000) num2(1000:end-1000)]) max([num1(1000:end-1000) num2(1000:end-1000)])])
            end
        else
            for k=1:length(Alldata(i).Datt) % for each channel
                subplot(length(Alldata(i).Datt),1,k);hold on;
                num1=filtfilt(b,a,median(abs(Alldata(i).Datt(k).data(Alldata(i).isESC,:))));
                num2=filtfilt(b,a,median(abs(Alldata(i).Datt(k).data(Alldata(i).isHIT,:))));
                plot(num1,'b');
                plot(num2,'r');
                ylim([min([num1(1000:end-1000) num2(1000:end-1000)]) max([num1(1000:end-1000) num2(1000:end-1000)])])
            end
        end
    end
    %
    
    figure;hold on;
%     subplot(length(currdata)+1,1,1)
%     %imagesc(t,f,log(avA));syn;ylim([0 1e4])
%     xlim([-0.2 0.4])
    for i=1:length(currdata)
        subplot(length(currdata)+1,1,i+1);hold on;
%         plot([avgwin/fs+1/fs:1/fs:size(currdata(i).data,2)/fs],ravesc(i,:))
%         plot([avgwin/fs+1/fs:1/fs:size(currdata(i).data,2)/fs],ravhit(i,:),'r')  
        plot([avgwin/fs+1/fs:1/fs:size(currdata(i).data,2)/fs],rav1(i,:),'b')          
%         plot([36000/fs 36000/fs],[min([ravhit(i,:) ravesc(i,:)])*0.9 max([ravhit(i,:) ravesc(i,:)])*1.1],'k')
       plot([36000/fs 36000/fs],[min([rav1(i,:)])*0.9 max([rav1(i,:)])*1.1],'k')        
%         xlim([0.8 1.4]); ylim([min([ravhit(i,:) ravesc(i,:)])*0.9 max([ravhit(i,:) ravesc(i,:)])*1.1])
        xlim([0.8 1.4]); ylim([min([rav1(i,:)])*0.9 max([rav1(i,:)])*1.1])        
    end
%
pitchAll2=pitchcontour(fv,2000,3000);
pts=[-0.8:0.01:0.8];
clear ga
for i=1:length(pts)
    cc=corrcoef(pitchvals,mean(abs(currdata(1).data(:,1.12*fs+pts(i)*fs-0.1*fs:1.12*fs+pts(i)*fs))'));
    ga(i)=cc(2);
end
ind3=find(pitchAll(900,:)>2220 & pitchAll(900,:)<2480);
ind3=ind3(1:120);
pitchvals=mean(pitchAll2(900:1100,:));
figure;plot(pitchAll(900,ind3),mean(abs(Datt422(2).data(ind3,1.12*fs-0.1*fs:1.12*fs))'),'.')
corrcoef(pitchAll(900,ind3),mean(abs(Datt422(2).data(ind3,1.12*fs-0.1*fs:1.12*fs))'))
figure;plot(pitchAll(900,ind3),mean(abs(Datt422(3).data(ind3,1.12*fs-0.1*fs:1.12*fs))'),'.')
corrcoef(pitchAll(900,ind3),mean(abs(Datt422(3).data(ind3,1.12*fs-0.1*fs:1.12*fs))'))
