ind2=[1:21];
for i=1:length(ind2)
    i
    if ~isempty(Exp(i).fv2)
        pt2=pitchcontour(Exp(i).fv2,Exp(i).FF(1),Exp(i).FF(2));
        Exp(i).t2=timing3(Exp(i).fv2);
        if length(Exp(i).window)>1
            Exp(i).pt2=mean(pt2(Exp(i).window,:));
        else
            Exp(i).pt2=(pt2(Exp(i).window,:));
        end
    end
end
% for i=1:length(EPC)
%     if length(ExperimentPC(i).time)>1
%         EPC(i).baseline=mean(mean(ExperimentPC(i).pitchPre(ExperimentPC(i).time,:)));
%     else
%         EPC(i).baseline=mean(ExperimentPC(i).pitchPre(ExperimentPC(i).time,:));
%     end
% end


%%%%
load /bulbul2/ImmediacyALL.mat

        clear INTvals TIMEvalsCTL
        for i=1:length(EPC)
            direction=1-2*(isequal(EPC(i).DIR,'down'));
            n=1/EPC(i).samprate;
            vls=[1:n:length(EPC(i).FF)*n];
            X=vls;
            Y=[direction*(EPC(i).FF-(EPC(i).baseline))/(EPC(i).baseline)];
            INTvals(i,:)=interp1(X,Y,[1:1:550],'linear');
        end
        clear mnPC sizePC sdPC mnTIMEPC
        for i=1:size(INTvals,2)
            mnPC(i)=mean(INTvals(find(~isnan(INTvals(:,i))),i));
            sizePC(i)=length(find(~isnan(INTvals(:,i))));
            sdPC(i)=jcstd(INTvals(find(~isnan(INTvals(:,i))),i))/sqrt(sizePC(i));    
        end
        figure;plot(runningaverage(mnPC-sdPC,50));hold on;plot(runningaverage(mnPC+sdPC,50));plot(runningaverage(mnPC,50),'Linewidth',3)
%
for i=1:length(EPC)
    ttime(i)=max(EPC(i).time)-min(EPC(i).time);
    tnotes(i)=length(EPC(i).FF)/EPC(i).samprate;
end
mean(ttime./tnotes) % One syllable every 0.0287 hours = 1.7243 minutes
        %%%%
clear INTvals2
for i=1:length(Exp)
    direction=1-2*(isequal(Exp(i).DIR,'down'));
    LN1=length(Exp(i).pt1);
    INTvals2(i,1:LN1)=(direction*(Exp(i).pt1-Exp(i).baseline))/Exp(i).baseline;
    vls=[LN1+1:5:LN1+length(Exp(i).pt2)*5];
    if ~isempty(Exp(i).pt2)
        X=vls;
        Y=[direction*(Exp(i).pt2-(Exp(i).baseline))/(Exp(i).baseline)];
        INTvals2(i,LN1+1:1:700)=interp1(X,Y,[LN1+1:1:700],'linear');
    end
end
        clear mn2 size2 sd2
        for i=1:size(INTvals2,2)
            mn2(i)=mean(INTvals2(find(INTvals2(:,i)~=0 & ~isnan(INTvals2(:,i))),i));
            size2(i)=length(find(INTvals2(:,i)~=0 & ~isnan(INTvals2(:,i))));
            sd2(i)=jcstd(INTvals2(find(INTvals2(:,i)~=0 & ~isnan(INTvals2(:,i))),i))/sqrt(size2(i));    
        end
%        
clear ttime2 tnotes2
for i=1:length(Exp)
    if ~isempty(Exp(i).t2)
        ttime2(i)=[max(Exp(i).t2)-min(Exp(i).t1)];
        tnotes2(i)=length(Exp(i).pt1)+length(Exp(i).pt2)*5;
    else
        ttime2(i)=[max(Exp(i).t1)-min(Exp(i).t1)];
        tnotes2(i)=length(Exp(i).pt1);        
    end
end
mean(ttime2./tnotes2) % One syllable every 0.0573 hours = 3.44 minutes
    
        
        
     % CV recovery
         indX=[3:5 8 11:12 14:21]; % n=14
        INTvals3=INTvals2(indX,:);
        clear mn3 size3 sd3
        for i=1:size(INTvals3,2)
            mn3(i)=mean(INTvals3(find(INTvals3(:,i)~=0 & ~isnan(INTvals3(:,i))),i));
            size3(i)=length(find(INTvals3(:,i)~=0 & ~isnan(INTvals3(:,i))));
            sd3(i)=jcstd(INTvals3(find(INTvals3(:,i)~=0 & ~isnan(INTvals3(:,i))),i))/sqrt(size3(i));    
        end
        
  %%   Good comparison 
        figure;hold on;
        subplot(121);
        plot(runningaverage(mn2-sd2,50),'color','r');hold on;plot(runningaverage(mn2+sd2,50),'color','r');plot(runningaverage(mn2,50),'Linewidth',3,'color','r');
        xlim([0 445]);ylim([-0.01 0.02]);plot([0 445],[0 0]);plot([0 480],[mean(mean(INTvals2(:,1:20))) mean(mean(INTvals2(:,1:20)))])
        subplot(122);hold on;plot(runningaverage(mnPC-sdPC,50));hold on;plot(runningaverage(mnPC+sdPC,50));plot(runningaverage(mnPC,50),'Linewidth',3)
        xlim([0 445]);ylim([-0.01 0.02]);plot([0 445],[0 0]);plot([0 480],[mean(mean(INTvals2(:,1:20))) mean(mean(INTvals2(:,1:20)))])
  %%  BEST - less information
        index1=[1 1 51 101 151 201 251 301 351 401 451];
        index2=[5 index1(2:end)+49];
        a=[];b=[];c=[];a1=[];b1=[];c1=[];
        for i=1:length(index1)
            a=[a 100*mean(mnPC(index1(i):index2(i)))];
            b=[b 100*mean(mnPC(index1(i):index2(i)))-100*mean(sdPC(index1(i):index2(i)))];
            c=[c 100*mean(mnPC(index1(i):index2(i)))+100*mean(sdPC(index1(i):index2(i)))];    
            a1=[a1 100*mean(mn2(index1(i):index2(i)))];
            b1=[b1 100*mean(mn2(index1(i):index2(i)))-100*mean(sd2(index1(i):index2(i)))];
            c1=[c1 100*mean(mn2(index1(i):index2(i)))+100*mean(sd2(index1(i):index2(i)))];    
        end
        figure;hold on;subplot(121);hold on;
        plot(mean([index1;index2]),a);plot(mean([index1;index2]),b);plot(mean([index1;index2]),c);
        plot(mean([index1;index2]),a,'.','Markersize',15);plot(mean([index1;index2]),b,'.','Markersize',15);plot(mean([index1;index2]),c,'.','Markersize',15);
        plot([-20 500],[0 0],'k');ylim([-0.5 2]);xlim([-20 500]);plot([-20 500],[0.99 0.99],'k')        
        subplot(122);hold on;
        plot(mean([index1;index2]),a1,'r');plot(mean([index1;index2]),b1,'r');plot(mean([index1;index2]),c1,'r');
        plot(mean([index1;index2]),a1,'r.','Markersize',15);plot(mean([index1;index2]),b1,'r.','Markersize',15);plot(mean([index1;index2]),c1,'r.','Markersize',15);
        plot([-20 500],[0 0],'k');ylim([-0.5 2]);xlim([-20 500]);plot([-20 500],[0.99 0.99],'k')
 %% 
  % How long until it gets to the first n hours?
  
  
  %%
figure;hold on;plot(mean(INTvals2(:,1:50)),'Linewidth',3)
plot(mean(INTvals2(:,1:50))-std(INTvals2(:,1:50))/sqrt(21));
plot(mean(INTvals2(:,1:50))+std(INTvals2(:,1:50))/sqrt(21));
plot([0 20],[0 0],'k')
xlim([1 20])

figure;hold on;plot(mean(INTvals(:,1:50)),'Linewidth',3)
plot(mean(INTvals(:,1:50))-std(INTvals(:,1:50))/sqrt(14));
plot(mean(INTvals(:,1:50))+std(INTvals(:,1:50))/sqrt(14));
plot([0 20],[0 0],'k')
xlim([1 20])
