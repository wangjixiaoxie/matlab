% bk63
% bk61
% bk50 (ctl)
% pu56
% pu57


% TSA720
% /cardinal/bk63w43/wnon1205newtemp/batch06jc
pitchTSAhar1=jc_pitchmat1024(shiftedTSA,1024,1020,1,1800,2700,[1],'obs0',1);
figure;plot(median(pitchTSAhar1(390:425,:)),'*')
% Look for a bimodal distribution
figure;hist(std(pitchTSAhar1(550:700,:)))
indE=find((std(pitchTSAhar1(550:700,:))<100));
figure;plot(pitchTSAhar1(:,indE))

figure;plot(indE,median(pitchTSAhar1(420:460,indE)),'Color','r')
hold on;plot(indH,median(pitchTSAhar1(420:460,indH)),'Color','k')
%%%%
medianpitch=median(pitchTSAhar1(420:460,:));
%%%%
countHit=zeros(1,50);
countEsc=zeros(1,50);
facHit=[];
facEsc=[];
for i=51:size(pitchTSAhar1,2)
    res=medianpitch(i)-median(medianpitch(i-50:i-30));
    for k=1:50
        j=51-k;
        index=i-j;
        if isempty(find(indE==index))
            %Hit=median((pitchTSAhar1(420:460,index)-(pitchTSAhar1(420:460,i-2)))');
            % Does it go away from the hits?  If so, these should be
            % negative values.
            countHit(j)=countHit(j)+1;
            facHit(j,countHit(j))=res;
            
        else
            %Esc=median((pitchTSAhar1(420:460,index)-(pitchTSAhar1(420:460,i-2)))');
            countEsc(j)=countEsc(j)+1;
            facEsc(j,countEsc(j))=res;
        end
    end
end
figure;subplot(211);plot(median(facHit(:,1:199)'))
hold on;plot(median(facEsc(:,1:78)'),'r')
hold on;plot([0 50],[0 0],'k')
hold on;plot(median(facHit(:,1:199)'),'.','MarkerSize',15)
hold on;plot(median(facEsc(:,1:78)'),'.','MarkerSize',15,'Color','r')
subplot(312);hold on;
plot(median(pitchTSAhar1(420:460,indH)'))
plot(median(pitchTSAhar1(420:460,indE)'),'r')
subplot(313);hold on;
plot(indE,median(pitchTSAhar1(420:460,indE)),'Color','r')
hold on;plot(indH,median(pitchTSAhar1(420:460,indH)),'Color','b')

%%%%%%%%
countHit=zeros(1,50);
countEsc=zeros(1,50);
facHit=[];
facEsc=[];
for i=51:size(pitchTSAhar1,2)
    res=medianpitch(i)-median(medianpitch(i-50:i-30));
    for k=1:50
        j=51-k;
        index=i-j;
        if isempty(find(indE==index))
            Hit=medianpitch(index)-median(medianpitch(i-50:i-30));
            % Does it go away from the hits?  If so, these should be
            % negative values.
            countHit(j)=countHit(j)+1;
            facHit(j,countHit(j))=abs((res-Hit)/Hit);
            
        else
            Esc=medianpitch(index)-median(medianpitch(i-50:i-30));
            countEsc(j)=countEsc(j)+1;
            facEsc(j,countEsc(j))=abs((res-Esc)/Esc);
        end
    end
end
figure;plot(median(facHit(:,1:199)'))
hold on;plot(median(facEsc(:,1:78)'),'r')
hold on;plot([0 40],[0 0],'k')
hold on;plot(median(facHit(:,1:199)'),'.','MarkerSize',15)
hold on;plot(median(facEsc(:,1:78)'),'.','MarkerSize',15,'Color','r')


%%%%%%%%%
%%%%%%%%%
% bk50 no signal - random
% /cardinal/bk50w18_trialbytrial/Reinforcement/batchNOTE
pitchX=jc_pitchmat1024(shiftedX,1024,1020,1,2000,2700,[1],'obs0',1);
figure;plot(median(pitchX(650:750,:)))
figure;hist((std(pitchX(900:970,:))))
indE=find((std(pitchX(900:970,:))<100));
indH=find((std(pitchX(900:970,:))>100));
medianpitchX=median(pitchX(700:750,:));
%%%%
countHit=zeros(1,50);
countEsc=zeros(1,50);
facHit=[];
facEsc=[];
for i=51:size(pitchX,2)
    res=medianpitchX(i)-mean(medianpitchX(:));
    for k=1:50
        j=51-k;
        index=i-j;
        if isempty(find(indE==index))
            Hit=medianpitchX(index)-mean(medianpitchX(:));
            % Does it go away from the hits?  If so, these should be
            % negative values.
            countHit(j)=countHit(j)+1;
            facHit(j,countHit(j))=abs((res-Hit)/Hit);
            
        else
            Esc=medianpitchX(index)-mean(medianpitchX(:));
            countEsc(j)=countEsc(j)+1;
            facEsc(j,countEsc(j))=abs((res-Esc)/Esc);
        end
    end
end
figure;plot(median(facHit(:,1:833)'))
hold on;plot(median(facEsc(:,1:244)'),'r')
hold on;plot([0 40],[0 0],'k')
hold on;plot(median(facHit(:,1:833)'),'.','MarkerSize',15)
hold on;plot(median(facEsc(:,1:244)'),'.','MarkerSize',15,'Color','r')
%%%%%%%%
%%%
% pu56
%%% /cardinal5/pu56w26/Timshift/222tmp2/batchJC2
%%%%%%%%
indE=find((std(pitchX(260:320,:))<100));
medianpitch=median(pitchX(180:200,:));
%%%%
countHit=zeros(1,50);
countEsc=zeros(1,50);
facHit=[];
facEsc=[];
for i=51:size(pitchX,2)
    res=medianpitch(i)-median(medianpitch(i-20:i-2));
    for k=1:50
        j=51-k;
        index=i-j;
        if isempty(find(indE==index))
            %Hit=median((pitchTSAhar1(420:460,index)-(pitchTSAhar1(420:460,i-2)))');
            % Does it go away from the hits?  If so, these should be
            % negative values.
            countHit(j)=countHit(j)+1;
            facHit(j,countHit(j))=res;
            
        else
            %Esc=median((pitchTSAhar1(420:460,index)-(pitchTSAhar1(420:460,i-2)))');
            countEsc(j)=countEsc(j)+1;
            facEsc(j,countEsc(j))=res;
        end
    end
end
figure;plot(median(facHit(:,1:150)'))
hold on;plot(median(facEsc(:,1:80)'),'r')
hold on;plot([0 50],[0 0],'k')
hold on;plot(median(facHit(:,1:150)'),'.','MarkerSize',15)
hold on;plot(median(facEsc(:,1:80)'),'.','MarkerSize',15,'Color','r')
% subplot(312);hold on;
% plot(median(pitchTSAhar1(420:460,indH)'))
% plot(median(pitchTSAhar1(420:460,indE)'),'r')
% subplot(313);hold on;
% plot(indE,median(pitchTSAhar1(420:460,indE)),'Color','r')
% hold on;plot(indH,median(pitchTSAhar1(420:460,indH)),'Color','b')
% 
% 
%%%%%%%%
countHit=zeros(1,50);
countEsc=zeros(1,50);
facHit=[];
facEsc=[];
for i=51:size(pitchX,2)
    res=medianpitch(i)-median(medianpitch(i-50:i-40));
    for k=1:50
        j=51-k;
        index=i-j;
        if isempty(find(indE==index))
            Hit=medianpitch(index)-median(medianpitch(i-50:i-40));
            % Does it go away from the hits?  If so, these should be
            % negative values.
            countHit(j)=countHit(j)+1;
            facHit(j,countHit(j))=abs((res-Hit)/Hit);
            
        else
            Esc=medianpitch(index)-median(medianpitch(i-50:i-40));
            countEsc(j)=countEsc(j)+1;
            facEsc(j,countEsc(j))=abs((res-Esc)/Esc);
        end
    end
end
figure;plot(median(facHit(:,1:150)'))
hold on;plot(median(facEsc(:,1:80)'),'r')
hold on;plot([0 40],[0 0],'k')
hold on;plot(median(facHit(:,1:150)'),'.','MarkerSize',15)
hold on;plot(median(facEsc(:,1:80)'),'.','MarkerSize',15,'Color','r')


%%%%%%%%%
%%%%%%%%%% 
% /cardinal5/pu57w52/21209_ACSF
figure;hist(median(abs(shiftedX(:,2000:3000))'))
indE=find((median(abs(shiftedX(:,2000:3000))'))<2500)
indH=find((median(abs(shiftedX(:,2000:3000))'))>2500)
figure;plot(indE,median(pitchX(240:260,indE)),'*','Color','r')
hold on;plot(indH,median(pitchX(240:260,indH)),'*','Color','k')
%%%%
medianpitch=median(pitchX(240:260,:));
%%%%
countHit=zeros(1,50);
countEsc=zeros(1,50);
facHit=[];
facEsc=[];
for i=51:280%size(pitchX,2)
    res=medianpitch(i)-median(medianpitch(i-50:i-30));
    for k=1:50
        j=51-k;
        index=i-j;
        if isempty(find(indE==index))
            %Hit=median((pitchTSAhar1(420:460,index)-(pitchTSAhar1(420:460,i-2)))');
            % Does it go away from the hits?  If so, these should be
            % negative values.
            countHit(j)=countHit(j)+1;
            facHit(j,countHit(j))=res;
            
        else
            %Esc=median((pitchTSAhar1(420:460,index)-(pitchTSAhar1(420:460,i-2)))');
            countEsc(j)=countEsc(j)+1;
            facEsc(j,countEsc(j))=res;
        end
    end
end
figure;plot(median(facHit(:,1:90)'))
hold on;plot(median(facEsc(:,1:134)'),'r')
hold on;plot([0 50],[0 0],'k')
hold on;plot(median(facHit(:,1:90)'),'.','MarkerSize',15)
hold on;plot(median(facEsc(:,1:134)'),'.','MarkerSize',15,'Color','r')

%%%%
countHit=zeros(1,50);
countEsc=zeros(1,50);
facHit=[];
facEsc=[];
for i=51:280%size(pitchX,2)
    res=medianpitch(i)-median(medianpitch(i-50:i-30));
    for k=1:50
        j=51-k;
        index=i-j;
        if isempty(find(indE==index))
            Hit=medianpitch(index)-median(medianpitch(i-50:i-30));
            % Does it go away from the hits?  If so, these should be
            % negative values.
            countHit(j)=countHit(j)+1;
            facHit(j,countHit(j))=abs((res-Hit)/Hit);
            
        else
            Esc=medianpitch(index)-median(medianpitch(i-50:i-40));
            countEsc(j)=countEsc(j)+1;
            facEsc(j,countEsc(j))=abs((res-Esc)/Esc);
        end
    end
end
figure;plot(median(facHit(:,1:90)'))
hold on;plot(median(facEsc(:,1:134)'),'r')
hold on;plot([0 40],[0 0],'k')
hold on;plot(median(facHit(:,1:90)'),'.','MarkerSize',15)
hold on;plot(median(facEsc(:,1:134)'),'.','MarkerSize',15,'Color','r')

%%%%%%%%%%
% /cardinal/bk61w42/wnon1205newtemp
% 12/6/08 - batchTSAnotes
        figure;hist(abs(median(shiftedX(:,2000:6000)')))
        indE=find((median(abs(shiftedX(:,2000:6000))'))<300)
        indH=find((median(abs(shiftedX(:,2000:6000))'))>300)
        medianpitch=median(pitchX(380:400,:));
        figure;plot(indE,medianpitch(indE),'r')
        hold on;plot(indH,medianpitch(indH),'k')

        %%%%%%%%%%

        scalarFF(indE,medianpitch)
        %%%%
        %%%%%%%%%
        scalarFFattract(indE,medianpitch)
        %%%%%%%%%%
        
% /cardinal2/PitchShiftData/bk20bk45/wnon713
% bk20bk45 note B - wnon
medianpitch=median(pitchX(250:270,:));
indE=find(mean(abs(shiftedX(:,1500:4000))')<2000)
indH=find(mean(abs(shiftedX(:,1500:4000))')>2000)
        figure;plot(indE,medianpitch(indE),'r')
        hold on;plot(indH,medianpitch(indH),'k')
scalarFF(indE,medianpitch)
        scalarFFattract(indE,medianpitch)