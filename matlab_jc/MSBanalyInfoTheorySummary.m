load /cardinal/FiguresA/MSBinfoApril2010.mat

% Rationale - bird learns to shift upwards/downwards at the times when there are the
    % highest positive/negative correlation between behavior (i.e. pitch)
    % and outcome (i.e. probability of white noise).  In single contingency
    % experiments, there is only one time.  In dual contingency experiments
    % there are two times.
% Single contingency experiments
        for j=1:10
            for i=1:28
                if isequal(Predict(i).direction,'up')
                    valu=60;
                else
                    valu=40;
                end
                if i<23
                    [SCpoint(i,j)]=ContingSimMSBSingle(Predict(i).Targeting-Predict(i).onset,Predict(i).ResidAC(Predict(i).onset:Predict(i).onset+400,:),valu);
                else
                    [SCpoint(i,j)]=ContingSimMSBSingle(Predict(i).Targeting-Predict(i).onset,Predict(i).ResidAC(Predict(i).onset:Predict(i).offset,:),valu);
                end
                SCactual(i)=round(median(Predict(i).Targeting-Predict(i).onset));
            end
        end
        
        
   % processing...     
   halves=[40 42 44 46 48 50 52 54 56 58 60 62 64 66 68 70 40];
   for k=1:length(halves)
       
        SCpointAVG=mean(SCpoint');
        % mean(SCpointAVG-SCactual)=-21.6 (less than 3ms);
        % mean absolute value is 34.2 - which is ~4ms - within window
       % figure;hold on;
       halfwidth=halves(k);
        allpts=zeros(28,1400);
        for i=1:28
            point1=round(median(Predict(i).Targeting)-Predict(i).onset);
            point2=round(Predict(i).offset-median(Predict(i).Targeting));
            %plot([-1*point1:1:point2],gaussian([-1*point1:1:point2],SCpointAVG(i)-SCactual(i),halfwidth*8))
            allpts(i,700-point1:700+point2)=gaussian([-1*point1:1:point2],SCpointAVG(i)-SCactual(i),halfwidth*8);
        end
    for i=1:1400
        mnSCpred(i)=mean(allpts(allpts(:,i)>0,i));
    end
    mnSCpred=mnSCpred/max(mnSCpred(500:900));
    dist(k)=sum(abs(SSactual(500:900)-mnSCpred(500:900)));
  
   end
   
    
    
% Dual contingency experiments
% FRONT & REAR
    for i=1:8
        if isequal(DShifts(i).dirA,'up')
            valu=51;
        else
            valu=49;
        end
        [ITfirst(i) ITsecond(i)]=ContingSimMSBDual(DShifts(i).toffset-DShifts(i).onset,jc_residuals(DShifts(i).pitchBaseline(DShifts(i).onset:DShifts(i).onset+400,:)),valu);
        Actualfirst(i)=round(median(DShifts(i).toffset-DShifts(i).onset));
        Actualsecond(i)=(Actualfirst(i)+192);
    end

  % processing...
  halves=[30:0.1:35];
  for k=1:length(halves)
    halfwidth=halves(k);
    gaussShape=gaussian([1:1:1400],700,halfwidth*8);
    allptsDC=zeros(8,1400);
    for i=1:8
        lengthnote=DShifts(i).offset-DShifts(i).onset+1;
        ptB=round(median(DShifts(i).toffset)+192-DShifts(i).onset);
        allptsDCend=gaussShape(700-ITsecond(i):700-ITsecond(i)+lengthnote);
        allptsDCbegin=gaussShape(700-ITfirst(i):700-ITfirst(i)+lengthnote);
        apDC=(allptsDCend-allptsDCbegin);
        allptsDC(i,700-ptB:700-ptB+lengthnote)=apDC/max(apDC(200:end));
    end
    for i=1:1400
        mnDCpred(i)=mean(allptsDC(allptsDC(:,i)~=0,i));
    end
    %cd=corrcoef(DSactual,mnDCpred(400:800));
    dist(k)=sum(abs(DSactual-mnDCpred(400:800)));
  end
figure;plot(halves,dist)
hold on;plot([0 60],[0.97 0.97])