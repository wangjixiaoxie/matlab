function prediction=jccrossco2(residuals,targdistn)
% Predict(1).crosscorrAC=jccrossco(Predict(1).ResidAC,Predict(1).Targeting);

radius=32; % half of width of targeting windwo
targdistn=targdistn-radius;
clear crossco;

    Targetrandomdraw=round(targdistn);
    mm=mean(residuals(Targetrandomdraw-radius:Targetrandomdraw+radius,:));
    for j=1:size(residuals,1)
        nn=corrcoef(mm,residuals(j,:));
        crossco(j+1)=nn(2);
    end

prediction=(crossco);
