function prediction=jccrossco(residuals,targdistn)
% Predict(1).crosscorrAC=jccrossco(Predict(1).ResidAC,Predict(1).Targeting);

radius=32; % half of width of targeting windwo
targdistn=targdistn-radius;
clear crossco;
for k=1:100
    Targetrandomdraw=round(rand*(length(targdistn)-1)+1);
    mm=mean(residuals(round(targdistn(Targetrandomdraw))-radius:round(targdistn(Targetrandomdraw))+radius,:));
    for j=1:1700
        nn=corrcoef(mm,residuals(j,:));
        crossco(k,j+1)=nn(2);
    end
end
prediction=median(crossco);
