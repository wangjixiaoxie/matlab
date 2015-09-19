function val=specent(fdat);
%

fdat=abs(fdat);

NN=length(fdat);
tmp1=exp((1.0./NN).*sum(log(fdat)));
tmp2=mean(fdat);
val=tmp1./tmp2;
return;
