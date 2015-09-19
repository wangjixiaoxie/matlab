%summary anova stats;
clear muvls stvls anphvls anstvls
muvls=plotvlsph.mupct;
lnph=length(muvls(:,1));
stvls=plotvlsst.mupct;
lnst=length(stvls(:,1));
for ii=1:3
    anphvls(:,ii)=muvls(:,ii)-muvls(:,1);
    anstvls(:,ii)=stvls(:,ii)-stvls(:,1);
end

y=[anphvls(:,1);anphvls(:,2);anphvls(:,3);anstvls(:,1);anstvls(:,2);anstvls(:,3)]
tm=[zeros(lnph,1);ones(lnph,1);2*ones(lnph,1);zeros(lnst,1);ones(lnst,1);2*ones(lnst,1)];
type=[zeros(3*lnph,1);ones(3*lnst,1)]

P=anovan(y,{tm type})

clear anphvls anstvls

muvls=phpre.mu./phpre.ac;
lnph=length(muvls(:,1));
stvls=stpre.mu./stpre.ac;
lnst=length(stvls(:,1));
for ii=1:2
    anphvls(:,ii)=muvls(:,ii)-muvls(:,1);
    anstvls(:,ii)=stvls(:,ii)-stvls(:,1);
end

y=[anphvls(:,1);anphvls(:,2);anstvls(:,1);anstvls(:,2)]
tm=[zeros(lnph,1);ones(lnph,1);zeros(lnst,1);ones(lnst,1)];
type=[zeros(2*lnph,1);ones(2*lnst,1)]

Ppre=anovan(y,{tm type})
