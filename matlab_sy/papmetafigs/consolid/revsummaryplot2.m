
function [statvls]=revsummaryplot(sv,ph,bs)
bascomb=[];
shmucomb=[];
shaccomb=[];
conaccomb=[];
conmucomb=[];
recaccomb=[];
recmucomb=[];
finaccomb=[];
finmucomb=[];
%indices with natural recovery
NAT=[1 7]
for ii=1:length(sv)
   %first define the direction of the shift.... 
   %easiest way to do this is to check the shift value
   %and find which of the shiftruns it is in.
   bsnum=sv(ii).bsnum
   
   phcr=ph(bsnum)
   sfact=bs(bsnum).sfact(1)/1000;
   shruns=phcr.shiftruns
   shind=sv(ii).sh
   clear shnum
   for jj=1:length(shruns)
       if(ismember(shind,shruns{jj}))
           shnum=jj
           drxn=phcr.drxn{shnum};
           if drxn=='up'
               FAC=1;
           else
               FAC=-1;
           end
       
       end
   end
   basind=sv(ii).bas;
   conind=sv(ii).con;
   recind=sv(ii).rec;
   finind=sv(ii).fin;
   basoutac=0;
   basvl=phcr.acprez(phcr.ntind,basind);
   
   crsd=phcr.initsd(phcr.ntind);
   
      con_ac=((FAC*(phcr.acprez(phcr.ntind,conind)-basvl)*crsd)./sfact);
   con_mu=((FAC*(phcr.muz(phcr.ntind,conind)-basvl)*crsd)./sfact);
   NFAC=con_ac;
   
   
   conaccomb=[conaccomb con_ac./NFAC]
   conmucomb=[conmucomb con_mu./NFAC]
   
   
   bas_mu=(FAC*(phcr.muz(phcr.ntind,basind)*crsd)./sfact);
   bascomb=[bascomb bas_mu./NFAC];
   
   sh_ac=((FAC*(phcr.acprez(phcr.ntind,shind)-basvl)*crsd)./sfact);
   sh_mu=((FAC*(phcr.muz(phcr.ntind,shind)-basvl)*crsd)./sfact);
   shmucomb=[shmucomb sh_mu./NFAC];
   shaccomb=[shaccomb sh_ac./NFAC];
   
   
   con_ac=((FAC*(phcr.acprez(phcr.ntind,conind)-basvl)*crsd)./sfact);
   con_mu=((FAC*(phcr.muz(phcr.ntind,conind)-basvl)*crsd)./sfact);
   
   
   rec_ac=mean(((FAC*(phcr.acprez(phcr.ntind,recind)-basvl)*crsd)./sfact));
   rec_mu=mean(((FAC*(phcr.muz(phcr.ntind,recind)-basvl)*crsd)./sfact));
   recaccomb=[recaccomb rec_ac./NFAC]
   recmucomb=[recmucomb rec_mu./NFAC]
     if(~isempty(finind))
         fin_ac=((FAC*(phcr.acprez(phcr.ntind,finind)-basvl)*crsd)./sfact);
        fin_mu=((FAC*(phcr.muz(phcr.ntind,finind)-basvl)*crsd)./sfact);
        finaccomb=[finaccomb fin_ac./NFAC];
        finmucomb=[finmucomb fin_mu./NFAC];
     end
end
 
%modified so that plotting diff...
figure
colvl=[0.5 0.5 0.5]
%baseline
df=-bascomb
bar(2,mean(df),0.6,'EdgeColor',colvl,'FaceColor',colvl);
hold on;
stvl=std(df)./sqrt(length(df))
plot([2 2],[mean(df)-stvl mean(df)+stvl],'k');
% plot(2*ones(length(bascomb),1),bascomb,'ro');

%shift
crdf=-shaccomb+shmucomb
bar(4,mean(crdf),0.6,'EdgeColor',colvl,'FaceColor',colvl)
stvl=std(crdf)./sqrt(length(crdf));
plot([4 4],[mean(crdf)-stvl mean(crdf)+stvl],'k');
% plot(4*ones(length(shaccomb),1),shaccomb,'ko');

% bar(5,mean(shmucomb),0.6,'EdgeColor','r','FaceColor','none')
% stvl=std(shmucomb)./sqrt(length(shmucomb));
% plot([5 5],[mean(shmucomb)-stvl mean(shmucomb)+stvl],'r');
% plot(5*ones(length(shmucomb),1),shmucomb,'ro');

%consolid
crdf=-conaccomb+conmucomb;
bar(6,mean(crdf),0.6,'EdgeColor',colvl,'FaceColor',colvl)
stvl=std(crdf)./sqrt(length(crdf));
plot([6 6],[mean(crdf)-stvl mean(crdf)+stvl],'k');
% plot(7*ones(length(conaccomb),1),conaccomb,'ko');

% bar(8,mean(conmucomb),0.6,'EdgeColor',colvl,'FaceColor',colvl)
% stvl=std(conmucomb)./sqrt(length(conmucomb));
% plot([8 8],[mean(conmucomb)-stvl mean(conmucomb)+stvl],'r');
% plot(8*ones(length(conmucomb),1),conmucomb,'ko');

%rec
%get notNAT
allind=1:length(recaccomb);
NOTNAT=setdiff(allind,NAT);
crdiffnat=-recaccomb(NAT)+recmucomb(NAT);
crdiffnotnat=-recaccomb(NOTNAT)+recmucomb(NOTNAT);
bar(7.5,mean(crdiffnat),0.6,'FaceColor',colvl,'EdgeColor',colvl)
stvl=std(crdiffnat)./sqrt(length(crdiffnat));
plot([7.5 7.5],[mean(crdiffnat)-stvl mean(crdiffnat)+stvl],'k');

bar(8.5,mean(crdiffnotnat),0.6,'FaceColor',colvl,'EdgeColor',colvl)
stvl=std(crdiffnotnat)./sqrt(length(crdiffnotnat));
plot([8.5 8.5],[mean(crdiffnotnat)-stvl mean(crdiffnotnat)+stvl],'k');


% plot(10*ones(length(recaccomb),1),recaccomb,'ko');
% plot(10*ones(length(NAT),1),recaccomb(NAT),'ko','MarkerEdgeColor','k','MarkerFaceColor','k');
% 
% 
% bar(11,mean(recmucomb),0.6,'FaceColor','none','EdgeColor','r')
% stvl=std(recmucomb)./sqrt(length(recmucomb));
% plot([11 11],[mean(recmucomb)-stvl mean(recmucomb)+stvl],'r');
% plot(11*ones(length(recmucomb),1),recmucomb,'ro');
% plot(11*ones(length(NAT),1),recmucomb(NAT),'ro','MarkerFaceColor','r');
% 
% plot([10 11],[recaccomb(1) recmucomb(1)]);
% plot([10 11],[recaccomb(7) recmucomb(7)]);
%fin
crdiff=-finaccomb+finmucomb;
bar(10,mean(crdiff),0.6,'FaceColor',colvl,'EdgeColor',colvl)
stvl=std(crdiff)./sqrt(length(crdiff));
plot([10 10],[mean(crdiff)-stvl mean(crdiff)+stvl],'k');
% plot([13 13],[mean(finaccomb)-stvl mean(finaccomb)+stvl],'k');
% plot(13*ones(length(finaccomb),1),finaccomb,'ko');
% 
% bar(14,mean(finmucomb),0.6,'FaceColor','none','EdgeColor','r')
% stvl=std(finmucomb)./sqrt(length(conmucomb));
% plot([14 14],[mean(finmucomb)-stvl mean(finmucomb)+stvl],'r');
% plot(14*ones(length(finmucomb),1),finmucomb,'ro');
% 

% daspect([1 40 1])
%  axis([-.5 15 -1 1])
% 

