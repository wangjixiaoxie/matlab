function [tm,mnpt,stdpt]=plotdailymeanstdv(avls,CRNT)
combvls=[];
mnpt=[]
stdpt=[]
tm=[];

%create a master structure of combvls
for ii=1:length(avls.PTIND)
    crind=avls.PTIND(ii);
    crpt=avls.pvls{crind}
    crbt=avls.cvl{crind}
     if(isfield(avls,'baspath'))
        cmd=['cd ' avls.baspath crpt]
    else
        cmd=['cd ' crpt]
    end
    eval (cmd);  
    
    cmd =['load ' crbt '.mat']
    eval(cmd);
    combvls=[combvls; vls{CRNT}]
    
end

%now go through the unique days
flrvls=floor(combvls(:,1));
unqvls=unique(flrvls);
for ii=1:length(unqvls)
    crind=find(flrvls==unqvls(ii))
    mnpt=[mnpt median(combvls(crind,2))];
    tm=[tm median(combvls(crind,1))];
    stdpt=[stdpt std(combvls(crind,2))];
end

figure
plot(tm,mnpt,'ko');
hold on;
plot([tm;tm],[mnpt-stdpt;mnpt+stdpt],'k')