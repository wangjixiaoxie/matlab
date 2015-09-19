
%run /oriole5/evrenlesiondata/mkfigs.m to get prevls, postvls
%typically saved in /oriole
function[outpreind,outpstind] =lesionplot_tmcourse(prevls,postvls)

clear indct mxindprelesion
%tstlesionplotscript
specialcolorinds=[2]
subplot(221)
for ii=1:length(prevls.postmeanmat(:,1));
    ind=find(prevls.postmeanmat(ii,:)~=0);
    recind=find(prevls.recmeanmat(ii,:)~=0);
    preind=find(prevls.premeanmat(ii,:)~=0);
    
    offsetcr=(prevls.postmeanmat(ii,ind));
    offsetrec=(prevls.recmeanmat(ii,recind));
    
    initvl=prevls.trigmeanmat(ii,1);
    if initvl==0
        initvl=prevls.trigmeanmat(ii,2);
    end
    offset_trig=abs((prevls.trigmeanmat(ii,ind)-initvl));
    if(~ismember(ii,specialcolorinds))
        plot(ind,offsetcr,'ko');
        hold on;
        plot(recind+ind(end),offsetrec,'ko');
        subplot(222)
        plot(ind,offset_trig,'r')
        hold on;
        subplot(221)
        plot(ind,offsetcr,'k');
        plot(recind+ind(end),offsetrec,'k');
    else
        plot(ind,offsetcr,'co');
        hold on;
         plot(recind+ind(end),offsetrec,'co');
        subplot(222)
        
        plot(ind,offset_trig,'m');
        hold on;
        subplot(221)
        plot(ind,offsetcr,'c');
         plot(recind+ind(end),offsetrec,'c');
    end
    hold on;
    mxindprelesion(ii)=max(ind);
end
subplot(224);
%first find values of mxind equal to or above a threshold
PRETHRESH=4;
threshruns=find(mxindprelesion>=PRETHRESH);
[outmn,outstd,pvlpre]=calcmeanstder2(prevls.postmeanmat(threshruns,:))
threshruns=threshruns(1:8)
outpreind=threshruns;
for indtmp=1:length(threshruns)
    crindtmp=threshruns(indtmp);
    nozeroind=find(prevls.premeanmat(crindtmp,:)~=0)
    prevls.premeanmat(crindtmp,nozeroind)=prevls.premeanmat(crindtmp,nozeroind)-prevls.premean(crindtmp);
end
[outmnpre,outstdpre]=calcmeanstder2(prevls.premeanmat(threshruns,:));
% plot(1:PRETHRESH,outmn(1:PRETHRESH),'ko');
hold on;
plot([0 1:PRETHRESH],[0 outmn(1:PRETHRESH)],'k');
plot([1:PRETHRESH;1:PRETHRESH],[outmn(1:PRETHRESH)-outstd(1:PRETHRESH);outmn(1:PRETHRESH)+outstd(1:PRETHRESH)],'k')
text(5,100,['n=' num2str(length(threshruns))]) 
plot(-1:0,outmnpre(1:2),'k');
plot([-1:0;-1:0],[outmnpre(1:2)-outstdpre(1:2);outmnpre(1:2)+outstdpre(1:2)],'k');
clear mxindprelesion

x=[0 2.5 4.5]
y1(1)=0;
y1(2)=mean(outmn(2:3))*.6277;
y1(3)=mean(outmn(4))*.766;

plot(x,y1,'r--','Marker','o');


subplot(223)
for ii=1:length(postvls.postmeanmat(:,1));
    ind=find(postvls.postmeanmat(ii,:)~=0);
    
    preind=find(postvls.premeanmat(ii,:)~=0);
    
    offsetcr=postvls.postmeanmat(ii,ind);
    plot(ind,offsetcr,'ko')
    hold on;
    plot(ind,offsetcr,'k')
    mxindprelesion(ii)=max(ind);
end

subplot(224);
%first find values of mxind equal to or above a threshold
PRETHRESH=6;
threshruns=find(mxindprelesion>=PRETHRESH);
[outmn,outstd,pvlpst]=calcmeanstder2(postvls.postmeanmat(threshruns,:))
%excluding two runs where there was no out 6-2, 9-2, and excluding 5, where
%prelesion shift was only for three days.
indvl=[1:8 10 12:14 ]
outpstind=indvl;
for indtmp=indvl
    crindtmp=threshruns(indtmp);
    nozeroind=find(postvls.premeanmat(crindtmp,:)~=0)
    postvls.premeanmat(crindtmp,nozeroind)=postvls.premeanmat(crindtmp,nozeroind)-postvls.premean(crindtmp);
end
% plot(1:PRETHRESH,outmn(1:PRETHRESH),'ko');
hold on;
plot([0 1:PRETHRESH],[0 outmn(1:PRETHRESH)],'r');
plot([1:PRETHRESH;1:PRETHRESH],[outmn(1:PRETHRESH)-outstd(1:PRETHRESH);outmn(1:PRETHRESH)+outstd(1:PRETHRESH)],'r')

ind=find(mean(postvls.premeanmat,2)>1);
for kk=1:length(ind)
    crind=ind(kk);
    postvls.premeanmat(crind,:)=postvls.premeanmat(crind,:)-mean(postvls.premeanmat(crind,:))
end


[outmnpre,outstdpre]=calcmeanstder2(postvls.premeanmat(threshruns,:));
plot(-1:0,outmnpre(1:2),'r');
plot([-1:0;-1:0],[outmnpre(1:2)-outstdpre(1:2);outmnpre(1:2)+outstdpre(1:2)],'r');




text(5,50,['n=' num2str(length(threshruns(indvl)))],'Color','r') 

