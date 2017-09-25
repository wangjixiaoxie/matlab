function lv_labelpsds(mdlfile,batchname)

load(mdlfile) %have: smallpcmdl, basismtx,newclusterlabels,normby
pcmdl = smallpcmdl;

load(['psddata_5parts_' batchname '.mat'])


%% identify and take out WN syllables
maxvol = cellfun(@(x) mean(x.^2),SylOut);
wnidx = maxvol>6000000;

oldsylout = SylOut;
oldpsdsylout = PSDSylOut;
SylOut = SylOut(~wnidx);
PSDSylOut = PSDSylOut(~wnidx,:);



%% calc distance matrix to basis
 allpsds = PSDSylOut;

distmtx = lv_calcdistmatrix(basismtx,allpsds,normby);

%% cluster
[clusnum,~,p,~,m] = cluster(pcmdl,distmtx);

% allposterior = posterior(pcmdl,distmtx);
% [bestpost, clusnum2] = max(allposterior,[],2); %clusnum ist gleich wie out
% clusnum==clusnum2

unique(clusnum)

% p

%% look at
[sortclus, sortidx] = sort(clusnum);
figure
imagesc(allpsds(sortidx,:))
hold on
trans = find(diff(sortclus));
for i =1:length(trans)
    line(get(gca,'xlim'),[trans(i) trans(i)],'color','r')
    text(1,trans(i)-20,newclusterlabels(sortclus(trans(i))),'color','r','fontsize',12)
    
end

%% look at:
newclusnum = clusnum;
% % 
[~,plotbasis] = princomp(allpsds);
plotbasis = plotbasis(:,1:3);

figure
hold on
dim1 = 1;
dim2 = 2;
dim3 = 3;
idx = newclusnum==1;
plot3(plotbasis(idx,dim1),plotbasis(idx,dim2),plotbasis(idx,dim3),'k+')
idx = newclusnum==2;
plot3(plotbasis(idx,dim1),plotbasis(idx,dim2),plotbasis(idx,dim3),'bo')
idx = newclusnum==3;
plot3(plotbasis(idx,dim1),plotbasis(idx,dim2),plotbasis(idx,dim3),'ro')
idx = newclusnum==4;
plot3(plotbasis(idx,dim1),plotbasis(idx,dim2),plotbasis(idx,dim3),'go')
idx = newclusnum==5;
plot3(plotbasis(idx,dim1),plotbasis(idx,dim2),plotbasis(idx,dim3),'yo')
idx = newclusnum==6;
plot3(plotbasis(idx,dim1),plotbasis(idx,dim2),plotbasis(idx,dim3),'mo')
idx = newclusnum==7;
plot3(plotbasis(idx,dim1),plotbasis(idx,dim2),plotbasis(idx,dim3),'co')
idx = newclusnum==8;
plot3(plotbasis(idx,dim1),plotbasis(idx,dim2),plotbasis(idx,dim3),'o','color',[.3 .2 .7])
idx = newclusnum==9;
plot3(plotbasis(idx,dim1),plotbasis(idx,dim2),plotbasis(idx,dim3),'m*')
idx = newclusnum==10;
plot3(plotbasis(idx,dim1),plotbasis(idx,dim2),plotbasis(idx,dim3),'c*')

legend(newclusterlabels)

% plot3(oldbasismtx(wnidx,dim1),oldbasismtx(wnidx,dim2),oldbasismtx(wnidx,dim3),'k*')
% legend([newclusterlabels 'w'])

%% add WN syllables back in

ucl = unique(newclusnum);
wnclusnum = ucl(end)+1;

newclusnumwithWN = nan(size(oldsylout));
newclusnumwithWN(~wnidx) = newclusnum;
newclusnumwithWN(wnidx) = wnclusnum;

newclusterlabels = [newclusterlabels 'w'];
newclusnum = newclusnumwithWN;

%% label
lv_givelabels(batchname,newclusterlabels,newclusnum)

