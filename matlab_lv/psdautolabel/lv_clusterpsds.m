close all
%lena modified from dave/hamish

npsds = 5;
nsyls = 3000;
batchname = 'batch.keep.lvrand.keep';
% batchname = 'batch.handselect';
load(['psddata_' num2str(npsds) 'parts_' batchname '.mat'])

PSDSylOut = PSDSylOut(1:nsyls,:); %lena hier limit syls; % +51


myfs = Fs./1000;
allpsds = PSDSylOut;

% Parameters
nlower = 400;
nupper = 3500; %size(allpsds,2)./npsds;
ncomponents = 15;
fuzzyline = 0.95; %minimum posterior probability of belonging to cluster in order to stay in cluster
fixeddurcutoff = 25;%15; %which window around median syllable duration is used to calculate std of syllable duration
sdcutoff = 5; %how far out in each dimension to be removed as an outlier?



% 1. Remove lowest freq: (war 200:2:8000)
% size for 4: 15604
psdsize = size(allpsds,2)./npsds;
% toremove = [];
% for p = 1:npsds
%     toremove = [toremove (p-1)*psdsize+1:(p-1)*psdsize+nlower];
%     toremove = [toremove (p-1)*psdsize+nupper:(p-1)*psdsize+psdsize];
% end
% allpsds(:,toremove) = [];


%dave method
% basisidx = 13:127;
% basisidx = [12 16 21 29 38 46 49 53 57 59 63 66 68 70 79 82 86 104 116 120]; %br br br43
basisidx =  [4 7 8 14 15 16 23 24 25 26 27 34 35 36 37 38 39 46 47 50]; %br43 new for 0528
% basisidx = [1 3 7 10 11 14 15 16 32 33 34 35 37 38 40 42 44 47 53 54]; % rd66?
% basisidx = [4 6 11 13 14 15 16 17 18 19 20 25 26 28 29 30 31 32 33 34];% rd49
% basisidx = 1:20; %bu50 %wh25 wh09 %
% basisidx = [12 13 21 23 24:39]; %rd82


basismtx = allpsds(basisidx,:);

[distmtx, normby] = lv_calcdistmatrix(basismtx,allpsds,0);

        
fprintf('done \n')
%% identify and take out WN syllables - do before calc distance matrix!!!
% maxvol = cellfun(@(x) max(x),SylOut);
% maxvol = maxvol(1:nsyls+51);
% wnidx = maxvol>4000;
% 
% oldbasismtx = basismtx;
% basismtx = basismtx(~wnidx,:);
% oldsylout = SylOut;
% oldpsdsylout = PSDSylOut;
% SylOut = SylOut(~wnidx);
% PSDSylOut = PSDSylOut(~wnidx,:);


%% 3 fold cross validation for number of components:
clear BIC_notcv BIC_cv nloglik

cvidx = randperm(nsyls);
nfolds = 3;
nfoldsyls = nsyls/nfolds;
ndims = size(distmtx,2);
for f = 1:nfolds
    testidx = (f-1)*nfoldsyls+1:f*nfoldsyls;
    testidx(1) 
    testidx(end)
    testdata = distmtx(cvidx(testidx),:);
    traindata = distmtx;
    traindata(cvidx(testidx),:) = [];

    for k = 1:15
        options = statset('MaxIter',5000,'Display','final');

        pcmdl=fitgmdist(traindata,k,'CovarianceType','Full', 'Replicates',25,'SharedCovariance',false,'Start','plus','Option',options);
        BIC_notcv(k,f) = pcmdl.BIC;
        [~,nloglik(k,f)] = cluster(pcmdl,testdata);
        df = (ndims*ndims - ndims)/2 + ndims+ ndims +1;
        mparams = df*k - 1;
        BIC_cv(k,f) = 2*nloglik(k,f) + mparams*log(nfoldsyls);
    end
end
figure
plot(BIC_cv)
[~,lala] = min(BIC_cv,[],1);
lala
[~,numclust] = min(mean(BIC_cv,2));
%% 5. cluster
%often works better with 1-2 clusters more than determined in previous step
options = statset('MaxIter',5000,'Display','final');

i = numclust; 

pcmdl = fitgmdist(distmtx,i,'CovarianceType','Full', 'Replicates',25,'SharedCovariance',false,'Start','plus','Option',options);
% BICOut(i) = pcmdl.BIC;
% pcmdl.BIC


%% eval model

allposterior = posterior(pcmdl,distmtx);
[bestpost, clusnum] = max(allposterior,[],2); 

figure
plot(sort(bestpost))
title('posterior probability of best cluster')


%plot sorted PSDs
[sortclus, sortidx] = sort(clusnum);
figure
imagesc(PSDSylOut(sortidx,:))
hold on
trans = find(diff(sortclus));
for i =1:length(trans)
    line(get(gca,'xlim'),[trans(i) trans(i)],'color','r')
    text(1,trans(i)-20,num2str(sortclus(trans(i))),'color','r','fontsize',12)
end


%% look at and label cluster
% adjust parameters for removal of outliers
% do not close figures after labeling
% 


close all

distbasis = distmtx;

fuzzysylls = bestpost<fuzzyline;
mean(fuzzysylls)

figure
hist(clusnum,0:numclust)
clusterlabel = [];

m = mean(distbasis,1);
for i = 1:numclust
    fh(i) = figure('position',[50 50 1500 1000]);
    set(gcf,'UserData',{clusterlabel})
    set(gcf,'KeypressFcn',@assignlabel)
    hold on

    idx = clusnum==i&~fuzzysylls;

    clusnum(fuzzysylls) = 0;
    
    %within vs between cluster similarity
    cisize = sum(idx);
    m_i = mean(distbasis(idx,:),1);
    thisclustwss = sum((distbasis(idx,:)-repmat(m_i,cisize,1)).^2,1);
    thisclustbss = cisize*(m - m_i).^2;
%     fprintf('cluster %d: WSS/BSS: %.2f \n',i,sum(thisclustwss)./sum(thisclustbss))

    ahandle = subplot(4,4,1);
    hold on
    text(0.1,0.3,['n not fuzzy= ' num2str(cisize)])
    text(0.1,0.1,['Proportion: ' num2str(pcmdl.ComponentProportion(i).*100)])
    text(0.1,0.2,['WSS/BSS: ' num2str(sum(thisclustwss)./sum(thisclustbss))])
    text(0.1,0.7,num2str(i),'fontsize',14)

    subplot(4,4,2)
    plot(sort(bestpost(idx)))
    
    subplot(4,4,3)
    hold on
    plot(thisclustwss,'r')
    plot(thisclustbss,'b')
    
    subplot(4,4,4)
    plot(mean(allposterior(clusnum==i,:)))
    

    sylidx = find(idx);
    try
        allsyls = cellfun(@(x) x(500:end-500),SylOut(sylidx(1:50)),'UniformOutput',false);
    catch
        allsyls = cellfun(@(x) x(500:end-500),SylOut(sylidx),'UniformOutput',false);
    end
    SylSpectData = cell2mat(allsyls);
    [~,F,T,P] = spectrogram(SylSpectData,512,256,(200:2:8000),Fs);
    subplot(4,4,6:7)
    T=T*1000;
    imagesc(T,F,log10(P))
    set(gca,'clim',[-2 3])
    axis xy

    
%     subplot(4,4,10:11)
%     plot(allpsds(idx,:)')
    
    subplot(4,4,8)
    syldur = cellfun(@(x) size(x,2),SylOut(idx))./myfs;
    hist(syldur,1000)
    line([median(syldur) median(syldur)],get(gca,'ylim'),'color','r')
    thisstd = std(syldur(syldur>median(syldur)-fixeddurcutoff&syldur<median(syldur)+fixeddurcutoff));
    durcutoff = 2.2*thisstd
    line([median(syldur)-durcutoff median(syldur)-durcutoff],get(gca,'ylim'),'color','r')
    line([median(syldur)+durcutoff median(syldur)+durcutoff],get(gca,'ylim'),'color','r')

    %remove duration outliers
    removeidx = (syldur>median(syldur)+durcutoff | syldur<median(syldur)-durcutoff);
    clusnum(sylidx(removeidx)) = 0;
    nremoved = sum(removeidx);

       
    %remove outliers in first 20 dimensions
    for k = 1:size(basismtx,1)
    removeidx = basismtx(k,idx)>median(basismtx(k,idx))+sdcutoff*std(basismtx(k,idx)) | basismtx(k,idx)<median(basismtx(k,idx))-sdcutoff*std(basismtx(k,idx));
    stillincluster = sum(clusnum(sylidx(removeidx))== i);
    clusnum(sylidx(removeidx)) = 0;
    nremoved = nremoved+stillincluster;    
    end

    subplot(4,4,1)
    text(0.1,0.5,['n removed = ' num2str(nremoved)])
    
    subplot(4,4,14:15)
    
    cleansylidx = find(clusnum==i);
    try
        allsyls = cellfun(@(x) x(500:end-500),SylOut(cleansylidx(1:50)),'UniformOutput',false);
    catch
        allsyls = cellfun(@(x) x(500:end-500),SylOut(cleansylidx),'UniformOutput',false);
    end
    SylSpectData = cell2mat(allsyls);
    [~,F,T,P] = spectrogram(SylSpectData,512,256,(200:2:8000),Fs);
    T=T*1000;
    imagesc(T,F,log10(P))
    set(gca,'clim',[-2 3])
    axis xy
    

end

%% collect labels

clusterlabels = {};
for i = 1:numclust
   ud = get(fh(i),'UserData');
    clusterlabels{i} = ud{1};
end
figure
hist(clusnum,0:numclust)

%%  clusters label / merge clusters with same label:
close all
ucl = unique(clusterlabels);
nclus = length(ucl);
newclusterlabels = {};
for cc = 1:nclus
    %find clusters with this label
    cidx = find(strcmp(clusterlabels,ucl(cc)));
    for fc = 1:length(cidx)
        newclusnum(clusnum==cidx(fc)) = cc;
    end
    newclusterlabels(cc) = ucl(cc);
end
zeroidx = find(strcmp(newclusterlabels,'0'));
if ~isempty(zeroidx)
    newclusnum(clusnum==0) = zeroidx;
else
    ucl = [ucl '0'];
    newclusterlabels = [newclusterlabels '0'];
    zeroidx = nclus+1;
    newclusnum(clusnum==0) = zeroidx;
    
end



%check again
[sortclus, sortidx] = sort(newclusnum);

figure
imagesc(PSDSylOut(sortidx,:))
hold on
trans = find(diff(sortclus));
for i =1:length(trans)
    line(get(gca,'xlim'),[trans(i) trans(i)],'color','r')
    text(1,trans(i)-20,newclusterlabels{i},'color','r','fontsize',12)
end
text(1,trans(i)+40,newclusterlabels{i+1},'color','r','fontsize',12)


zeroidx = find(strcmp(newclusterlabels,'0'));
unassignedidx = newclusnum==zeroidx;
fprintf('unassigned proportion: %.2f \n',sum(unassignedidx)./length(newclusnum))

figure
hist(newclusnum,0:nclus)

%% look at assigned
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

% plot3(oldbasismtx(wnidx,dim1),oldbasismtx(wnidx,dim2),oldbasismtx(wnidx,dim3),'k*')
legend([newclusterlabels 'w'])



%% give labels
lv_givelabels(batchname,newclusterlabels,newclusnum)


%% make small model with actual number of clusters
clusterlabelstring = [newclusterlabels{:}];
clusterstokeep = unique(newclusnum);
clusterstokeep = clusterstokeep(clusterlabelstring~='0');

%take out clusters labeled 0
zeroclusters = find(clusterlabelstring == '0');
keepidx = ones(size(newclusnum));
for i = 1:length(zeroclusters)
    keepidx(newclusnum==zeroclusters(i)) = 0;
end
keepidx = keepidx>0;
assert(all(cellfun(@(x) strcmp(x,'0'),newclusterlabels(newclusnum(keepidx==0)))))

%calc parameters of new model (could be done directly from labels if have
%enough labeled songs...
newmu = grpstats(distmtx(keepidx,:),newclusnum(keepidx));
newsigma = grpstats(distmtx(keepidx,:),newclusnum(keepidx),'cov');
newsigma = shiftdim(newsigma,1);
newnclusters = length(clusterstokeep);
cp = histc(newclusnum(keepidx),clusterstokeep);
cp = cp./sum(cp);
starti = struct('mu',newmu,'Sigma',newsigma,'ComponentProportion',cp);

% smallpcmdl=fitgmdist(basismtx,newnclusters,'CovarianceType','Full','RegularizationValue',0.01,'SharedCovariance',false,'Option',options,'Start',starti);
% smallbic = smallpcmdl.BIC;

smallpcmdl = gmdistribution(newmu,newsigma,cp);

newclusterlabels = newclusterlabels(clusterstokeep);

save('smallmdl7','newclusterlabels','smallpcmdl','basismtx','normby')

%% check all labels

% jc_chcklbl('batch.handselect','d',.01,.01,'','',[],[])
% jc_chcklbl('batch.handselect','e',.01,.01,'','',[],[])
% jc_chcklbl('batch.handselect','s',.01,.01,'','',[],[])
% jc_chcklbl('batch.handselect','r',.01,.01,'','',[],[])
% jc_chcklbl('batch.handselect','p',.01,.01,'','',[],[])
% jc_chcklbl('batch.handselect','u',.01,.01,'','',[],[])
% 
% jc_chcklbl('batch.handselect','0',.01,.01,'','',[],[])