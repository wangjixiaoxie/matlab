%repeat_stabil_anal - combining JC 4 day and 30 day data source.

figure;
cd /cardinal9/Repeat_Analysis/
load stationaryFinal.mat;
load stationaryMONTH.mat
subplot(4,4,9:10)
plot(avgrplength(:,1),avgrplength(:,4),'ko')
axis square
axis([0 12 0 12])
hold on;
plot([0 12], [0 12],'k--')
hold on;
for ii=1:length(Rep)
    plot(median(Rep(ii).rplengthPre), median(Rep(ii).rplengthPost),'ro')
    
end

%THE COMMENTED BIT BELOW WILL PLOT SCATTER OF TRANSITION PROBABILITIES.
%NEED TO CALCULATE "probstay", which is done in 
load /cardinal9/SyntaxLearningFigs/Repeats5.mat
[probstay,numstay]=repeatsTOptrans(Exp);


subplot(4,4,11:12)

NUMTOKEEP=15
nn=9
             for i=1:nn
                 for j=1:length(probstay(i).day(:,1))
                     % - note that I only keep the ones with a decently large sample size
                     if (numstay(i).day(j,1)>NUM2KEEP&numstay(i).day(j,4)>NUM2KEEP)
                        plot(probstay(i).day(j,1),probstay(i).day(j,4),'ko')
                        hold on;
                     end
                 end
             end
             
             axis square;
             
             axis([0 1 0 1])
             plot([0 1], [0 1],'k--')
             
             
             for ii=1:length(Rep)
                 inds=find(Rep(ii).num_in_ptransPre>15&Rep(ii).num_in_ptransPost>15)
                    plot(Rep(ii).ptransPre(inds),Rep(ii).ptransPost(inds),'ro')
             end
 
             
subplot(4,4,1:4)

         
    ax=gca();

    exsong.ax=ax;
    exsong.path='/cardinal9/Repeat_Analysis/pk30bk79/ampoff'
    exsong.fn='pk30bk79_030108_1236.5025.cbin';
    exsong.bnds=[6.5 9]
    plotcbin(exsong);    

  %get example data
  for ii=1:4
    pt{ii}='/cardinal9/Repeat_Analysis/pk30bk79/screen'
  end
  for ii=5:8
      pt{ii}='/cardinal9/Repeat_Analysis/pk30bk79/ampoff'
  end
  bt{1}='bat09'
  bt{2}='batch10'
  bt{3}='batch11'
  bt{4}='batch12'
  bt{5}='batch31.rand'
  bt{6}='batch01.rand'
  bt{7}='batch02.rand'
  bt{8}='batch03'
  
   
  subplot(4,4,5)
  for ii=1:8
      cmd=['cd ' pt{ii}];
      eval(cmd);
      cmd=['load ' bt{ii} '.mat'];
      eval(cmd);
      ind=find(distnt>2);
      cr_rps{ii}=distnt(ind);
      for kk=1:500 % resampling to obtain CIs
             nrsamp=ind(ceil(length(ind)*rand(1,length(ind))));
             indcheck=find(nrsamp>length(ind));
             if(~isempty(indcheck))
                 nrsamp(indcheck)=length(ind);
             end
             msmn(kk)=mean(cr_rps{ii}(nrsamp));
       end
                            prchigh(ii)=prctile(msmn,97.5); % 95prct
                            prclow(ii)=prctile(msmn,2.5);  
                            prcmn(ii)=prctile(msmn,50);
  end
  
  hstinds=[1 4 8]
  edges=1:12
  col={'k' 'r' 'c'}
  for ii=1:length(hstinds)
     hstout=histc(cr_rps{hstinds(ii)},edges);
     stairs(edges,hstout./sum(hstout),col{ii})
     hold on;
     crmn=prcmn(hstinds(ii));
     plot(crmn,.3+ii*.03,'Marker','^','Color',col{ii});
  end


%PLOT OF EXAMPLE OF STABILITY ACROSS DAYS

  ax(1)=subplot(4,4,6:8)
  
  dayind{1}=[1:4]
  dayind{2}=[5:8]
  dys=[1 2 3 4 25 26 27 28]
  for crdy=1:2
        dysind=dayind{crdy}
       
        lndys=length(dysind)
        plot(dys(dysind),prcmn(dysind),'k');
        hold on;
        plot([dys(dysind) ;dys(dysind)],[prclow(dysind); prchigh(dysind)],'k');
        hold on;
        plot(dys(dysind),prcmn(dysind),'ko');
    end



