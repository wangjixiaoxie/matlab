clear pathvl vals
pathvl{1}='/doyale1/twarren/w23pk100/directed/'
catchvl{1}='batch'
pathvl{2}='/doyale1/twarren/w23pk100/initscreen/'
catchvl{2}='batchcomb.keep.rand'

numnotes=2
pathvl{3}='/doyale1/twarren/w23pk100/ac91607/'
pathvl{4}='/doyale1/twarren/w23pk100/mu91607/'
% pathvl{2}='/doyale1/twarren/w23pk100/mu91707/'
pathvl{5}='/doyale1/twarren/w23pk100/ac91707/'
pathvl{6}='/doyale1/twarren/w23pk100/mu91707/'
pathvl{7}='/doyale1/twarren/w23pk100/mu918/'
pathvl{8}='/doyale1/twarren/w23pk100/ac918/'
pathvl{9}='/doyale1/twarren/w23pk100/mu92005/'

catchvl{3}='batch.keep.rand'
catchvl{5}='batchcomb.rand'
catchvl{4}='batch.keep'
catchvl{6}='batch.keep'
catchvl{7}='batch.keep'
catchvl{8}=catchvl{7}
catchvl{9}='batch.keep.rand'

muon=[1 4 6 7 9]
acon=[2 3 5 8]

%here is the 
 valsemu=[]
valseac=[]
valsbmu=[]
valsbac=[]




for ii=1:length(muon)
    ind=muon(ii)
    valsemu=[valsemu ;vals{2*ind-1}]
end

for ii=1:length(acon)
    ind=acon(ii)
    valseac=[valseac; vals{2*ind-1}]
end

for ii=1:length(muon)
    ind=muon(ii)
    valsbmu=[valsbmu ;vals{2*ind}]
end

for ii=1:length(acon)
    ind=acon(ii)
    valsbac=[valsbac; vals{2*ind}]
end

figure
hold on;
plot (valsbac(:,1),valsbac(:,2),'k.','MarkerSize',5.3)
hold on;
plot (valsbmu(:,1),valsbmu(:,2),'r.','MarkerSize',5.3)


figure
hold on;
plot (valseac(:,1),valseac(:,2),'k.','MarkerSize',5.3)
hold on;
plot (valsemu(:,1),valsemu(:,2),'r.','MarkerSize',5.3)



for ii=1:length(pathvl)
    for jj=1:numnotes

        
          strcmd=['load ' pathvl{ii} catchvl{ii} '.mat']
        eval(strcmd);
        
        vals{(ii*2)+jj-2}=getvals(fv{jj},1,'TRIG');
    end
end

edges{1}=[3000:50:3700]
edges{2}=[4200:50:5000]

histbdir=histc(vals{2}(:,2),edges{2})
histbun=histc(vals{4}(:,2),edges{2})
histbdir=normhist(histbdir)
histbun=normhist(histbun)

histedir=histc(vals{1}(:,2),edges{1})
histeun=histc(vals{3}(:,2),edges{1})
histedir=normhist(histedir)
histeun=normhist(histeun)





histedir=histc/sum(histbpre)
histbmu=histbmu/sum(histbmu)
histbpost=histbpost/sum(histbpost)

figure
stairs(edges{2},histbdir)
hold on;stairs(edges{2},histbun,'r--')

figure
stairs(edges{1},histedir)
hold on;stairs(edges{1},histeun,'r--')


figure
stairs(edges,histbpost)
hold on;stairs(edges,histbmu,'r--')


indb=find(valsb(:,2)>3400&valsb(:,2)<4200)
stdb=std(valsb(indb,2))
indmb=find(valsmb(:,2)>3400&valsmb(:,2)<4200)
stdmb=std(valsmb(indmb,2))
