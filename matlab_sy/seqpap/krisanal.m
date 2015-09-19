%quick script to analyze kris' questions about syllable gaps.

sylls=[1 2]
nsyll=2;
nfv=6;
hstbnds=100:5:200

%for each fvstruct;
% fvanal{1}=fvcombpre;
% fvanal{2}=fv26
% fvanal{3}=fv27
% fvanal{4}=fv30
% fvanal{5}=fv01;
% fvanal{6}=fv02;







%output is 
%dfot{fvnum}{syllnum}
for ii=1:nfv
    for jj= 1:nsyll
    otdiff{ii}{jj}=[]
    end
end
    

for fvnum=1:length(fvanal)
    crfvst=fvanal{fvnum};
    for crsyll=1:length(sylls)
       crfv=crfvst{crsyll}
       for ii=1:length(crfv)
           crons=crfv(ii).ons;
           croffs=crfv(ii).offs;
           crind=crfv(ii).ind;
           crdiff=crons(crind)-crons(crind-1);
           
           otdiff{fvnum}{crsyll}=[otdiff{fvnum}{crsyll} crdiff];
           
       end
        
    end
  
end


for fvnum=1:length(fvanal)
    for crsyll=1:length(sylls)
        
        %exclude outliers
        inds=find(otdiff{fvnum}{crsyll}>hstbnds(1)&otdiff{fvnum}{crsyll}>hstbnds(2))
        
        mnout{fvnum}{crsyll}=mean(otdiff{fvnum}{crsyll}(inds));
        hstout{fvnum}{crsyll}=histc(otdiff{fvnum}{crsyll},hstbnds);
        sterout{fvnum}{crsyll}=std(otdiff{fvnum}{crsyll}(inds))./sqrt(length(inds));
    end
end

%make histograms
% e.g.  hstb1=histc(otdiff{2}{1},hstbnds)

mo=mnout
so=sterout

figure;

for ii=1:length(fvanal)

subplot(length(fvanal),1,ii)
stairs(hstbnds,hstout{ii}{1}./sum(hstout{ii}{1}),'r');
hold on;
stairs(hstbnds,hstout{ii}{2}./sum(hstout{ii}{2}),'k');
f=ii;s=1;
plot([mo{f}{s}-so{f}{s} mo{f}{s}+so{f}{s}],[.6 .6],'r')
f=ii;s=2;
plot([mo{f}{s}-so{f}{s} mo{f}{s}+so{f}{s}],[.7 .7],'k')

%base lines
plot([mo{1}{1} mo{1}{1}],[0 1],'r--');
plot([mo{1}{2} mo{1}{2}],[0 1],'k--');

axis square
end
