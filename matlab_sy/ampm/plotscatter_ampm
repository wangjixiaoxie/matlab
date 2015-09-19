%ampmscatter.m
ps.asympcol='k';
ps.shiftcol='r';

crbs=phsumbs(5);
%find the am_pairs go that are asympruns
%plot some value
crasympruns=crbs.allasympruns;
ampmplotind.asymp=[];
ampmplotind.shift=[];
ammuplotind.asymp=[];
ammuplotind.shift=[];


for ii=1:length(crasympruns)
    crasympvls=crasympruns{ii};
    plotind=find(ismember(am_pmout.amruns,crasympvls))
    ampmplotind.asymp=[ampmplotind.asymp plotind']
    plotind=find(ismember(am_muout.amruns,crasympvls));
    ammuplotind.asymp=[ammuplotind.asymp plotind']
    
end


%find the am_pm pairs that are shiftruns
%plot through in one color.
crshiftruns=crbs.allshiftruns;

for ii=1:length(crshiftruns)
    crshiftvls=crshiftruns{ii};
    plotind=find(ismember(am_pmout.amruns,crshiftvls))
    ampmplotind.shift=[ampmplotind.shift plotind']
    plotind=find(ismember(am_muout.amruns,crshiftvls));
    ammuplotind.shift=[ammuplotind.shift plotind']
    
end

