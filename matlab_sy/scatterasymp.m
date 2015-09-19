%7.15 analysis.
diffcomb=[];
offcomb=[];

for ii=1:length(sumdyn)
    crsd=sumdyn(ii)
    if(length(crsd.tms)>1)
        for jj=2:length(crsd.tms)
            cracdiff=crsd.acasympdist(jj-1)-crsd.acasympdist(jj)
            croff=crsd.off(jj);
            if(crsd.acz(jj)>0)
                croff=-croff
            end
            diffcomb=[diffcomb cracdiff];
            offcomb=[offcomb croff];
        end
        
    end
end