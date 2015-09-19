mnsum=[]
stsum=[]
daysum=[]
flrvls=floor(vals(:,1));
unqvls=unique(flrvls);
for ii=1:length(unqvls);
    vlsinds=find(flrvls==unqvls(ii));
    mnsum=[mnsum mean(vals(vlsinds,2))]
    stsum=[stsum std(vals(vlsinds,2))]
    daysum=[daysum vals(vlsinds(1),1)]
    
end
    