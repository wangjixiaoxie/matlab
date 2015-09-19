%list of stim times

for bsvl=1:length(sumbs)
    crbs=sumbs(bsvl);
    bsoutsh(bsvl).early=[];
    bsoutsh(bsvl).late=[];

    for shiftnum=1:length(crbs.shiftruns)
        shiftinds=crbs.shiftruns{shiftnum};
        
        ind_early=find(mod(crbs.tmvec(shiftinds,1),1)<(10/24));
        ind_late=find(mod(crbs.tmvec(shiftinds,2),1)>(17/24));
        
        bsoutsh(bsvl).early=[bsout(bsvl).early shiftinds(ind_early)];
        bsoutsh(bsvl).late=[bsoutsh(bsvl).late shiftinds(ind_late)];
    
     
    end
    
end

for bsvl=1:length(phsumbs)
    crbs=phsumbs(bsvl);
    phbsout(bsvl).early=[];
    phbsout(bsvl).late=[];
    for shiftnum=1:length(crbs.shiftruns)
        shiftinds=crbs.shiftruns{shiftnum};
        ind_early=find(mod(crbs.rawtimes(shiftinds,1),1)<(10/24));
        ind_late=find(mod(crbs.rawtimes(shiftinds,2),1)>(17/24));
        phbsout(bsvl).early=[phbsout(bsvl).early shiftinds(ind_early)];
        phbsout(bsvl).late=[phbsout(bsvl).late shiftinds(ind_late)];
    end
    
end