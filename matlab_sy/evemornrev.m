%list of stim times

for bsvl=1:length(sumbs)
    crbs=sumbs(bsvl);
    bsoutshrev(bsvl).early=[];
    bsoutshrev(bsvl).late=[];
%     bsoutshrev(bsvl).early=[];
%     bsoutshrev(bsvl).late
    for shiftnum=1:length(crbs.revruns)
        shiftinds=crbs.revruns{shiftnum};
        
        ind_early=find(mod(crbs.tmvec(shiftinds,1),1)<(11/24));
        ind_late=find(mod(crbs.tmvec(shiftinds,2),1)>(15/24));
        if(~isempty(ind_early))
            bsoutshrev(bsvl).early=[bsoutshrev(bsvl).early makerow(shiftinds(ind_early))];
        end
        if(~isempty(ind_late))
            bsoutshrev(bsvl).late=[bsoutshrev(bsvl).late makerow(shiftinds(ind_late))];
        end
     
    end
    
end

for bsvl=1:length(phsumbs)
    crbs=phsumbs(bsvl);
    phbsoutrev(bsvl).early=[];
    phbsoutrev(bsvl).late=[];
    for shiftnum=1:length(crbs.revruns)
        shiftinds=crbs.revruns{shiftnum};
        ind_early=find(mod(crbs.rawtimes(shiftinds,1),1)<(11/24));
        ind_late=find(mod(crbs.rawtimes(shiftinds,2),1)>(15/24));
        if(~isempty(ind_early))
            phbsoutrev(bsvl).early=[phbsoutrev(bsvl).early makerow(shiftinds(ind_early))];
        end
        if(~isempty(ind_late))
            phbsoutrev(bsvl).late=[phbsoutrev(bsvl).late makerow(shiftinds(ind_late))];
        end
    end
    
end