function [spkwid]=Find_Spk_Width(inarray,thresh);
    sfct=(1e3/2^15);
    TH=-sfct*thresh;
    inarray=inarray-thresh/2;
    crs_ind=find(inarray>0);
    spkwid=0;
    if(isempty(crs_ind))
        spkwid=0;
        return;
    end
    cnscvc=getconsec(crs_ind);
    lnvec=length(cnscvc(:,1));
    
    if(lnvec>1)
        for ii=1:lnvec
            testvec=crs_ind(cnscvc(ii,1):cnscvc(ii,2));
            ind=find(testvec==20);
            if(isempty(ind)==0)
                spkwid=(cnscvc(ii,2)-cnscvc(ii,1)+1);
            end
        end

    else
        
        spkwid=(cnscvc(1,2)-cnscvc(1,1)+1);
    end    