function [shiftoffcomb,revoffcomb,bsnumcomb,shiftnumcomb]=plotstim_revpairs(shiftplot)
shiftoffcomb=[];
revoffcomb=[];
bsnumcomb=[];
shiftnumcomb=[];
for bsnum=1:length(shiftplot)
    crsplot=shiftplot(bsnum);
    [outsize]=size(crsplot.aczpr_vls);
    totrows=outsize(1);
    totcol=outsize(2);
    for rnum=1:totrows
        for cnum=1:totcol
            if(~isempty(crsplot.aczpr_vls{rnum,cnum}))
                if(crsplot.drxn{rnum}=='do')
                    mult=-1;
                else
                    mult=1;
                end
         
                acd=mean(crsplot.aczpr_vls{rnum,cnum});
                mud=mean(crsplot.muzpr_vls{rnum,cnum});
                
                acrevd=mean(crsplot.revacpr_vls{rnum,cnum})
                murevd=mean(crsplot.revmupr_vls{rnum,cnum})
                shiftoff=-mult*mean(acd-mud);
                revoff=-mult*(acrevd-murevd);
                shiftoffcomb=[shiftoffcomb shiftoff];
                revoffcomb=[revoffcomb revoff];
                bsnumcomb=[bsnumcomb bsnum]
                shiftnumcomb=[shiftnumcomb rnum]
            end
        end
    end
    
            
end
        
    
    
