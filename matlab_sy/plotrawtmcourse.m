%plot all upshifts 
%different subplot, include bas runs.

%
ct=1;
totct=1;
inputbs='phar'
if(inputbs=='stim')
    inbs=sumbs
else
    inbs=phsumbs
end

for ii=1:length(inbs)
    crbs=inbs(ii);
    for jj=1:length(crbs.shiftruns)
     
        if(inputbs=='stim')
        
                    [out,sortind]=sort(crbs.tmvec(crbs.STANRUNS));
                    sortruns=crbs.STANRUNS(sortind);
                    basind=find(ismember(sortruns,crbs.basruns));
                    shind=find(ismember(sortruns,crbs.shiftruns{jj}));
                    outind=find(basind==(shind(1)-1));
                    if(~isempty(outind))
                        basruns=sortruns(basind(outind));
                        inbs(ii).allruns{jj}=[basruns sortruns(shind)];
                    else
                        inbs(ii).allruns{jj}=[sortruns(shind)];
                    end
        else
        
        
        
            cr_shiftruns=crbs.shiftruns{jj};
            basruns=crbs.basruns;
            nonzeroind=find(crbs.mulist~=0);
            mulistnonzero=crbs.mulist(nonzeroind);
            shind=find(ismember(mulistnonzero,cr_shiftruns'));
%             shindloc=find(nonzeroind==shind(1));
            basind=find(ismember(mulistnonzero(shind(1)-1),basruns))
                if(~isempty(basind))
               
                    inbs(ii).allruns{jj}=nonzeroind([shind(1)-1 shind]);
                      
                else
                       inbs(ii).allruns{jj}=nonzeroind(shind);
                end
        end
        end
end

for ii=1:length(inbs)
    for jj=1:length(inbs(ii).allruns)
        ax(totct)=subplot(7,2,ct);
        ps.flip=1;
        ps.adjx=1;
        ps.axin=ax(totct);
        ps.plotz=1;
        ps.plotind=inbs(ii).allruns{jj}
        ps.printvl=ii;
        if (inputbs=='stim')
            plot_tmcourse3stim(inbs(ii),ps)
        else
            plot_tmcourse3(inbs(ii),ps)
        end
            
            hold on;
        
        totct=totct+1;
        ct=ct+1;
    end
    
    if(ct==12)
        figure
        ct=1;
       
    end

end




ps.flip=1