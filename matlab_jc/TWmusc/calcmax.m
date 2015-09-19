%called by inactivanal3. 
%designed to calculate the maximum pitch average pitch shift across a
%series of acsf runs.
function avls=calcmax(avls)

for ntind=1:length(avls.muanal)
    strvl=['avls.muanal{' num2str(ntind) '}']
%     if(exist(strvl))
        for shiftnum=1:length(avls.muanal{ntind})
            muind=avls.muanal{ntind}{shiftnum}
            acvls=avls.acmean(ntind,muind);
            indnan=isnan(acvls)
            if(indnan(1)<1)
                diffvl=acvls-avls.initmean{ntind};
                [maxvl, maxind]=max(abs(diffvl));
                
                acmax=acvls(maxind);
                avls.acmax(ntind,shiftnum)=acmax;
                avls.acmaxz(ntind,shiftnum)=(acmax-avls.initmean{ntind})/avls.initsd{ntind};
            end
        end
%     end
end

      