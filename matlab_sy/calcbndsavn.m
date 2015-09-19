function [colbnds]=calcbndsavn(avna,ps)
       mnfrow=floor(ps.freqlimbound(1)/ps.freqspacing);
       mxrow=ceil(ps.freqlimbound(2)/ps.freqspacing);
       for ii=1:length(avna)
        maxvl(ii)=max(max(avna{ii}(mnfrow:mxrow,:)))
        minvl(ii)=min(min(avna{ii}(mnfrow:mxrow,:)))
       end
       colbnds=[min(minvl) max(maxvl)];
