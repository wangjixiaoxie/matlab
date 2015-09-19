function [toff]=get_trgoff(trigs);


dts=[];
tvl=[];
for ii=1:length(trigs)
    
	
    if isempty(trigs(ii).toffset)
        trigs(ii).toffset=-50;
    end
    dn=fn2datenum(trigs(ii).fn);
	datesadd=dn*ones(length(trigs(ii).toffset),1)
    
    dts=[dts;datesadd];
    

    tvl=[tvl;trigs(ii).toffset];
    
end

toff=[tvl dts];