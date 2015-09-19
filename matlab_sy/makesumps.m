%makesumps

clear ps
for ii=1:length(sh)
    for jj=1:length(sh(ii).subruns)
        ps.exclude{ii}(jj)=0;
    end
end

