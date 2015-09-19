for crind=1:length(avls.pvls)
    if(~isempty(avls.pvls{crind}))
    pathvl=avls.pvls{crind}
    
    if(isfield(avls,'baspath'))
        cmd=['cd ' avls.baspath pathvl]
    else
        cmd=['cd ' pathvl]
    end
    eval(cmd);
    btfile=[avls.cvl{crind}];
    cmd=['changelabel(''' btfile ''',''a'', ''0'', ''0'', ''z'');'] 
    eval(cmd);
    end
end