function lblarray=getlabels(bt)

    files=load_batchf(bt);
    
    for ii=1:length(files)
        try
            fn=files(ii).name;
            
            matfn=[fn '.not.mat'];
            strcmd2=['load ' matfn ';']
            eval(strcmd2)
            lblarray{ii}=labels;
        catch err
            continue
        end
    end