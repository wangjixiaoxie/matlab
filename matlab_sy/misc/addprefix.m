function [pathvl]=addprefix(pathvl, prefix)
for ii=1:length(pathvl)
    if(~isempty(pathvl{ii}))
        pathvl{ii}=[prefix pathvl{ii}]
    end
end
