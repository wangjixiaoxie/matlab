%function get_dirtimes called by pitchsyntaxanal2 to pick out the times
%from the rawtimes which are dirtimes.

function [dirtimes] =get_dirtimes(avls)
dirtimes=[];    
if(~isempty(avls.diron))
        for ii=1:length(avls.diron)
            indvl=avls.diron(ii)
            dirtimes=[avls.rawtimes(indvl,:);avls.diron]
            
        end
end