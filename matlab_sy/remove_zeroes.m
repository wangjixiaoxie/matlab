function outarray=remove_zeroes(lblarray)

indkeep=[];
%first test whether all the labels are zeroes
for ii=1:length(lblarray)
    lbls=lblarray{ii}
    ind=find(lbls~='0')
    
    lblarray{ii}=lbls(ind);
    
    if (isempty(ind)==0)
        indkeep=[indkeep ii];
    end
    
end
for ii=1:length(indkeep)

	outarray{ii}=lblarray{indkeep(ii)}
end
    
        
    
