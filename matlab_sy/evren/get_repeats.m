function [vals2,vals]=get_repeats(labels);
% vals=get_repeats(labels);
% given labels returns the repeat counts for the notes
%

df = diff(fix(labels));
pp = find(df~=0);
if (length(pp)==0)
    vals=[1,length(labels),fix(labels(1))];
else
    vals=zeros([length(pp)+1,3]);
    vals(1,:)=[1,pp(1),fix(labels(1))];
    if (length(pp)>1)
        for ii = 2:length(vals)-1
            vals(ii,:)=[vals(ii-1,2)+1,pp(ii),fix(labels(pp(ii)))];
        end
    end
    vals(end,:)=[vals(end-1,2)+1,length(labels),fix(labels(end))];
end

vals2=zeros([1,length(labels)]);
for ii = 1:size(vals,1)
    ind1=vals(ii,1);
    ind2=vals(ii,2);
    vals2(ind1:ind2)=[ind1:ind2]-ind1+1;
end

return;
