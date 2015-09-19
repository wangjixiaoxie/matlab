function h=novelty(labels15a,note)
k=0;
j=0;
for i=1:length(labels15a)
    if labels15a(i)=='/'
        k=k+1;
    else 
        if labels15a(i)=='a'
            j=j+1;
        end
    end
    if j==note
        h=k;
    end
end
            