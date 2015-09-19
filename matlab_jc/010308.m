function 010308(trigs)
% hits vs. escapes
allescapes=[]; % 1 means escape
for n=1:length(trigs)

    for j=1:length(trigs(n).ntmintmp)
        g(j)=isempty(find(trigs(n).trigmintmp==trigs(n).ntmintmp(j)));
    end
    for j=1:length(trigs(n).trigmintmp)
        h(j)=isempty(find(trigs(n).ntmintmp==trigs(n).trigmintmp(j)));
    end

    % Get feedback (H/E) info from rec file
    rd=jcreadrecf2(trigs(n).fn,0);

    % For each note, classify as hit or escape --- 1 means escape
    rdcounts=rd.catch(find(h==0));
    count=1;
    for k=1:length(trigs(n).ntmintmp)
        if g(k)==0
            escapes(k)=rdcounts(count);
            count=count+1;
        else
            escapes(k)=1;
        end
    end
    allescapes=[allescapes escapes];
end