function Indest=smoother(FFOUT)
m=size(FFOUT);
numberOfNotes=m(1);
lengthOfNote=m(2);
for j=1:numberOfNotes
    for i=2:lengthOfNote-2
        Indest(j,i)=pinterp([i-1;i;i+1], [FFOUT(j,i-1);FFOUT(j,i);FFOUT(j,i+1)]);
    end
end

