  function []=plotarrows(origarrow,shiftarrow,ps)
%         axes(axhist);
        
        origlen=origarrow(3)-origarrow(1);
        shiftlen=shiftarrow(3)-shiftarrow(1);
        headwidth(1)=.01/abs(origlen);
        headwidth(2)=.01/abs(shiftlen);
        headht(1)=.01/abs(origlen)
        headht(2)=.01/abs(shiftlen)
        if(ps.plotshiftarrow)
            h=arrow([origarrow(1) origarrow(2)],[origarrow(3) origarrow(4)],'Length',5);
            set(h,'FaceColor',ps.shift_col)
            set(h,'EdgeColor','none')
        end
        if(ps.plotrevarrow)
            h=arrow([shiftarrow(1) shiftarrow(2)],[shiftarrow(3) shiftarrow(4)],'Length',5)
            set(h,'FaceColor',ps.rev_col)
            set(h,'EdgeColor','none')
        end
        