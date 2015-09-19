function []=plotarrows2(xvec,outacmn,outmumn,ps)
%         axes(axhist);
        arrow_col=ps.col;
        
%         shiftlen=shiftarrow(3)-shiftarrow(1);
%         headwidth(1)=.8/abs(origlen);
% %         headwidth(2)=.01/abs(shiftlen);
%         headht(1)=.01/abs(origlen)
%         headht(2)=.01/abs(shiftlen)
if (isfield(ps,'axbnds'))
    axis(ps.axbnds)
end
if(isfield(ps,'pct'))
    pct=ps.pct
else
    pct=100
end
if(isfield(ps,'arrowind'))
    arrowind=ps.arrowind
else
    arrowind=1:length(xvec);
end

    
axis square
        for ii=arrowind
            crind=ii
            xvl=xvec(crind);
            y1=outacmn(crind);
            y2=outmumn(crind);
            arrowh([xvl xvl],[y1 y2],arrow_col,[pct pct]);
            plot([xvl xvl],[y1 y2],'Color',arrow_col,'Linewidth',2)
        end

