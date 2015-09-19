function []=plothorizline(pt1, pt2,ls,colval)
if(~exist('ls'))
    ls='-'
end
if(~exist('colval'))
    colval='k';
end
plot([pt1(1) pt2(1)], [pt1(2) pt2(2)],colval,'LineStyle',ls);
