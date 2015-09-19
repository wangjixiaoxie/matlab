function [ax]=plotbox(vls,col, ax)
    
   
    xbnds=vls(1:2)
    ybnds=vls(3:4);
%     axes(ax);
    %plot 4 lines;
    plot([vls(1) vls(2)],[vls(3) vls(3)],'Color',col)
    plot([vls(1) vls(2)],[vls(4) vls(4)],'Color',col)
    plot([vls(1) vls(1)],[vls(3) vls(4)],'Color',col)
    plot([vls(2) vls(2)],[vls(3) vls(4)],'Color',col)