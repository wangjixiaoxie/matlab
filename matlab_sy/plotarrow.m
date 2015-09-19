%plotarrow.m
function [ax]=plotarrow(as)
for ii=1:length(as)
    arrow(as(ii).plotpts(1,:), as(ii).plotpts(2,:), 'EdgeColor',as(ii).col,'FaceColor',as(ii).col)
end