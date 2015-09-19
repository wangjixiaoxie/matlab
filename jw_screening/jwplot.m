function jwplot(var,scale)
%function jwplot(var,scale)
%var is the variable to plot with surf. scale is a string, 'log' or 'lin' for color/z axis

%load matfile;

figure
newplot;

if scale=='log'
        args = {10*log10(abs(var)+eps)};
end
if scale=='lin'
	args = {((var)+eps)};
end
    % Axis labels
    xlbl = '';
    ylbl = '';

hndl = surf(args{:},'EdgeColor','none');

axis xy; axis tight;
colormap(jet);

caxis([3.5203,85.9716]);
if scale=='lin'
	caxis([0,1]);
end
% AZ = 0, EL = 90 is directly overhead and the default 2-D view.
view(0,90);

ylabel(ylbl);
xlabel(xlbl);
return;
