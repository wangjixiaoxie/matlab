function h=lt_plot_stem3(X, Y, Z, color, plotmesh)
%% plots stem3, and adds mesh at z=0;
% automatically plots +, - values diff shades

indstmp=Z>0;
h=stem3(X(indstmp), Y(indstmp), Z(indstmp));
set(h, 'Color',color, 'MarkerFaceColor',color)

indstmp=Z<0;
h=stem3(X(indstmp), Y(indstmp), Z(indstmp));
set(h, 'Color',color)

if plotmesh==1
x=xlim; x=linspace(x(1), x(2), 10);
y=ylim; y=linspace(y(1), y(2), 10);
z=zeros(length(x), length(y));
hmesh=mesh(x, y, z);

set(hmesh, 'FaceAlpha', 0)
set(hmesh, 'EdgeColor', [0.5 0.5 0.5])
set(hmesh, 'LineWidth', 0.3)

end




