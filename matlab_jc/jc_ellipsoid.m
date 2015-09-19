function jc_ellipsoid(scores)
xc=mean(scores(:,1));
yc=mean(scores(:,2));
zc=mean(scores(:,3));
xr=std(scores(:,1));
yr=std(scores(:,2));
zr=std(scores(:,3));

ellipsoid(xc,yc,zc,xr,yr,zr);
