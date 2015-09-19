function imavn(avn);
%

imagesc(([0:size(avn,2)-1]-10)*3.2,[0:size(avn,1)-1]*62.5e-3,log(avn));
syn;ylim([0,10]);grid on;
return;
