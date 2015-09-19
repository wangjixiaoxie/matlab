%plotscatterxy3
function []=plotscatterxy3(sumplot);
figure

colvec={'b' 'k' 'c'}


for ii=1:3
   
    plot(sumplot(ii).ac, sumplot(ii).mu, 'o','Color',colvec{ii});
          hold on;
            
plot([-10 10],[-10 10],'k-')
end








    
