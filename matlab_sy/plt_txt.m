
function h=plt_txt(lblarray)

   figure
   for ii=1:length(lblarray)
       text(1,ii*3, lblarray{ii})
       hold on;
   end