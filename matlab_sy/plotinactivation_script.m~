%plotinactivation_script

%the purpose of this is to make figures for a basic inactivation/pitchshift
%figure


 %1 
 %3 
 %4 is barplot of zscore offset.
 %5 is barplot of pct offset.
 
 figstoplot=[4 5]
 
 if(find(figstoplot==4)||find(figstoplot==5))
     modstruct.bsnum=7;
     modstruct.shftnum=1;
     modstruct.runs=20
     [sumplotout,sumplotrevout]=selectinbasruns3(sumbs,1, {1:5 1:5}, [1.5 0],modstruct)
 end
 
if (find(figstoplot==4))
     figure
     [vlsbar]=plotbarinitind5(sumplotout, 'bar',1,0)

end
     
   
     
     
if(find(figstoplot==5))
    [vlspct]=plotbarinitind5(sumplotout,'pct',1,0)
end