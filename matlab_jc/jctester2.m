function A=jctester2(B,lowpitch,highpitch) 
for i=1:length(B)
    A(i).pitchcurves=jc_pitchmat1024(B(i).rawdata,1024,1020,2,lowpitch,highpitch,[1],'obs0',1);
end
