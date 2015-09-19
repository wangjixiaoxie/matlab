function pplotter1(shiftedCOMPARISON,mini,maxi,window,step)

figure; hold on;
for i=11:20
    %[pitch_data]=jc_PitchDataVCO(shiftedCOMPARISON(i,:),1024,1020,0.4,mini,maxi);
    st=jc_evstftA(shiftedCOMPARISON(i,:),mini,maxi,window,step);
    %ll=resample(pitch_data,1,4);
    plot(st+180*i)
end