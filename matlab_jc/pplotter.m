function pplotter(shiftedCOMPARISON,mini,maxi,window,shift)

figure; hold on;
diffs=shift/2;
shifter=(512-window/2)/shift;
for i=11:20
    [pitch_data]=jc_PitchDataVCO(shiftedCOMPARISON(i,:),1024,1020,0.4,mini,maxi);
    st=jc_evstftA(shiftedCOMPARISON(i,:),mini,maxi,window,shift);
    ll=resample(pitch_data,1,diffs);
    vv=resample(st,1,2);
    plot(ll+i*180,'r');plot(vv(1+shifter/2:length(vv))+i*180)
end