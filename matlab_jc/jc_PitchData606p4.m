function med=jc_PitchData606p4(N,OVERLAP,sigma,F_low,F_high,data)
%F_low and F_high are the boundary frequencies in which to look for a peak
%in the autocorrelation function (i.e. in which to calculate the pitch at
%each point in time).
%Recommend jc_PitchData(data,1024,1010,2,2000,4000,'batch.labeled','d','d','d')


%This extends jc_PitchData by parabolically interpolating over adjacent
%frequency bins.

%Sets parameters and runs findwnote4 to extract all instances of the
%desired note in the specified context.
Timeshift=0;
Note_length=6500;
UseFit=1;
ChanSpec='obs0';
%[fvalsstr]=findwnote4(batchFWN,noteFWN,prenoteFWN,postnoteFWN,Timeshift,[F_low F_high],Note_length,UseFit,ChanSpec);

%Parameters that one may want to change.
SAMPLING=32000;       
ZOOMT=1;
ZOOMF=3;
TL=5;
FL=5;



    
    %Calls Tim Gardner's program, which takes the phase into
    %consideration and perhaps even uses consensus of different sampling rates
    %to determine the instantaneous fourier transform.
    %current_note=fvalsstr(found_note).datt;
    %chopped_note=current_note(1000:5500);
    [med]=ifdvCONplot1(data,SAMPLING,N,OVERLAP,sigma,ZOOMT,ZOOMF,TL,FL); 
    


        
        
        
        
        
        