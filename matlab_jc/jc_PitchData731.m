function [pitch_data,P]=jc_PitchData731(shifted,N,OVERLAP,sigma,F_low,F_high,window,shift)
%F_low and F_high are the boundary frequencies in which to look for a peak
%in the autocorrelation function (i.e. in which to calculate the pitch at
%each point in time).
%Recommend jc_PitchData(data,1024,1010,2,2000,4000,'batch.labeled','d','d','d')


%This extends jc_PitchData by parabolically interpolating over adjacent
%frequency bins.

m=size(shifted);
NumberOfNotes=m(1);

%Parameters that one may want to change.
SAMPLING=32000;      

for found_note=1:NumberOfNotes
    current_note=shifted(found_note,:);
    [sonogram,F,T,P]=ifdv(current_note,SAMPLING,N,OVERLAP,sigma); 
    %The following is some processing of the output from Tim Gardner's
    %program that facilitates indexing and calculations in the rest of the
    %program.
    a=size(sonogram);
    total_bins=a(2);
    Nyquist=SAMPLING/2;
    step=Nyquist/a(1); %This is the width (in ms) between bins.
    mini=round(F_high/step);
    maxi=round(F_low/step);
    highest_harmonic=3; 
    
    %This loop determines the pitch (fundamental frequency) at each time bin by
    %taking the peak of the autocorrelation function within a range specified
    %by the user.
    for currenttime_bin=1:total_bins
        slice=sonogram(:,currenttime_bin); %power at each frequency for the current time bin
        Powerful(1)=0;
        Freqbinest(1)=0;
        for i=1:highest_harmonic %:highest_harmonic %:highest_harmonic
            minimum=a(1)-mini*i;  
            maximum=a(1)-maxi*i;
            freq_window=slice(minimum:maximum);
            [Pow,Ind]=max(freq_window);  
            % Interpolation--subtraction of 0.5 because we care about the
            % central frequency of the bin, not the upper boundary.
            if Ind==1 || Ind==length(freq_window)
                Indest=Ind;
            else
                Indest=pinterp([Ind-1;Ind;Ind+1], [freq_window(Ind-1);freq_window(Ind);freq_window(Ind+1)]);
            end
            Real_Index=minimum+Indest-1.5; % -1 to account for window;(-1?!) to account for weird error in ifdv
            Freqbinest(i)=(a(1)-Real_Index)/i;
            Powerful(i)=Pow;
        end
        normalizer=sum(Powerful);
        normpower=Powerful/normalizer;
        freqbin_estimate(currenttime_bin)=dot(normpower,Freqbinest);
    end
    pitch_data(found_note).pitches=freqbin_estimate*step;
end

diffs=shift/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PLOTTING%%%%%%%%%%%%%%%%%%%%%%%
%Calculate mean to use for the red hashed line
for i=1:length(pitch_data(1).pitches)
    for j=1:length(pitch_data)
        matr(i,j)=pitch_data(j).pitches(i);
    end
    averaged(i)=mean(matr(i,:));
end
mmean=mean(averaged(500:800));
for i=1:length(pitch_data(1).pitches)/diffs
    t(i)=mmean;
end



%Plotting - like pplotter

shifter=(512-window/2)/shift;


figure; hold on;subplot(141); hold on;xlim([0 150]);ylim([F_low F_high+1800])
for i=1:10
    plotshiftedD=pitch_data(i).pitches+180*(i);
    %plot(plotshiftedD)
    plot(t+180*i,'LineStyle',':','Color',[1 0 0])
    st=jc_evstftA(shifted(i,:),F_low,F_high,window,shift);
    ll=resample(plotshiftedD,1,diffs);
    vv=resample(st,1,2);
    plot(ll,'r');plot(vv(1+shifter/2:length(vv))+i*180)
end
subplot(142); hold on; xlim([0 250])
for i=11:20
    plotshiftedD=pitch_data(i).pitches+180*(i);ylim([F_low+1800 F_high+3600])
    %plot(plotshiftedD)
    plot(t+180*i,'LineStyle',':','Color',[1 0 0])
    st=jc_evstftA(shifted(i,:),F_low,F_high,window,shift);
    ll=resample(plotshiftedD,1,diffs);
    vv=resample(st,1,2);
    plot(ll,'r');plot(vv(1+shifter/2:length(vv))+i*180)
end
        
        
        
        