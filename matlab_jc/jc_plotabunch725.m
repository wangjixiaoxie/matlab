function CVs=jc_plotabunch725(pitch_data,f,xmax,avgnote,alpha,beta)

%This program is especially designed to analyze the pitch of notes that
%receive white noise aversive reinforcement.  The problem with aligning
%these notes is that the loud white noise overlaps with the note and thus
%screws up the alignment.  To avoid this, I align using a part of the 
%smoothed rectified amplitude profile that occurs before the white noise 
%onset.
wn_onset=2300; %white noise onset in note - this parameter should be 
               %adjusted so that it is large enough to include enough 
               %information to allow alignment but small enough so that it
               %doesn't include the white noise burst.
               
               
%This program returns the median value of the note between the onset of the
%smooth stack and the onset of the white noise burst.  Choose these
%parameters after plotting the curves.
stack_start=600;
wn_start=750;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spacing=180;
note_length=6500;

sampling_rate=32000;

%Note that xmax should be around 800 or so typically.
xmin=0;
ymin=200;
ymax=10000; 


%Smooth the data
for i=1:length(pitch_data)
    [holder]=SmoothData(f(i).datt,sampling_rate,1);
    smooth(i).smoothed=holder;
end


figure; hold on
for i=1:length(pitch_data)
    h=xcov(avgnote(1:wn_onset),smooth(i).smoothed(1:wn_onset));
    [peak,index]=max(h);
    index=index+(note_length-wn_onset);
    shift(i)=(index-note_length)/4;
    %smooth
    processed=jc_pinterp(pitch_data(i).pitches);
    %shift them over to align
    k=1;
    for j=1:xmax
        align_shift=shift(i);
        shift_index=round(j-align_shift);
        if shift_index>0
            shifted(i,k)=processed(shift_index);
        else
            shifted(i,k)=ymin;
        end
        k=k+1;
    end
    plotshiftedD=shifted(i,:);
    %meds(i)=median(plotshiftedD(stack_start:wn_start));
    CVs(i)=mean(plotshiftedD(alpha:beta));
    plotshiftedD=plotshiftedD+spacing*(i);

    plot(plotshiftedD); xlim([xmin xmax]); ylim([ymin ymax]); title('1-10')
end
