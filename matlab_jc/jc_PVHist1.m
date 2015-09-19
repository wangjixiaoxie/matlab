function [HistoMatrix,normpd]=jc_PVHist(pitch_data,onset,offset)
%Run this after running jc_PitchData to get the pitch_data array and plotting a couple pitch curves to
%determine good values for onset and offset.  The onset and offsets should
%be the values on the x-axis (these are in units of sample number, not ms).


%Make a data structure chunk with just the points between the onset and the
%offset time.
for note_count=1:length(pitch_data)
    for point_count=onset:offset
        n=point_count-onset+1;
        chunk(note_count).pitches(n)=pitch_data(note_count).pitches(point_count);
    end
end
chunk_size=offset-onset+1;

%Calculate the mean pitch curve by summing over the sum of the columns
for pointer=1:chunk_size
    summater=0;
    for note_count=1:length(chunk)
        summater=summater+chunk(note_count).pitches(pointer);
    end
    mean_pc(pointer)=summater/length(chunk);
end
overall_mean=mean(mean_pc); %mean CM
norm_pc=mean_pc-overall_mean; %f=CM-meanCM --- this takes the average curve and centers it at zero


%Normalize the data
for n=1:length(chunk)
    
    %Center the note at zero by subtracting out its mean
    row_mean=mean(chunk(n).pitches);
    normpd(n).pitches=chunk(n).pitches-row_mean;
    
    %Subtract out the average cuve centered at zero
    normpd(n).pitches=normpd(n).pitches-norm_pc;
    
end

%Combine all the stuff in the area into a matrix so that you can do a
%histogram of the variation.
HistoMatrix=normpd(1).pitches;
for k=2:length(chunk)
    HistoMatrix=[HistoMatrix normpd(k).pitches];
end