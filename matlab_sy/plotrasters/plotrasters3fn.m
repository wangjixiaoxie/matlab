%This function takes an Rx2 cell array where R is number of rasters to be plotted.
%C1 is 2 column matrices with spike times and trials, and
%a cell array with titles of the matrices. C2 is the string with name of the stimulus It scales the size of the
%plotted output to the  largest number of trials SO THAT each spike is the
%same size....

%Rewritten 2.13.06 for diferrent daya format, one cell array, where trial
%is cellarray{trial), first column is number of spikes

function plotrasters3fn(inarray,numrast)

    for rast=1:numrast;
        ax(rast)=subplot(numrast,1,rast);
        plotrasters3(inarray);
        lenlist(rast)=length(inarray{rast}(:,1));
        %title(inarray{rast,2});
    end
    
    linkaxes(ax);
    axis([0 12 0 20])   

