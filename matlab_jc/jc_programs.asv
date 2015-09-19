%Reads in a batch file of .not.mat's and stores the data for the pitch
%curves.  On my slow computer, this takes a couple minutes per 100 notes.
function [Lesioned_pitch_data]=jc_PitchData(1024,1010,2,2000,4000,'batchbk93_507.labeled','d','d','d');

%Takes the pitch curves for a bunch of notes and makes a 1xn vector of all
%the variations after normalizing. Normalization: First, it centers the curve 
%at zero by subtracting the mean of the individual curve and second, it 
%subtracts the zero-centered mean curve from the zero-centered individual curve.
function [HistoMatrix,normpd]=jc_PVHist(pitch_data,onset,offset)

%Eliminates outliers - it's a good idea to run this before computing the
%standard deviation.
function [new_HistoMatrix]=jc_simple430(input_matrix,lower_limit,higher_limit)

%Compare the pitch curves for a bunch of notes from each data set
function jc_plotabunch501(arrayfileD,arrayfileUD,arrayfileL)   

%Make overlapping plots of the histogram of variations(red,blue,black)
function jc_histplot501(Dvector,UDvector,Lvector)   

%Compares the autocovariance plots for a bunch of notes from each data set
function jc_autocovplot501(arrayD,arrayUD,arrayL,first_note,last_note,start_bin,stop_bin)
    calls
function jc_autocov501(arrayD,arrayUD,arrayL,first_note,last_note,start_bin,stop_bin)

