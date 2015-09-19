function [wspect]=win_spect(spectin,f_step,f_low,f_high,f_band, order)

%takes input spectrogram and windows freqeuncies using a window that is:
%  0 at F_Low and F_high and rises  to one over F-band Hz.
%  the shape of the roll-off is half of hanning window raised to order power
%  higher order giving steeper roll-off
%  f_step is the spacing of frequency bins of the input spectrogram 
% input spect is assumed to start at 0 hz for lowest bin

%originally written with different conventions; here, convert input parameters
% for conventions of algorithm
F_low=f_low+f_band;
F_high=f_high-f_band;
F_band=f_band;

taper_pts=ceil(F_band/f_step);
dtaper=hanning(2*taper_pts);
taper=zeros(1,size(spectin,1));
ltaper=length(taper);
taper(ltaper-taper_pts+1:ltaper)=dtaper(1:taper_pts);
taper_win=ones(1,size(spectin,1));
low_point=round(F_low/f_step);
high_point=round(F_high/f_step);
taper_win(1:low_point)=taper(ltaper-low_point+1:ltaper);
taper=fliplr(taper);
taper_win(high_point:length(taper_win))=taper(1:length(taper_win)-high_point+1);
taper_win=taper_win.^order;


taper_mat=taper_win'*ones(1,size(spectin,2));
wspect=spectin.*taper_mat;
