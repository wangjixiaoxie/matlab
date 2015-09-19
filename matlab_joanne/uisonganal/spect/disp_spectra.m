function disp_spectra(spect,n,print_flag);

%displays every nth spectrum from spectrogram spect
%if pflag=1, prints spectra

figure
for i=1:n:size(spect,2)
  plot(spect(:,i))
  pause
  if print_flag==1
    print
  end
end
