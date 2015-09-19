function LMANacorr=jc_reconstruct(pitchLMAN,pitchNONE)
% Input the part of the pitch file you care about e.g. (pitch(450:850,:))


num_notes=min([size(pitchLMAN,2) size(pitchNONE,2)]);
[ac,psdLMAN]=jc_psd(pitchLMAN(:,1:num_notes));
[ac,psdNONE]=jc_psd(pitchNONE(:,1:num_notes));

% see help deconvwnr for code
LMANacorr=fftshift(real(ifft(psdLMAN-psdNONE)));