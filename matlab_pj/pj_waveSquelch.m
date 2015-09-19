function pj_waveSquelch(cbinFile,datChan,waveFile,fs,squelchThresh)
%% Take as input a cbin file and channel, zero all entries below the threshold 
% squelchThresh, and write a wave file named waveFile. 

d = ReadCbinFile(cbinFile);
chanDat = d(:,datChan);

chanDat(abs(chanDat) < squelchThresh) = 0;

wavwrite(chanDat,fs,waveFile);

end

