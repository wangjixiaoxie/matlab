[fvalsstr,pitch_data]=jc_PitchData611(1024,1022,1,1900,2700,'testbatch1','d','d','d');
%divide the last three by two if you use 1024/1020 instead of 1024/1022 
norm_resids=jc_getresiduals616(pitch_data,fvals,1200,200,1000)
%At this point I should plot a bunch and choose a good interval to chop.
%plotresids(norm_resids)

function totals=WaveletAnalysis(norm_resids_chopped)
norm=norm_resids_chopped;
norm=norm';
a=size(norm);
notes=a(1);
ppnote=a(2); %Number of points sampled per note
stepsize=20;
L=ppnote;
Window=2^(nextpow2(L)-1);
steps=floor((L-Window)/20)+1;

qmf=MakeONFilter('Symmlet',4); % The wavelet - I could try different ones.
for i=1:notes
    waveweight=zeros(1,7);
    for j=1:steps
        beginning=1+stepsize*(j-1);
        ending=Window+stepsize*(j-1);
        resid_chunk=norm(i,beginning:ending);
        wc=FWT_PO(resid_chunk,1,qmf);
        returnee=PlotWaveCoeff2(wc,1,0);
        for k=1:length(waveweight)
            waveweight(i)=waveweight(i)+sum(abs(returnee(k).bumps));
        end
    end
    waveweight=waveweight/length(waveweight);
    %Weigh the big variations by their actual weight in the real residual.
    maximum=max(resid_chunk);
    minimum=min(resid_chunk);
    weighting=maximum-minimum;
    waveweight=waveweight*weighting; 
    note(i).waves=waveweight;
end

%We will ultimately care about the total weight of 1 vs. 2 vs. 3 vs. ...
totals=zeros(1,3);
for i=1:length(note)
    totals(1)=totals(1)+note(i).waves(1);
    totals(2)=totals(2)+note(i).waves(2);
    totals(3)=totals(3)+note(i).waves(3);
end
    
    
