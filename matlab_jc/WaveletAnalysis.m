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
    waveweight=zeros(1,8);
    for j=1:steps
        beginning=1+stepsize*(j-1);
        ending=Window+stepsize*(j-1);
        resid_chunk=norm(i,beginning:ending);
        wc=FWT_PO(resid_chunk,1,qmf);
        returnee=PlotWaveCoeff2(wc,1,0);
        for k=1:length(waveweight)
            waveweight(k)=waveweight(k)+sum(abs(returnee(k).bumps));
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
totals=zeros(1,8);
for i=1:length(note)
    for j=1:8
        totals(j)=totals(j)+note(i).waves(j);
    end
end