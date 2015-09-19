function [xax,waves]=WaveletAnalysisJC2(norm_resids_chopped)
norm=norm_resids_chopped;
norm=norm';
a=size(norm);
notes=a(1);
ppnote=a(2); %Number of points sampled per note
stepsize=20;
L=ppnote;
Window=2^(nextpow2(L)-1);
steps=floor((L-Window)/20)+1;
numWavelets=nextpow2(L)-2;
qmf=MakeONFilter('Coiflet',1); % The wavelet - I could try different ones.
wc=zeros(1,256);
for i=1:notes
    waveweight=zeros(1,numWavelets);
    for j=1:steps
        beginning=1+stepsize*(j-1);
        ending=Window+stepsize*(j-1);
        resid_chunk=norm(i,beginning:ending);
        wc=wc+abs(FWT_PO(resid_chunk,1,qmf));
        index=PlotWaveCoeff3(wc,1,0);
    end
end
g=257-find(index<0.75);
wc(g)=0;
for k=1:7
    waves(k)=sum(wc(dyad(k)));
end
% %We will ultimately care about the total weight of 1 vs. 2 vs. 3 vs. ...
% totals=zeros(1,numWavelets);
% for i=1:length(note)
%     for j=1:numWavelets
%         totals(j)=totals(j)+note(i).waves(j);
%     end
% end


%This may be off by a factor of two now.
for i=1:numWavelets
     xax(i)=(2^(numWavelets+2-i))/(32000/4000);
end
