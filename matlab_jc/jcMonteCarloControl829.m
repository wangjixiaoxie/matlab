function [zScoreMatrix]=jcMonteCarloControl829(setsPerNote,NUMnotes,norm,startMatrix,window_width)

%sepnotes=setsPerNote*NUMnotes;
%long=floor(length(PITCHES)/sepnotes);

for ii=3 %1:setsPerNote-1
    for jj=2 %1:NUMnotes-1
        if size(startMatrix,2)==2
            firstpoint=startMatrix(jj,2)-window_width; %end of first
            nextpoint=startMatrix(jj+1,1); %beginning of second
        else %center of note
            firstpoint=startMatrix(1,jj);
            nextpoint=firstpoint;
        end
        a=1+(jj-1)*setsPerNote;
        b=1+jj*setsPerNote;
        c=a+ii;
        d=b+ii;
        zScoreMatrix(ii,jj)=jcMonteCarloResid(norm(a).residuals,norm(b).residuals,norm(c).residuals,norm(d).residuals,firstpoint,nextpoint,window_width);
    end
end


    
    