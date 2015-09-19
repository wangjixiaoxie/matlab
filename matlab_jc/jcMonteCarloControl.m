function [zScoreMatrix,norm]=jcMonteCarloControl(setsPerNote,NUMnotes,pitchA,pitchB,pitchC,pitchD,pitchE,pitchF,startMatrix,window_width)

%Make numerical indices
if NUMnotes==2;
    PITCHES=[pitchA pitchB];
end
if NUMnotes==5;
    PITCHES=[pitchA pitchB pitchC pitchD pitchE];
end
if NUMnotes==6;
    PITCHES=[pitchA pitchB pitchC pitchD pitchE pitchF];
end
sepnotes=setsPerNote*NUMnotes;
long=floor(length(PITCHES)/sepnotes);

%Calculate residuals
for i=1:sepnotes
    starter=(i-1)*long+1;
    ender=i*long;
    norm(i).residuals=jc_plotresiduals725b(PITCHES(starter:ender));
end

for ii=1:setsPerNote-1
    for jj=1:NUMnotes-1
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


    
    