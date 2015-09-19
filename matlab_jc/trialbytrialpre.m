function [Experiment]=trialbytrialpre(Experiment,expnumber,fvBaseline,fvWNon)

for i=1:length(fvBaseline)
    shiftedBaseline(i,:)=fvBaseline(i).datt;
end
pitchBaseline=jc_pitchmat1024(shiftedBaseline,1024,1020,1,7400,8800,[1],'obs0',1);
        recognized=[];
        hitrange=[];
        catchtrial=[];
        for i=1:length(fvWNon)
            shiftedWNon(i,:)=fvWNon(i).datt;
            recognized(i)=fvWNon(i).CATCH>-1;
            hitrange(i)=fvWNon(i).TEMPLATE==1;
            catchtrial(i)=fvWNon(i).CATCH==1;
        end
        indEscapeAbove=recognized & hitrange & catchtrial;
        indHitAbove=recognized & hitrange & ~catchtrial;
pitchWNon=jc_pitchmat1024(shiftedWNon,1024,1020,1,7400,8800,[1],'obs0',1);    
Experiment(expnumber).pitchBaseline=pitchBaseline;
Experiment(expnumber).pitchWNon=pitchWNon;
Experiment(expnumber).indEscapeAbove=indEscapeAbove;
Experiment(expnumber).indHitAbove=indHitAbove;
