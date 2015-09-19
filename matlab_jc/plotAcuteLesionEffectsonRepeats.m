function [mnpre,sdpre,mnpost,sdpost,indlesion,indours,indMU,indAPV,notetype]=plotAcuteLesionEffectsonRepeats(Experiment);
        countlesion=0;countMU=0;countAPV=0;countours=0;
        notetype=[];indours=[];indlesion=[];indMU=[];indAPV=[];mnpre=[];mnpost=[];sdpre=[];sdpost=[];
        for i=1:length(Experiment)
            mnpre(i)=mean(Experiment(i).distnt{1});
            sdpre(i)=std(Experiment(i).distnt{1});
            mnpost(i)=mean(Experiment(i).distnt{2});    
            sdpost(i)=std(Experiment(i).distnt{2});
            if isequal(Experiment(i).manipulation,'lesion') & isequal(Experiment(i).author,'tim')
                countours=countours+1;
                indours(countours)=i;
            end
            if isequal(Experiment(i).manipulation,'lesion')
                countlesion=countlesion+1;
                indlesion(countlesion)=i;
            else
                if isequal(Experiment(i).manipulation,'lidocaine') | isequal(Experiment(i).manipulation,'muscimol')
                    countMU=countMU+1;
                    indMU(countMU)=i;
                else
                    if isequal(Experiment(i).manipulation,'apv')
                        countAPV=countAPV+1;
                        indAPV(countAPV)=i;
                    end
                end
            end
            if isequal(Experiment(i).type,'sweepy') | isequal(Experiment(i).type,'stacky')
                    notetype(i)=1;
                else if isequal(Experiment(i).type,'sweep')
                    notetype(i)=2;
                else if isequal(Experiment(i).type,'stack') | isequal(Experiment(i).type,'high stack')
                    notetype(i)=3;
                else if isequal(Experiment(i).type,'low stack')
                    notetype(i)=4;
            end;end;end;end
        end
        figure;hold on;
        msize=15;
        plot(mnpre(indlesion),mnpost(indlesion),'k.','Markersize',msize)
        plot(mnpre(indours),mnpost(indours),'g.','Markersize',msize)
        plot(mnpre(indMU),mnpost(indMU),'b.','Markersize',msize)
        plot(mnpre(indAPV),mnpost(indAPV),'r.','Markersize',msize)
        plot([0 20],[0 20])
        xlim([0 14]);ylim([0 20])