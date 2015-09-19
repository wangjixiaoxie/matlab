            thisbird=1;
            wins=[-820:10:1420];
            numtbins=length(wins);
            figure;hold on;
            for sylls=1:3
                subplot(4,3,sylls);hold on;
                plot(wins,Bird(thisbird).varsUD.Tmn(sylls,:),'Linewidth',4)
                plot(wins,Bird(thisbird).varsUD.Tmn(sylls,:)+Bird(thisbird).varsUD.Tse(sylls,:),'Linewidth',2)
                plot(wins,Bird(thisbird).varsUD.Tmn(sylls,:)-Bird(thisbird).varsUD.Tse(sylls,:),'Linewidth',2)
                ylim([-2 3])
                plot([-820 1420],[0 0],'k','Linewidth',2)
                xlim([-150 150])
            end
            for sylls=1:3
                subplot(4,3,sylls+3);hold on;
                plot(wins,Bird(thisbird).varsD.Tmn(sylls,:),'Linewidth',4)
                plot(wins,Bird(thisbird).varsD.Tmn(sylls,:)+Bird(thisbird).varsD.Tse(sylls,:),'Linewidth',2)
                plot(wins,Bird(thisbird).varsD.Tmn(sylls,:)-Bird(thisbird).varsD.Tse(sylls,:),'Linewidth',2)
                ylim([-2 3])
                plot([-820 1420],[0 0],'k','Linewidth',2)
                xlim([-150 150])
            end
            thisbird=1;
            for sylls=1:3
                subplot(4,3,sylls+6);hold on;
                plot(wins,Bird(thisbird).varsUD.Rmn(sylls,:),'Linewidth',4)
                plot(wins,Bird(thisbird).varsUD.Rmn(sylls,:)+Bird(thisbird).varsUD.Rse(sylls,:),'Linewidth',2)
                plot(wins,Bird(thisbird).varsUD.Rmn(sylls,:)-Bird(thisbird).varsUD.Rse(sylls,:),'Linewidth',2)
                ylim([-0.5 0.5])
                plot([-820 1420],[0 0],'k','Linewidth',2)
                xlim([-150 150])
            end
            for sylls=1:3
                subplot(4,3,sylls+9);hold on;
                plot(wins,Bird(thisbird).varsD.Rmn(sylls,:),'Linewidth',4)
                plot(wins,Bird(thisbird).varsD.Rmn(sylls,:)+Bird(thisbird).varsD.Rse(sylls,:),'Linewidth',2)
                plot(wins,Bird(thisbird).varsD.Rmn(sylls,:)-Bird(thisbird).varsD.Rse(sylls,:),'Linewidth',2)
                ylim([-0.5 0.5])
                plot([-820 1420],[0 0],'k','Linewidth',2)
                xlim([-150 150])
            end