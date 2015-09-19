


for runvl=1:length(avls.pvls)
         strcmd=['cd '  avls.pvls{runvl}];
         eval(strcmd);
        strcmd=['load ' avls.pvls{runvl} avls.cvl{runvl} '.mat'];
            eval(strcmd);
        for ii=1:length(avls.NT)
            crnote=ii;
            vals=getvals(fv{crnote},1,'trig')
            mn=mean(vals(:,2))
            stdout=std(vals(:,2))
            mnlist(runvl,crnote)=mn;
            stdlist(runvl,crnote)=stdout;
        end
end

figure

%plotting
for runvl=1:36
for nt=1:3

    ax(nt)=subplot(3, 1 ,nt)
    timevls=avls.tmvc(:,1)-avls.tmvc(1,1)+1;
    plot([timevls(runvl) timevls(runvl) ],[mnlist(runvl,nt)+stdlist(runvl,nt) mnlist(runvl,nt)-stdlist(runvl,nt)])
    hold on;
end
end

linkaxes(ax,'x')
axes(ax(1))
plot([1 18],[mean(mnlist(1:2,1)) mean(mnlist(1:2,1))],'k--')

plot([7 18],[mean(mnlist(4:8,1)) mean(mnlist(4:8,1))],'c--')

plot([13 18],[mean(mnlist(9:12,1)) mean(mnlist(9:12,1))],'r--')

axes(ax(2))
plot([1 18],[mean(mnlist(1,2)) mean(mnlist(1,2))],'k--')
hold on;
plot([11 13],[mean(mnlist(8,2)) mean(mnlist(8,2))],'c--')
plot([13 18],[mean(mnlist(9:12,2)) mean(mnlist(9:12,2))],'r--')

axes(ax(3))
plot([1 18],[mean(mnlist(1:8,3)) mean(mnlist(1:8,3))],'k--')
plot([13 18],[mean(mnlist(9:12,3)) mean(mnlist(9:12,3))],'r--')

