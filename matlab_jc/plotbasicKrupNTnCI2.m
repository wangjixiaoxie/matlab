function plotbasicKrupNTnCI2(Post2NT,ind,finalpointPC,finalpoint,post1,semPCfinal,semfinal,sempost,semNT,pcCIbase,CIbaseCTL,CIbase,postPC1)
            plot([0.1 0.1],[pcCIbase],'b')
            sempostPC1=std(postPC1)/sqrt(length(postPC1));
            plot([3.0 3.0],[mean(finalpointPC)-semPCfinal mean(finalpointPC)+semPCfinal],'b')
            plot([3.5 3.5],[mean(postPC1)-sempostPC1 mean(postPC1)+sempostPC1],'b')            
            %plot([3.4 3.4],[pcCIpost1],'b'); % plot([4.9 4.9],[pcCIpost2],'b')
            plot([0 3.0 3.5],[0 mean(finalpointPC) mean(postPC1)],'b')
            plot([0 0],[CIbase],'r'); % plot([1 1],[CIfirst],'r');%plot([2 2],[CIsecond],'r');%plot([3 3],[CIthird],'r')
            plot([3 3],[mean(finalpoint(ind))-semfinal mean(finalpoint(ind))+semfinal],'r')
            plot([3.5 3.5],[mean(post1(ind))-sempost mean(post1(ind))+sempost],'r');%plot([5 5],[CIpost2],'r')
            plot([0 3 3.5],[0 mean(finalpoint(ind)) mean(post1(ind)) ],'r')
            plot([-1 6],[0 0],'k')
            
            plot([0.05 0.05],[CIbaseCTL],'k')
            plot([3.5 3.5],[mean(Post2NT)-semNT mean(Post2NT)+semNT],'k')
            plot([0 3.5],[0 mean(Post2NT)],'k')

