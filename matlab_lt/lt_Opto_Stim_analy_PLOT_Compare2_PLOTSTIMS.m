    % Plot Stim times
    % Figure out what y value to place marks
    Ylimits=ylim;
    Yrange=Ylimits(2)-Ylimits(1);
    Y1=Ylimits(1);
    Y2=Ylimits(1)+Yrange/20;
    
    X=StatsStruct.(fieldname).TimeSinceLastTrig; % matrix of stim times (preceding alignment time)
    X=PreDur-X; % time from start of plot, ms
    
            if StimDur==0; % i.e. plot a dot
                plot(X,Y2,'^r','MarkerFaceColor','r','MarkerSize',5);
            else
                plot(X,Y1,'^k','MarkerFaceColor','k','MarkerSize',5); % start
                plot(X+StimDur,Y1,'^k','MarkerFaceColor','k','MarkerSize',5); % end
                
                % put lines for all stims
                for kk=1:length(X);
                    Yrand=rand*(Yrange/15); % use rand to have lines not entirely overlap
                    line([X(kk) X(kk)+StimDur], [Y2+Yrand Y2+Yrand],'LineWidth',1,'Color','r'); 
                end
            end    
