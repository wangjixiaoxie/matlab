function ax=plot_conting(switchdts, firstval, thresh, dir, freqbnds,scale,scalebnds)
    
    %figure
    %firstval=vals(1,1)
    timevals=switchdts-firstval;
    for ii=1:length(timevals)-1
        y1=thresh(ii);
        if(dir(ii))
            y2=freqbnds(2)
            if(scale)
                y2=scalebnds(2);
            end
        else
                y2=freqbnds(1)
                if(scale)
                    y2=scalebnds(1);
                end
        end
        x1=timevals(ii)
        x2=timevals(ii+1);
        if(scale)
            y1=(thresh(ii)-freqbnds(1))/(freqbnds(2)-freqbnds(1))*(scalebnds(2)-scalebnds(1))+scalebnds(1);
        end  
            ax=fill([x1 x2 x2 x1],[y1 y1 y2 y2],[1 .806 .817]);
        hold on;
    end  
        %return;



