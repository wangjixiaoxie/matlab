function [minind]=Find_Spk_Width2(waves);
    waves2=waves(:,20:40,:);
    [minval,minind]=min(waves2,[],2);
    