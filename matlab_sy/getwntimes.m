function [wnoutreg, wnrevout] = getwntimes(wnnorm, wnrev)
%10.31.08 called by inactivanal3...
%outputs jesustime for wnon and wnoff
if (exist('wnrev'))
    numinputs=2;
else
    numinputs=1;
end
for it=1:numinputs
    clear wnout
    if it==1
        wn=wnnorm
    else
        wn=wnrev
    end




    for ii=1:length(wn)
        if ~isempty(wn(ii).tmon)
            for jj=1:length(wn(ii).tmon)
                for kk=1:length(wn(ii).tmon{jj})
                    wnout(ii).on(jj,kk)= datenum(wn(ii).tmon{jj}{kk},'yyyy-mm-dd HH')
                    wnout(ii).off(jj,kk)=datenum(wn(ii).tmoff{jj}{kk},'yyyy-mm-dd HH')
                end
            end
        else
            wnout=[];
        end
    end

if it==1
    wnoutreg=wnout;
else
    wnrevout=wnout;
end
end

