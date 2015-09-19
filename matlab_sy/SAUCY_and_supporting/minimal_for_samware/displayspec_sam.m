%
% cribbed whole from spectrogram.m, not altered
%
function displayspectrogram(t,f,Pxx,isFsnormalized,faxisloc)

if 0
id=f<10000;
f=f(id);Pxx=Pxx(id,:);

id=Pxx<10;
Pxx(id)=10;
end
imagesc(t,f,log(Pxx));set(gca,'YD','n');
if 0    %old version
    %
    % cribbed whole from spectrogram.m, not altered
    %
    % Cell array of the standard frequency units strings
    frequnitstrs = getfrequnitstrs;
    if isFsnormalized,
        idx = 1;
    else
        idx = 2;
    end

    newplot;
    if strcmpi(faxisloc,'yaxis'),
        if length(t)==1
            % surf requires a matrix for the third input.
            args = {[0 t],f,10*log10(abs([Pxx Pxx])+eps)};
        else
            args = {t,f,10*log10(abs(Pxx)+eps)};
        end

        % Axis labels
        xlbl = 'Time';
        ylbl = frequnitstrs{idx};
    else
        if length(t)==1
            args = {f,[0 t],10*log10(abs([Pxx' Pxx'])+eps)};
        else
            args = {f,t,10*log10(abs(Pxx')+eps)};
        end
        xlbl = frequnitstrs{idx};
        ylbl = 'Time';
    end
    hndl = surf(args{:},'EdgeColor','none');

    axis xy; axis tight;
    colormap(jet);

    % AZ = 0, EL = 90 is directly overhead and the default 2-D view.
    view(0,90);

    ylabel(ylbl);
    xlabel(xlbl);
end