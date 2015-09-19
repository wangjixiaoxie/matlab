function [vals,spec] = evtafsim(rsong,fs,templates,OLDWAY,USEERROR);
%[vals,spec] = evtafsim(rsong,fs,templates,OLDWAY,USEERROR);
% returns the threshold in case it was set in the program
% simtaf with the clocks and multiple templates put in
% templates is a matrix with each column is one template
%
if (~exist('USEERROR'))
    USEERROR=0;
end
if (~exist('OLDWAY'))
    OLDWAY=0;
end
blen = size(templates,1);
nfft = blen*2;
hammy  = hamming(2*blen);

ntempl = size(templates,2);
nrep = floor(length(rsong)/nfft)-1;

vals = zeros([nrep,ntempl]);
spec = zeros([nrep,blen]);

% templates are normalized just in case
for jj = 1:ntempl
    if (OLDWAY==0) % new way = subtract min then divide by max
        templates(:,jj) = templates(:,jj)-min(templates(:,jj));
        templates(:,jj) = templates(:,jj)./max(templates(:,jj));
    else
        normtmp = templates(:,jj).'*templates(:,jj);
        templates(:,jj) = templates(:,jj)./sqrt(normtmp);
    end
end

for ii = 1:nrep % for each fft window, NON OVERLAPPING (256 each)
	ind1 = (ii-1)*nfft+ 1;
	ind2 = ind1 + nfft - 1;
	datchunk = rsong(ind1:ind2) - mean(rsong(ind1:ind2));
	fdatchunk = abs(fft(hammy.*datchunk)); % gets 256 bins (symmetrical)

	sp = abs(fdatchunk(1:blen)); % takes 128 bins.
    if (USEERROR==0) % default is 0
        sp(1:6)=0.0;
    else
        %sp(2:end)=sp(1:end-1);
        %sp(1:6)=0.0;
        sp = [zeros([6,1]);sp(6:end-1)];
    end
    if (OLDWAY==1)
        normtmp = sqrt(sp.'*sp);
        sp = sp./normtmp;
    else
        sp = sp-min(sp);
        sp = sp./max(sp);
    end
    for jj = 1:ntempl
        if (OLDWAY==1)
            vals(ii,jj) = acos(sp.'*templates(:,jj));
        else
            vals(ii,jj) = (sp-templates(:,jj)).'*(sp-templates(:,jj));
        end
    end
    spec(ii,:) = sp.';
end
return;
