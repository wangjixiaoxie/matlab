function [spkinds] = FindSpikeTimes(rawdata,Fs,USECHANSPEC,TH,PKVAL,SpkSrch,SameSpkWin);
% rawdata is a matrix with all the data
% Fs - sampling freq in Hz
% USECHANSPEC is the collumn which has the neural data we want to plot
% TH is the level at which to threshold for spikes
% PKVAL = 0 look for max
% PKVAL = 1 look for min
% SPKSRCH is a 2 element vector which indicates in seconds 
%         how far backward and forward to look for a min/max 
%         from a threshold corssing event for a spike
% SAMESPKWIN (in seconds) - if two spike events are less than SAMESPKWIN
%                    apart they are considered the same event
%


%change to samples
SpkSrch = [floor(SpkSrch(1)*Fs),ceil(SpkSrch(2)*Fs)];
SameSpkWin = ceil(SameSpkWin*Fs);

spkvals_sv=[];
for iChan = 1:length(USECHANSPEC)
    tmpdat=rawdata(:,USECHANSPEC(iChan));
    if (PKVAL==0)
        pp=find((tmpdat(1:end-1)<=TH)&(tmpdat(2:end)>TH));
    else
        pp=find((tmpdat(1:end-1)>=TH)&(tmpdat(2:end)<TH));
    end
    
    spkvals=zeros([length(pp),2]);
    for ii=1:length(pp)
        strng=pp(ii)+SpkSrch(1);
        enrng=pp(ii)+SpkSrch(2);
        strng=max([1,strng]);
        enrng=min([length(tmpdat),enrng]);


        if (PKVAL==0)
            [tval,tind] = max(tmpdat(strng:enrng));
        else
            [tval,tind] = min(tmpdat(strng:enrng));
        end
        spkvals(ii,:)=[tind+strng-1,tval];
    end
    spkvals_sv=[spkvals_sv;spkvals];
end

% go through all the spikes and any spike times within SameSpkWin
% take the channel with the largest/smallest (PKVAL) amplitude as the 
% time of the spike

[y,i]=sort(spkvals_sv(:,1));
spkvals_sv=spkvals_sv(i,:);
Nspkinds=0;
spkinds=zeros(size(spkvals_sv));
ind1=1;ind2=2;
while (1)
    if ((spkvals_sv(ind2,1)-spkvals_sv(ind1,1))>SameSpkWin)|(ind2>=size(spkvals_sv,1))
        Nspkinds=Nspkinds+1;
        if (PKVAL==0)
            [y,i]=max(spkvals_sv(ind1:(ind2-1),2));
        else
            [y,i]=min(spkvals_sv(ind1:(ind2-1),2));
        end
        spkinds(Nspkinds,:)=spkvals_sv(ind1+i-1,:);
        ind1=ind2;ind2=ind1+1;
    else
        ind2=ind2+1;
    end
    if (ind1>=size(spkvals_sv,1))
        break;
    end
end
spkinds = spkinds(1:Nspkinds,:);

while (1)
    pp=find(diff(spkinds(:,1))<SameSpkWin);
    if (length(pp)==0)
        break;
    else
        for ii=length(pp):-1:1
            if (PKVAL==0)
                [y,i]=min(spkinds(pp(ii)+[0:1],2));
                spkinds(pp(ii)+i-1,:)=[];
            else
                [y,i]=max(spkinds(pp(ii)+[0:1],2));
                spkinds(pp(ii)+i-1,:)=[];
            end
        end
    end
end
spkinds=spkinds(:,1);
return;
