colordef white 
figure

ind=find(spkind((spksinclust),2)<5);
ind2=find(spkind(:,2)<5);
inspikes=spkamp(spksinclust(ind),:);
outspikes=spkamp(ind2,:);











for ii = 1:size(comb,1)
		
        scathandle(ii)=subplot(3,2,ii);
		%subplot(NCOL,ceil(size(comb,1)/NCOL),ii);
		plot(outspikes(:,comb(ii,1)),outspikes(:,comb(ii,2)),...
                                           'k.','MarkerSize',1);
                                        hold on;
		%Add a different color for the new trial
        plot(inspikes(:, comb(ii,1)),inspikes(:,comb(ii,2)),...
                                            'r.','MarkerSize',1);
        axis off;
                                            
        title([num2str(comb(ii,1)),' vs.' ,num2str(comb(ii,2))],'FontSize',14);
end


%%random code.
for i=1:length(pp)
    ind=find(sepmat(:,1)==pp(i));
    if(ind)
        sepmat(ind,2)=0;
    end
end

ind=find(sepmat(:,2)==1);
pp=sepmat(ind,1)

%This code is used for intersections

intmat(:,1)=zeros(length(spkind),1);
intmat(spksin{1},1)=1;





figure
for ii = 3:Nusechans+2
	ax(ii)=subplot(Nusechans+1,1,ii-1);
	hold on;grid on;
	%plot(spkt(pp),dat(spki,usechans(ii-1)),[clr,'o']);
	plot(0:1/32000:((length(data)-1)/fs),data(:,ii))
    indices=find(spkind(spksin{1},2)==5);
    x=[spkind(spksin{1}(indices,1))/fs spkind(spksin{1}(indices,1))/fs]
    y=[ones(length(x),1)-5000 ones(length(x),1)-7000]
    plot(x',y','LineWidth',2,'Color','r')
    
    
   % if(indices)
    %plot(spkind(spksin{1}(indices,1))/fs,-7000,'r.');
    %end
    
    %indices=find(spkind(spksin{3},2)==75);
    %plot(spkind(spksin{3}(indices,1))/fs,-5000,'g.');
end

linkaxes(ax);
axis([4 4.8 -5001 2500])



%This bit of code will be used to replace trials which were not assigned
%correctly...

%ind=find(spkind(:,3)==0)
%list=unique(spkind(ind,2));
for i=1:length(fnm)
    fn=fnm{i}
    rd=readrecf(fn);
    for ii=1:length(outfilefilter) 
        if (strcmp(rd.outfile,outfilefilter{ii}))
            inputnum=ii;
            break;
        end
    end
    ind=find(spkind(:,2)==i);
    spkind(ind,3)=inputnum;
end


ind=find(spkind(spksin{1},3)==3);
inmatrix=spkind(spksin{1}(ind),1:2);

%this will go through data, and assign trial values correctly*/


for stim=1:10
    meanvec(stim)=mean(spikedist{stim})
    errvec(stim)=std(spikedist{stim})/sqrt(length(spikedist{stim}))
end


%change the way of doing this instead of finding all the indices, 


maxtrial=max(spkind(:,2))


for clustnum=1
    
    for stim=[1 3 6]
        ind=find(spkind(spksin{clustnum},3)==stim);
        inmatrix{clustnum,stim}=spkind(spksin{clustnum}(ind),1:2)
    end
end

figure
subplot (2,1,1);
    plotrasters(inmatrix{1,1});
subplot (2,1,2)
    plotrasters(inmatrix{1,2});

figure
ifn=2;
totclusts=1
colorlist={'r' 'b'}
for ii = 3:Nusechans+2
     for jj=1:totclusts
                        
                        ax(1)=subplot(6,1,2)
                        [sm,sp,t,f]=evsmooth(data(:,song_chan),fs,100);
                        imagesc(t,f,log(abs(sp)));syn;ylim([0,1e4]);
                        
                        ax(ii)=subplot(Nusechans+2,1,ii);
                        hold on;grid on;
                        %plot(spkt(pp),dat(spki,usechans(ii-1)),[clr,'o']);
                        plot(0:1/32000:((length(data)-1)/fs),data(:,ii),'k', 'Linewidth', 1)
                        %indices=find(spkind(spksinclust{jj},2)==ifn);
                        %indices=find(spkind(:,2)==ifn);
                        indices2=find(spkind(spksin{jj},2)==ifn);
                        xrow=[spkind(spksin{jj}(indices2),1)/fs spkind(spksin{jj}(indices2),1)/fs];
                        yrow=[ones(length(xrow),1)-10000 ones(length(xrow),1)-14000];
                        %xrow=[6.95 6.95;7.2 7.2];
                        %yrow=[-7000 7000;-7000 7000];
                        
                        plot(xrow',yrow','LineWidth',1,'Color',colorlist{jj})
                        
                        %axis off;
                    end    
end

                



