function [spk_tms,ind]=analbird(upper,thr,snd,ch,samplerate);

% loading prelim
%load data_0715;
%snd = x371(1,:);
%ch = x371(2:4,:);

%clear x03 x152 x231 x28 x91 ;
%pack;

% define vars
%upper = TH(2);
%thr = TH(1);
%samplerate = 32000;

tvec = 0:length(snd)-1;
tvec = tvec * (1/samplerate);

% get spike times (stupid method)
for c = 1:1% all channels

  [val,ind] = findpeaks(ch(c,:),thr,upper);
  figure(1);
  plot(tvec,ch(c,:));
  hold on;
  spk_tms = tvec(ind);
  plot(tvec(ind),val,'.r');
  
%convolve trace
%cc = conv2(ch(c,:),mean(waves(1:10,:)),'same');
%figure(10);
%plot(tvec,cc);

%%%%%%%%%%%%%%%%%
% START CUTTING %
%%%%%%%%%%%%%%%%%

vt = 20;
happy = 0;
while (~happy)

  % arrange waveforms & get features
  figure(2);
  hold on;
  ww = 50; % waveform size
  offset = 10;
  waves = zeros(length(spk_tms),ww);
  energies = zeros(length(spk_tms),1);
  amps = zeros(length(spk_tms),1);
  vts = zeros(length(spk_tms),1);
  slopes = zeros(length(spk_tms),1);
  widths = zeros(length(spk_tms),1);
  tpeaks = zeros(length(spk_tms),1);
  baseline = mean(ch(c,:));

  for i = 1:length(ind)
    waves(i,:) = ch(c,ind(i)-ww/2+1+offset:ind(i)+ww/2+offset);
    waves(i,:) = waves(i,:) - baseline;
    dwave = diff(waves(1,:));
    energies(i) = sum(waves(i,:).^2);
    vts(i) = waves(i,vt);
    [maxamp,maxind] = max(waves(i,:));
    [minamp,minind] = min(waves(i,:));
    amps(i) = maxamp - minamp;
    widths(i) = abs(minind-maxind)+randn(1)/10;
    slopes(i) = dwave(vt)+randn(1)/1000;
    plot(waves(i,:));
  end
  
  % plot waveforms
  meanwave = mean(waves);
  stdwave = std(waves);
  plot(meanwave,'r','LineWidth',2);
  plot(meanwave+stdwave,'r--');
  plot(meanwave-stdwave,'r--');
  % plot some features
  hold off;
  figure(3);
  plot(energies,vts,'.');
  xlabel('energies');
  ylabel('vts');
  figure(4);
  plot(amps,vts,'.');
  xlabel('amps');
  ylabel('vts');
  figure(5);
  plot(energies,widths,'.');
  xlabel('energies');
  ylabel('widths');
  figure(6);
  plot(amps,slopes,'.');
  xlabel('amps');
  ylabel('slopes');
  
  % plot acorr
  figure;
  [corr,times] = spk_acorr2(spk_tms,1,300);
  bar(times,corr);  
  title(sprintf('unfiltered trace autocorrelation, nspks = %d',length(spk_tms)));
  v = axis;
  axis([times(1) times(end) v(3) v(4)]);
  xlabel('time lag (ms)');
  ylabel('spike probability');
  
  f = input('Happy? (y/n): ','s');
  
  switch(f)
    case('n')
      
       
      newvt = input('Select new Vt? (y/n): ','s');
      
      if(strcmp(newvt,'y'))
	
	figure(2);
	coords = ginput(1);
	vt = round(coords(1));
	continue;
	
      end
      
      fno = input('Cut in which figure? ');
      figure(fno);
      h = get(gca, 'xlabel');
      xval = eval(get(h,'String'));
      h = get(gca, 'ylabel');
      yval = eval(get(h,'String'));
      
      % zone drawing
      disp('Select to points to form cutting-out rectangle.');
      coords = ginput(2);
      miny = min(coords(:,2));
      maxy = max(coords(:,2));
      minx = min(coords(:,1));
      maxx = max(coords(:,1));
      
      hold on;
      
      plot([minx minx],[miny maxy],'k');
      plot([minx maxx],[miny miny],'k');
      plot([minx maxx],[maxy maxy],'k');
      plot([maxx maxx],[miny maxy],'k');
      
      % throw out selected stuff
      
      keep = ones(1,length(spk_tms));
      delcount = 0;
      for spk = 1:length(spk_tms)
	if ((xval(spk) >= minx) & (xval(spk) <= maxx) & (yval(spk) >= ...
	      miny) & (yval(spk) <= maxy))
	  keep(spk) = 0;
	  delcount = delcount + 1;
	end
      end
      
      disp(sprintf('Deleting %d spikes...',delcount));
      
      spk_tms = spk_tms(find(keep == 1));
      ind = ind(find(keep == 1));
      energies = energies(find(keep == 1));;
      vts = vts(find(keep == 1));
      %vts2 = vts2(find(keep == 1));
      slopes = slopes(find(keep == 1));
      amps = amps(find(keep == 1));
      widths = widths(find(keep == 1));
      waves = waves(find(keep == 1),:);
      
      close all;
    case('y')
      happy = 1;
  end % of switch

end % of while loop for cutting iterations

cut_spikes{c} = spk_tms;
isitimes = basicisi(cut_spikes{c},200,2);
figure;
hold on;
for i=2:length(isitimes)
  plot(isitimes(i),isitimes(i-1),'.');
end
title(sprintf('ISI vs. previous ISI plot for channel %d',c));
xlabel('ISI');
ylabel('preceding ISI');

end % of loop over traces




