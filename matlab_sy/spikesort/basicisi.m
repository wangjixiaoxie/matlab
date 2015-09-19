function isitimes = basicisi(spktimes,maxisi,binsize)

%binsize = 0.001; % in s

isitimes = [];
goodisi = 0;
badisi = 0;

for spk = 2:length(spktimes)
  
  isi = spktimes(spk)-spktimes(spk-1);
  
  if (isi <= maxisi)
    goodisi = goodisi + 1;
    isitimes = [isitimes isi];
  else
    badisi = badisi + 1;
  end
  
end

%figure;
nobins = ceil(maxisi / binsize);
[h,tvec] = hist(isitimes,nobins);
h = h / (goodisi+badisi);
bar(tvec,h);
axis([0 tvec(end) 0 max(h)+0.2*max(h)]);
