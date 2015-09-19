figure

for i=1:5
    ffout{i}=fft(pshift{i}(x(1):x(2)));
    
    subplot(5,1,i);
    
 
    edges=0:(fs/2/256):fs/2
    plot(edges,abs(ffout{1}(end:-1:256)),'b')
    hold on;
    plot(edges,abs(ffout{i}(end:-1:256)),'r')
    hold on;
    axis([0 22000 0 10])
    
    
end

figure
for i=1:5
    %ffout{i}=fft(pshift{i}(x(1):x(2)));
    
    subplot(5,1,i);
    
    plot(edges,angle(ffout{1}(end:-1:256)),'b')
    hold on;
    plot(edges,angle(ffout{i}(end:-1:256)),'r')
    hold on;
    axis([0 22000 -pi pi])
    
    
end