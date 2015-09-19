function jcplotblack(data)
map(1,:)=[1 1 1];
for i=2:30
    map(i,:)=[0 0 0];
end
figure; imagesc(data)
colormap(map)