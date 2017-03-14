function lt_test

x = 1:100;
y = 1:100;

plot(x,y, 'ButtonDownFcn', @lineCallback)
end


function lineCallback(src, event)
set(src, 'Color', 'red');
end
