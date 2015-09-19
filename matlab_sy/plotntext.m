function plotntext(day,synstruct, height)
    for ii=1:length(day)
        nval=synstruct{ii}.ntrans;
        strval=[ num2str(nval)]
        text(day(ii)-day(1), height,strval)
    end