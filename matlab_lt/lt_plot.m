function hfig=lt_plot(X,Y,Modifiers)
%% 11/30/14 - my version default of plot
% X, Y normal
% Modifiers = {'Color','k','LineWidth',2,'Errors',Yerr}; (enter as cell array, has to be even number)

%% Default params
plot_errors=0;

if nargin==2;
    % is last arg cell? then am missing X
    if iscell(Y);
        % reasign names
        Modifiers=Y;
        Y=X;
    else
        % then have X, Y, but no modifier
%            Modifiers={'LineStyle','none','Marker','o','Color','k','MarkerSize',6};
                      Modifiers={'LineStyle','none','Marker','o','Color','k','MarkerFaceColor','k','MarkerSize',5};

    end
elseif nargin==1;
    % then we just have a Y
    Y=X;
    X=1:length(Y);   
    Modifiers={'LineStyle','none','Marker','o','Color','k','MarkerFaceColor','k','MarkerSize',5};   
%     Modifiers={'LineStyle','none','Marker','o','Color','k','MarkerSize',6};   

end


% Assign default things to modifiers (apart from color) - linestyle and
% marker type and size, etc (i.e. things that I would not specify in my
% modifiers). first check if already in modifers, if not, then speciify
if ~any(strcmp(Modifiers,'LineStyle'));
    Modifiers=[Modifiers, 'LineStyle','none'];
end

if ~any(strcmp(Modifiers,'Marker'));
    Modifiers=[Modifiers, 'Marker','o'];
end
    
if ~any(strcmp(Modifiers,'Color'));
    Modifiers=[Modifiers, 'Color','k'];
end

if ~any(strcmp(Modifiers,'MarkerFaceColor')) && any(strcmp(Modifiers,'Color')); % set marker face color to equal color
    Modifiers=[Modifiers, 'MarkerFaceColor',Modifiers{find(strcmp(Modifiers,'Color'))+1}];
end

if ~any(strcmp(Modifiers,'MarkerSize'));
    Modifiers=[Modifiers, 'MarkerSize',6];
end

if any(strcmp(Modifiers,'Errors'));
    Yerr=Modifiers{find(strcmp(Modifiers,'Errors'))+1};
    plot_errors=1;
end




%% RUN

% then plot
if plot_errors==0;
hfig=plot(X,Y);

else
    if length(X)==1;
        X=ones(length(Y), 1)*X;
    end
       
    hfig=errorbar(X, Y, Yerr);
    errorbar_tick(hfig, 1000)
end

% apply modifiers
for i=1:length(Modifiers)/2;
    if any(strcmp(Modifiers{2*i-1}, 'Errors'));
        continue;
    end
    set(hfig,Modifiers{2*i-1},Modifiers{2*i});
end

% format color and font size
lt_plot_format;

set(hfig,'LineWidth',2)


