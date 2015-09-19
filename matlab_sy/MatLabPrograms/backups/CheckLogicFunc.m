function Match=CheckLogicFunc(SingleND,MatchVals);
%
%
NTemplates=size(SingleND.Templ,2);

% string with the arbitrary logic function
LogicFunc=SingleND.CntLog;
for itempl=1:NTemplates
    
    % the variable name of this counter in the logic function
    % by default the first counter is 'a' the second 'b' etc.
    VarName=SingleND.CntRng(itempl).VarName;
    
    % make a variable in matlab with the VarName with the corresponding value from MatchVals 
    eval([VarName,'=',num2str(MatchVals(itempl)),';']);
end

Match=eval(LogicFunc);
Match=(Match>=1);
return;
