%ST.days N rows (number of days)
%ST.colors V rows (number of transitions)
%ST.prob N rows and V columns
%ST.err N rows and V columns

function [outmat,ss]=plotseqtm(st)

for ii=1:length(st.prob(1,:))
    plot(st.days,100*st.prob(:,ii),st.col{ii});
    hold on;
    plot(st.days,100*st.prob(:,ii),'Color',st.col{ii},'Marker','o');

%     for ii=1:length(crdt.days)
%         text(st.days(ii),80, num2str(st.targn(ii)),'Color','k');
%         hold on;
%         text(st.days(ii),40,num2str(st.alln(ii)-st.targn(ii)),'Color','r')
%     end
end 
    box off;