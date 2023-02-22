% ANOVA-violin helper script

[P,ANOVATAB,STATS] = anova1(testme);
close
myfigure
myviolin(testme,'medc',[])
xticks([1 2 3])
xticklabels({'Before','During','After'})
xlim([0.5 3.5])
makefigpretty
ylabel(ystr)
title(sprintf('F(%i,%i) = %1.3f, p = %1.3f',...
    ANOVATAB{2,3},ANOVATAB{3,3},ANOVATAB{2,5},P),'fontsize',18)