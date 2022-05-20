%
% This is the setupfig.m used only for mySST
%
% 2011-09-09 Hau-tieng Wu
%

    % setup the handle of the final display figure
for QQQQ=1:Yno*Xno
    eval(['h',num2str(QQQQ),'=subplot(',num2str(Yno),',',num2str(Xno),',',num2str(QQQQ),');']);
end
