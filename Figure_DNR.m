function [h]=Figure_DNR(New_Fig);
% Usage : [h]=Figure_DNR(New_Fig);
% Creates a bigger figure with white background
% Either a new figure (New_Fig=1)
% Or an existing figure (New_Fig<>1)

if New_Fig==1
    h= figure;
else
    h= gcf;
end

figpos = get(h,'Position');
rect = [figpos(1)*.25, figpos(2) *.25, figpos(3)*1.5, figpos(4)*1.5];
set(h,'Position',rect);
set(h,'color','w'); % pour les bordures

