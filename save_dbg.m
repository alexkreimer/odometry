function save_dbg(filename)
set(gcf,'units','normalized','outerposition',[0 0 1 1])
set(gcf,'PaperPositionMode','auto');
set(gcf,'InvertHardcopy','off');
%print(filename, '-dpng', '-r864');
print(filename, '-dpng', '-r300');
end
