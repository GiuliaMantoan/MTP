function styleAxis(ax, dates, step_yrs)
%STYLEAXIS  Zero line, grid, annual tick labels.
    yline(ax, 0,'k-','LineWidth',0.75,'HandleVisibility','off');
    grid(ax,'on');
    ticks = datetime(year(dates(1)),1,1):calyears(step_yrs):dates(end);
    ax.XTick = ticks;
    ax.XAxis.TickLabelFormat = 'yyyy';
end