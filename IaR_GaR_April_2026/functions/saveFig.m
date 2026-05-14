function saveFig(fig, folder, filename, dpi)
%SAVEFIG  Save figure as PNG and close it.  Default dpi = 200.
    if nargin < 4, dpi = 200; end
    print(fig, fullfile(folder, filename), '-dpng', sprintf('-r%d',dpi));
    close(fig);
end
