function ensure_folder(p)
if ~exist(p,'dir'), mkdir(p); end
end
