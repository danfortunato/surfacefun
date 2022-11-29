function latexfonts
% LATEXFONTS plots several fonts avalible in the LaTeX interpreter for a sample string.
%
% 	See also: tex,  texlabel

% Author:   Yaroslav Don
% Created:  Jun-2008
% Review:   08-Jul-2009
% Version:  1.0.2
% Copyright 2008-2010

% set figure
figure('Units','normalized','Position',[0 0 1 1]);
%

% define various fonts and labels
family     = {'\rmfamily', '\sffamily',  '\ttfamily'};
familyname = {'Roman',     'Sans Serif', 'Typewriter'};
%
shape      = {'\upshape', '\slshape', '\itshape', '\scshape'};
shapename  = {'Upright',  'Slanted',  'Italics',  'Sm Caps'};
%
series     = {'\mdseries', '\bfseries'};
seriesname = {'Medium',    'Bold'};

% other definitions and preallocations
[lfa, lsh, lse]  = deal(length(family), length(shape), length(series));
h1               = zeros(lfa,lsh,lse);
h2               = zeros(lfa,lse);
h3               = zeros(lsh,lse);
%
txt              = 'WALT bla fi f{}i ff';
[x0, y0, dx, dy] = deal(0.10, 1.00, 0.35, 0.1);

regular_text = @(fa,sh,se) sprintf('%s%s%s %s',               family{fa},shape{sh},series{se},txt);
title_text   = @(se,fa)    sprintf('\\underline{%s%s %s %s}', series{se},family{fa},seriesname{se},familyname{fa});
label_text   = @(sh)       sprintf('%s %s:',                  shape{sh},shapename{sh});

% main text insertion
for se = 1:lse,
	for sh = 1:lsh,
		for fa = 1:lfa,
			lbk = (se-1)*(lsh+2);
			h1(fa,sh,se) = text(x0+(fa-1)*dx, y0-(sh+lbk)*dy,  regular_text(fa,sh,se));
			if sh == 1, % write titles
				h2(fa,se) = text(x0+(fa-1)*dx, y0-lbk*dy,       title_text(se,fa));
			end
			if fa == 1, % write titles
				h3(sh,se) = text(-x0,          y0-(sh+lbk)*dy,  label_text(sh));
			end
		end
	end
end

% set LaTeX interpreter
set([h1(:);h2(:);h3(:);], 'Interpreter', 'LaTeX', 'FontSize', 20);
axis off;
