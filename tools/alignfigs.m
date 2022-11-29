function alignfigs(varargin)
%alignfigs   Arrange all open figure windows nicely.
%   ALIGNFIGS() arranges all the open current figure windows to fit nicely on the
%   screen in a grid. If they won't fit, then the size is reduced uniformly.
%
%   ALIGNFIGS(N) specifies that there should be N columns.
%
%   Note that if not all windows have the same size, they will be adjusted to be
%   so (to the smallest open window).
%
%   Examples:
%       close all, for k = 1:8, figure, end, alignfigs()
%       close all, for k = 1:15, figure, end, alignfigs(5)

% Nick Hale - Dec 2014.

h = findobj('Type', 'figure');
nh = numel(h);

% Windows are taller than they claim to be?
fudge1 = 182; % Also the top on is different...
fudgeN = 158;

% Use these values if removing the toolbar.
% fudge1 = 130;
% fudgeN = 100;

if ( nh < 2 )
    return
end

% Sort into the right order:
[~, idx] = sort(cell2mat(get(h, 'Number')));
h = h(idx);

% Get the sizes:
for k = 1:numel(h)
    p(k,:) = get(h(k), 'position');
%     set(h(k), 'toolbar', 'none')
end
% Adjust to be the same size:
p(:,3:4) = repmat(min(p(:,3:4)),nh, 1);

% Get the screen size:
ss = get(0, 'ScreenSize');
ss = ss(3:4);

% Determine the grid size:
nxMax = floor(ss(1)/p(1,3));
if ( nargin == 0 )
    nx = nxMax;
else
    nx = varargin{1};
end

% Screen isn't large enough! Shrink the figures:
ny = floor( (ss(2)+fudgeN) / (p(1,4)+fudgeN) );
if ( nh > nx*ny || nx > nxMax )
    for k = 1:nh
        set(h(k), 'position', [p(k,1:2), .9*p(k,3:4)]);
    end
    alignfigs(varargin{:});
    return
end

% Add a margin at the left to centre:
leftMargin = (ss(1) - nx*p(1,3))/2;

% Place the first figure in the top lefT:
p(1,:) = [leftMargin ss(2) - p(1,4) p(1,3:4)];
set(h(1), 'position', p(1,:));
p(1,:) = get(h(1),'position');

% Initialise px and py:
if ( nx > 1 )
    px = p(1,1) + p(1,3);
    py = p(1,2);
else
    % Need to treat single column as a special case.
    px = leftMargin;
    py = p(1,2) - p(1,4) - fudge1;
end

% Loop over the remaining windows:
for k = 2:numel(h)
    % Set to the new px and py values
    p(k,1:2) = [px, py];
    set(h(k), 'position', p(k,:));
    if ( mod(k, nx) )
        % Move to the next column:
        px = px + p(k-1,3);
    else
        % Move to the next row. 
        px = leftMargin;
        % Update py:
        if ( k == nx )
            % First row is a special case (not sure why!)
            py = py - p(k-1,4) - fudge1;
        else
            py = py - p(k-1,4) - .5*fudgeN;
        end
    end        
end

% Bring all the windows to the front:
for k = 1:nh
    figure(h(k))
end

end
