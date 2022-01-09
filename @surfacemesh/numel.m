function N = numel(dom)
%NUMEL   Number of degrees of freedom in a SURFACEMESH.
%   N = NUMEL(DOM) returns the total number of degrees of freedom in the
%   SURFACEMESH object DOM. If DOM is an array of objects, then N is the
%   number of objects in the array.

numObj = builtin('numel', dom);
if ( numObj > 1 )
    % This is an array of objects. Don't overload.
    N = numObj;
    return
end

N = 0;
for k = 1:numel(dom.x)
    N = N + numel(dom.x{k});
end

end
