function val = get(f, prop)
%GET   Get properties of a SURFACEFUN.
%   VAL = GET(F, PROP) returns the value of the property specified in the
%   string PROP from the SURFACEFUN object F. Valid entries for PROP are:
%
%      'DOMAIN'
%      'VALS'
%
%   See also SUBSREF.

% Get the properties.
switch ( prop )
    case 'domain'
        val = f.domain;
    case {'vals', 'values'}
        val = f.vals;
    otherwise
        error('SURFACEFUN:get:propName', ...
            [prop,' is not a valid surfacefun property.'])
end

end
