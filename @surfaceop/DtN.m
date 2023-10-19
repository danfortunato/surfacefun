function DtN = DtN(S)

if ( ~isInitialized(S) )
    error('The object has not been initialized.');
end

if ( ~isBuilt(S) )
    build(S);
end

switch S.method
    case 'DtN'
        DtN = S.patches{1}.BtB;
    case 'ItI'
        ItI = S.patches{1}.BtB;
        I = eye(size(ItI));
        % TODO: Why is the prefactor not -1i*eta?
        DtN = 1i/S.eta*(ItI-I)\(ItI+I);
end

end
