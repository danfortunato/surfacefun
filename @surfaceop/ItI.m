function ItI = ItI(S, eta)

if ( ~isInitialized(S) )
    error('The object has not been initialized.');
end

if ( ~isBuilt(S) )
    build(S);
end

switch S.method
    case 'DtN'
        DtN = S.patches{1}.BtB;
        I = eye(size(DtN));
        ItI = (DtN-1i*eta*I)/(DtN+1i*eta*I);
    case 'ItI'
        ItI = S.patches{1}.BtB;
end

end
