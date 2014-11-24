"""
algorithm for simple fast biharmonic equation
"""

function = simplefastbiharmonic(F)
    m = length(F);
    L = 1./(m+1);
    p = pi*L*(1:m);
    sig = sin(p/2.).^2;
    S  = sin(p*(1:m));
    G = S*F*S;
    X = (L^2)*G./(4*(sig*ones(1,m)+ones(m,1)*sig )      ).^2);
    
    U = zeros(m+2,m+2)
    U(2:m+1,2:m+1) = S*X*S;