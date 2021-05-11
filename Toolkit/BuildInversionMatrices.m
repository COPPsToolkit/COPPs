function [Inverter,P,Q,R,C] = BuildInversionMatrices(m,params)
%--------------------------------------------------------------------------
%
% NOTES:
%   Let x = [x_{1} ; ... ; x_{H}] be a (m.nx*H x 1) vector of states
%   and u = [u_{1} ; ... ; u_{H}] be a (m.ns*H x 1) vector of unant shocks
%   and a = [a_{1} ; ... ; a_{H}] be a (m.ns*H x 1) vector of antic shocks
%
%   Then x = P\C*x_{0} + P\Q*a + P\R*u where, for H = 4:
%
%    P = [  I  0  0  0 ]   C = [ B ]   Q = [ PHI  F*PHI  F^2*PHI  F^3*PHI ]
%        [ -B  I  0  0 ]       [ 0 ]       [  0    PHI    F*PHI   F^2*PHI ]
%        [  0 -B  I  0 ]       [ 0 ]       [  0     0      PHI     F*PHI  ]
%        [  0  0 -B  I ]       [ 0 ]       [  0     0       0       PHI   ]
%
%    R = [ PHI   0    0    0  ]
%        [  0   PHI   0    0  ]
%        [  0    0   PHI   0  ]
%        [  0    0    0   PHI ]
%
% Acknowledgements: Modified from code provided by J. Barrdear
%--------------------------------------------------------------------------

AnticShocksToUse = {m.PolicyInstrumentsAndShocks{:,2}}';
UnantShocksToUse = m.ShockNames;

nv  = length({m.PolicyInstrumentsAndShocks{:,1}}');
nsa = length(AnticShocksToUse);
nsu = length(UnantShocksToUse);
    

    P = eye(m.nx*params.T_loss);
    for i = 1:(params.T_loss-1)
        P(m.nx*i+1:m.nx*(i+1),m.nx*(i-1)+1:m.nx*i) = -m.B;
    end

% Mitigating the Forward Guidance Puzzle, O. de Groot and F. Mazelis
m.ZZ = MitigateFGP(params.T_loss,params.TYPE,params);

    Q = zeros(m.nx*params.T_loss,m.ns*params.T_loss);
    for q = 1:params.T_loss
        for j = q:params.T_loss
            %Q(m.nx*(q-1)+1:m.nx*q,m.ns*(j-1)+1:m.ns*j) = m.F^(j-q)*m.PHI;
            
            % Mitigating the Forward Guidance Puzzle, O. de Groot and F. Mazelis
            Q(m.nx*(q-1)+1:m.nx*q,m.ns*(j-1)+1:m.ns*j) = m.F^(j-q)*m.PHI*m.ZZ(q,j);
        end
    end

    R = kron(eye(params.T_loss),m.PHI);
    C = kron(eye(params.T_loss,1),m.B);
    
    VarSelectorVec = zeros(1,m.nx);
    for i = 1:nv
        idx = find(strcmp({m.PolicyInstrumentsAndShocks{i,1}},m.VarNames));
        VarSelectorVec(idx) = 1;
    end
    VarSelector = diag(VarSelectorVec);
    VarSelector = VarSelector(logical(VarSelectorVec),:);
    VarSelector = kron(eye(params.T_loss),VarSelector);

    UnantSelectorVec = zeros(1,m.ns);
    for i = 1:nsu
        idx = find(strcmp(UnantShocksToUse{i},m.ShockNames));
        UnantSelectorVec(idx) = 1;
    end
    UnantSelector = diag(UnantSelectorVec);
    UnantSelector = UnantSelector(:,logical(UnantSelectorVec));
    UnantSelector = kron(eye(params.T_loss),UnantSelector);

    
    AnticSelectorVec = zeros(1,m.ns);
    for i = 1:nsa
        idx = find(strcmp(AnticShocksToUse{i},m.ShockNames));
        AnticSelectorVec(idx) = 1;
    end
    AnticSelector = diag(AnticSelectorVec);
    AnticSelector = AnticSelector(:,logical(AnticSelectorVec));
    AnticSelector = kron(eye(params.T_loss),AnticSelector);
    
    
    W = [VarSelector*(P\R)*UnantSelector , VarSelector*(P\Q)*AnticSelector];
    
    ShockSelector = blkdiag(UnantSelector , AnticSelector);
    
    if nsu + nsa >= nv
        Inverter = ShockSelector*W'/(W*W');
    else
        Inverter = ShockSelector*((W'*W)\W);
    end
end