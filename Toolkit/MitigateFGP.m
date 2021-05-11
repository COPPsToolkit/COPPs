function ZZ = MitigateFGP(H,TYPE,P)
%--------------------------------------------------------------------------
% O. de Groot and F. Mazelis (2020)
% "Mitigating the Forward Guidance Puzzle" (2020)
% ECB working paper #2426
%
% NOTES:
%   Let x = [x_{1} ; ... ; x_{H}] be a (m.nx*H x 1) vector of states
%   and a = [a_{1} ; ... ; a_{H}] be a (m.ns*H x 1) vector of antic shocks
%
% IRFs:
%   P x = Q a, where, for H = 4:
%
%   P = [  I  0  0  0 ]   QT =[ PHI  F*PHI  F^2*PHI  F^3*PHI ]   Q = ZZ .* QT
%       [ -B  I  0  0 ]       [  0    PHI    F*PHI   F^2*PHI ]
%       [  0 -B  I  0 ]       [  0     0      PHI     F*PHI  ]
%       [  0  0 -B  I ]       [  0     0       0       PHI   ]
%
%--------------------------------------------------------------------------

ZZ   = zeros(H,H);
for i = 1:H
    for j = i:H
        if strcmp(TYPE,'I') % Inattention
            ZZ(i,j) = P.alpha^(j-i);
        elseif strcmp(TYPE,'II') % Credibility
            if j==i
                ZZ(i,j) = 1;
            else
                ZZ(i,j) = ZZ(i,j-1)*P.alpha^(j-i);
            end
        elseif strcmp(TYPE,'III') % Finite planning horizon
            if j>i+P.N_planning
                ZZ(i,j) = 0;
            else
                ZZ(i,j) = 1;
            end
        elseif strcmp(TYPE,'IV') % Learning
            if j==i
                ZZ(i,j) = 1;
            else
                ZZ(i,j) = ZZ(i,j-1)/(1+exp(-P.beta1*(j-P.beta2)));
            end
        end
    end
end