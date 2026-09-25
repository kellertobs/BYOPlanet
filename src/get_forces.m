%--------------------------------------------------------------
% Gravitational force at new positions
%--------------------------------------------------------------

function [Fj] = get_forces(M,X,N)

DX = reshape(X,[N 1 3]) - reshape(X,[1 N 3]);

D2 = sum(DX.^2,3);
D  = sqrt(D2);

D3 = D2 .* D;

W = (M .* M.') ./ D3;

% Remove self-interaction
W(1:N+1:end) = 0;

Fj = -squeeze(sum(W .* DX,2));

end