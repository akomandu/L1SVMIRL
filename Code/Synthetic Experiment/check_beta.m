function [flag] = check_beta(P,A,S,gamma,beta_tar)
I1 = pinv(eye(S)-gamma*P{1});
Fcol = [];
for a=2:A; F{a}=(P{1}-P{a})*I1; end
for a=2:A; Fcol=[Fcol;F{a}]; end
try
    cvx_begin
    variables R(S)
    % minimize( 1/2*R.'*R + C * sum((-F{2}-F{3}-F{4}-F{5})*R) )
    minimize( norm(R,1))
    % minimize( 1/2*R.'*R -rho)
    subject to
    Fcol*R >= ones(S*(A-1),1)
    
    
    
    
    
    cvx_end
    beta = 1./norm(R,1);
catch
    flag = false;
    return
end

flag = beta>=beta_tar;
end