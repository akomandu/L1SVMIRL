function cvx_str = make_cvx_str_test(A)
%Fcol is all Fai in column where each row sums to 1 in the corresponding
%Pai
cvx_str_temp = ['cvx_begin\nvariables R(S) lam(S*(',num2str(A-1),')) rho\nminimize (sum(abs(R)) -rho)\nsubject to\n'];
for k = 2:A
    cvx_str_temp = [cvx_str_temp, '\tF{',num2str(k),'}*R >= rho*ones(S,1)\n'];
end
cvx_str_temp = [cvx_str_temp, '\tR == Fcol.''*lam \n'];
cvx_str_temp = [cvx_str_temp, '\tlam >= zeros(size(lam))\n'];
 cvx_str_temp = [cvx_str_temp, '\tlam >= zeros(size(lam))\n'];
 cvx_str_temp = [cvx_str_temp, '\tsum(lam) == 1\n'];
cvx_str_temp = [cvx_str_temp,'cvx_end'];
cvx_str = sprintf(cvx_str_temp);