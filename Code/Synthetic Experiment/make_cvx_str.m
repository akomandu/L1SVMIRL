function cvx_str = make_cvx_str(A)
cvx_str_temp = 'cvx_begin\nvariables R(S) rho\nminimize((CC/2)*R.''*R -rho)\nsubject to\n';
for k = 2:A
    cvx_str_temp = [cvx_str_temp, '\tF{',num2str(k),'}*R >= rho*ones(S,1)\n'];
end
cvx_str_temp = [cvx_str_temp,'cvx_end'];
cvx_str = sprintf(cvx_str_temp);