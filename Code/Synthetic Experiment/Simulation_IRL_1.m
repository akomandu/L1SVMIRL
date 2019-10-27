A = 4;
S = 4;
n1 = 50;
n2 =100;
gamma = 0.1;
m = 100;

RNg = {};
RhatNg = {};
RMe = {};
RhatMe = {};
C = 5;
lambda = 0.1;
for kk = 1:m
    [P,Phat] = f_make_data(A,S,n1,n2);
    %% Ng solutions
I1 = inv(eye(S)-gamma*P{1});
clear F
for a=2:A; F{a}=(P{1}-P{a})*I1; end

cvx_begin
variables R(S) rho
% minimize( 1/2*R.'*R + C * sum((-F{2}-F{3}-F{4}-F{5})*R) )
minimize( lambda*norm(R,1) + sum(max(max(-F{2}*R,-F{3}*R),-F{4}*R)) )
% minimize( 1/2*R.'*R -rho)
subject to
  F{2}*R >= zeros(S,1)
  F{3}*R >= zeros(S,1)
  F{4}*R >= zeros(S,1)
  
  R>= zeros(S,1)-1
  R<= zeros(S,1)+1
cvx_end
RNg{kk} = R;

% hat problem
I1 = inv(eye(S)-gamma*Phat{1});
clear F
for a=2:A; F{a}=(Phat{1}-Phat{a})*I1; end

cvx_begin
variables R(S) rho
% minimize( 1/2*R.'*R + C * sum((-F{2}-F{3}-F{4}-F{5})*R) )
minimize( lambda*norm(R,1) + sum(max(max(-F{2}*R,-F{3}*R),-F{4}*R)) )
% minimize( 1/2*R.'*R -rho)
subject to
  F{2}*R >= zeros(S,1)
  F{3}*R >= zeros(S,1)
  F{4}*R >= zeros(S,1)
  
  R>= zeros(S,1)-1
  R<= zeros(S,1)+1
cvx_end
RhatNg{kk} = R;


%% My method 
I1 = inv(eye(S)-gamma*P{1});
clear F
for a=2:A; F{a}=(P{1}-P{a})*I1; end

cvx_begin
variables R(S) rho
% minimize( 1/2*R.'*R + C * sum((-F{2}-F{3}-F{4}-F{5})*R) )
% minimize( lambda*norm(R,1) + sum(max(max(-F{2}*R,-F{3}*R),-F{4}*R)) )
minimize( 1*R.'*R -rho)
subject to
  F{2}*R >= rho*ones(S,1)
  F{3}*R >= rho*ones(S,1)
  F{4}*R >= rho*ones(S,1)
  
cvx_end
RMe{kk} = R;
% hat problem

I1 = inv(eye(S)-gamma*Phat{1});
clear F
for a=2:A; F{a}=(Phat{1}-Phat{a})*I1; end

cvx_begin
variables R(S) rho
% minimize( 1/2*R.'*R + C * sum((-F{2}-F{3}-F{4}-F{5})*R) )
% minimize( lambda*norm(R,1) + sum(max(max(-F{2}*R,-F{3}*R),-F{4}*R)) )
minimize( 1*R.'*R -rho)
subject to
  F{2}*R >= rho*ones(S,1)
  F{3}*R >= rho*ones(S,1)
  F{4}*R >= rho*ones(S,1)
  
cvx_end
RhatMe{kk} = R;
end


for kk = 1:m
    RdiffNg(kk) = max(abs(RNg{kk} - RhatNg{kk}));
    RdiffMe(kk) = max(abs(RMe{kk} - RhatMe{kk}));
end

