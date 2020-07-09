function xdot = fcn_K_ode(t,x,K)

% xdot=zeros(size(K,1),1);
xdot = K*x;

