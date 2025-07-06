function [c, ceq] = nonlcon(x)
    p = x(1);
    R = x(2);

    % Nonlinear constraint condition
    c = [-p; p-1; -R; -p]; %0<p<1,R>0
    ceq = [];
end
