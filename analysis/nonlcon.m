function [c, ceq] = nonlcon(x)
    p = x(1);
    R = x(2);

    % 非線形制約条件
    c = [-p; p-1; -R; -p]; %0<p<1,R>0
    ceq = [];
end
