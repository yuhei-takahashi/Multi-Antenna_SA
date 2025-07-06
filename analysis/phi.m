function ret=phi(N,c,n,K,A)

    J=0;
    for j=1:N
        zeta_j=Zeta_Rice(j,K,A);
        J1=zeta_j.*gamma(j).*c.^(-j);
        J=J+J1;
    end

        zeta_n=Zeta_Rice(n,K,A);
        ret=zeta_n./J;
end
