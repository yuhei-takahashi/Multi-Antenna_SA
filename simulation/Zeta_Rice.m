function ret=Zeta_Rice(n,K,A)

ret=(1+K)./(exp(K).*factorial(n-1).^2.*A).*(K.*(1+K)./A).^(n-1);

end
