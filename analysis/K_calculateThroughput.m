function Throughput = K_calculateThroughput(ch, steadyState, p, R, L)
ch.eta_0 = 2^R-1;

% calculate the (1,1) entry
tr_throughput_matrix(1) = Dec11(p, ch, L) + Dec21(p, ch, L) + 2.0 * Dec22(p, ch, L);

% calculate the (2,1) entry
tr_throughput_matrix(2) = 2.0 * Rec_nop(p, ch, L) + Rec_pot(p, ch, L);
% tr_throughput_matrix(2) = 0; %without inter-SIC

% calculate the (3,1) entry
tr_throughput_matrix(3) = 2.0 * Dec11(p, ch, L) + 2.0 * Dec22(p, ch, L) + 2.0 * Dec21(p, ch, L);
% tr_throughput_matrix(3) = 0; %without inter-SIC

% calculate throughput
result = steadyState .* tr_throughput_matrix';
Throughput = sum(result);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% P_{\mathsf{A}} or P_{\breve{\mathsf{A}}}
function prob = P_and(ch)
temp_22=0;
N=length(ch.MG_a);
eta_0=ch.eta_0;
for n=1:N
    a_n=ch.MG_a(n);
    b_n=ch.MG_b(n);
    c_n=ch.MG_c(n);
    temp_22_k1=0;
    for k_1=0:ch.MG_b(n)-1
        temp_22_k2=0;
        for k_2=0:k_1
            sum_m_term=0;
            for m=1:N
                a_m=ch.MG_a(m);
                b_m=ch.MG_b(m);
                c_m=ch.MG_c(m);
                sum_m_term=sum_m_term+a_m*(c_m+c_n*eta_0)^(-(b_m+k_2))*gamma(b_m+k_2)*gammainc(eta_0*(c_m+c_n*eta_0),b_m+k_2,"upper");
            end
            temp_22_k2=temp_22_k2+nchoosek(k_1,k_2)*sum_m_term;
        end
        temp_22_k1=temp_22_k1+((c_n*eta_0)^k_1/factorial(k_1))*temp_22_k2;
    end
    temp_22=temp_22+a_n*c_n^(-b_n)*gamma(b_n)*exp(-c_n*eta_0)*temp_22_k1;
end
prob=temp_22;
if eta_0<1
    prob=prob-calculatePrF(ch);
end
end

% P_{\mathsf{B}} or P_{\breve{\mathsf{B}}}
function prob = P_1(ch)
N=length(ch.MG_a);
eta_0=ch.eta_0;
temp_21=0;
for n=1:N
    a_n=ch.MG_a(n);
    b_n=ch.MG_b(n);
    c_n=ch.MG_c(n);
    temp_21_k1=0;
    for k_1=0:b_n-1
        temp_21_k2=0;
        for k_2=0:k_1
            sum_m2_term=0;
            for m=1:N
                a_m=ch.MG_a(m);
                b_m=ch.MG_b(m);
                c_m=ch.MG_c(m);
                sum_m2_term=sum_m2_term+a_m*((c_m+c_n*eta_0)^(-(b_m+k_2)))*gamma(b_m+k_2)*gammainc(eta_0*(c_m+c_n*eta_0),b_m+k_2);
            end
            temp_21_k2=temp_21_k2+nchoosek(k_1,k_2)*sum_m2_term;
        end
        temp_21_k1=temp_21_k1+(((c_n*eta_0)^k_1)/factorial(k_1))*temp_21_k2;
    end
    temp_21=temp_21+a_n*c_n^(-b_n)*gamma(b_n)*exp(-c_n*eta_0)*temp_21_k1;

end

prob=temp_21;
end

% P_{\mathsf{C}} or P_{\breve{\mathsf{C}}}
function prob = P_2(ch)
N=length(ch.MG_a);
eta_0=ch.eta_0;
if eta_0>=1
    temp_col2_1th_term=0;
    temp_col2_2th_term=0;
    % the first term
    for n=1:N
        a_n=ch.MG_a(n);
        b_n=ch.MG_b(n);
        c_n=ch.MG_c(n);
        sum_k_1=0;
        for k_1=0:b_n-1
            sum_m_1=0;
            for m=1:N
                a_m=ch.MG_a(m);
                b_m=ch.MG_b(m);
                c_m=ch.MG_c(m);
                sum_m_1=sum_m_1+a_m*(c_n+c_m)^(-(b_m+k_1))*gamma(b_m+k_1)*gammainc(eta_0*(c_n+c_m),b_m+k_1,"upper");
            end
            sum_k_1=sum_k_1+(c_n^k_1/factorial(k_1))*sum_m_1;
        end
        temp_col2_1th_term=temp_col2_1th_term+a_n*c_n^(-b_n)*gamma(b_n)*sum_k_1;
    end

    % the second term
    for n=1:N
        a_n=ch.MG_a(n);
        b_n=ch.MG_b(n);
        c_n=ch.MG_c(n);
        sum_k_1=0;
        for k_1=0:b_n-1
            sum_k_2=0;
            for k_2=0:k_1
                sum_m_1=0;
                for m=1:N
                    a_m=ch.MG_a(m);
                    b_m=ch.MG_b(m);
                    c_m=ch.MG_c(m);
                    sum_m_1=sum_m_1+a_m*(c_m+c_n*eta_0)^-(b_m+k_2)*gamma(b_m+k_2)*gammainc(eta_0*(c_n*eta_0+c_m),b_m+k_2,"upper");
                end
                sum_k_2=sum_k_2+nchoosek(k_1,k_2)*sum_m_1;
            end
            sum_k_1=sum_k_1+((c_n*eta_0)^k_1/factorial(k_1))*sum_k_2;
        end
        temp_col2_2th_term=temp_col2_2th_term+a_n*c_n^(-b_n)*gamma(b_n)*exp(-c_n*eta_0)*sum_k_1;
    end

    prob=temp_col2_1th_term-temp_col2_2th_term;
else
    prob=calculatePrC(ch);
end

end

% P_{\mathsf{D}} or P_{\breve{\mathsf{D}}}
function prob = Col_1(ch)
N=length(ch.MG_a);
eta_0=ch.eta_0;
temp_col1_1th_term=0;
temp_col1_2th_term=0;

% the first term
for n=1:N
    a_n=ch.MG_a(n);
    b_n=ch.MG_b(n);
    c_n=ch.MG_c(n);
    sum_k_1=0;
    for k_1=0:b_n-1
        sum_m_1=0;
        for m=1:N
            a_m=ch.MG_a(m);
            b_m=ch.MG_b(m);
            c_m=ch.MG_c(m);
            sum_m_1=sum_m_1+a_m*((c_m)^-(b_m))*gamma(b_m)*gammainc(eta_0*(c_m),b_m);
        end
        sum_k_1=sum_k_1+((c_n*eta_0)^k_1/factorial(k_1))*sum_m_1;
    end
    temp_col1_1th_term=temp_col1_1th_term+a_n*c_n^(-b_n)*gamma(b_n)*exp(-c_n*eta_0)*sum_k_1;
end

% the second term
for n=1:N
    a_n=ch.MG_a(n);
    b_n=ch.MG_b(n);
    c_n=ch.MG_c(n);
    sum_k_1=0;
    for k_1=0:b_n-1
        sum_k_2=0;
        for k_2=0:k_1
            sum_m_1=0;
            for m=1:N
                a_m=ch.MG_a(m);
                b_m=ch.MG_b(m);
                c_m=ch.MG_c(m);
                sum_m_1=sum_m_1+a_m*(c_m+c_n*eta_0)^-(b_m+k_2)*gamma(b_m+k_2)*gammainc(eta_0*(c_n*eta_0+c_m),b_m+k_2);
            end
            sum_k_2=sum_k_2+nchoosek(k_1,k_2)*sum_m_1;
        end
        sum_k_1=sum_k_1+((c_n*eta_0)^k_1/factorial(k_1))*sum_k_2;
    end
    temp_col1_2th_term=temp_col1_2th_term+a_n*c_n^(-b_n)*gamma(b_n)*exp(-c_n*eta_0)*sum_k_1;
end

prob=temp_col1_1th_term-temp_col1_2th_term;

end

% P_{\mathsf{E}} or P_{\breve{\mathsf{E}}}
function prob = Col_0(ch)
N=length(ch.MG_a);
eta_0=ch.eta_0;
temp_col0_1th_term=0;
temp_col0_2th_term=0;

% the first term
for n=1:N
    a_n=ch.MG_a(n);
    b_n=ch.MG_b(n);
    c_n=ch.MG_c(n);
    sum_k_1=0;
    for k_1=0:b_n-1
        sum_m_1=0;
        for m=1:N
            a_m=ch.MG_a(m);
            b_m=ch.MG_b(m);
            c_m=ch.MG_c(m);
            sum_m_1=sum_m_1+a_m*(c_n+c_m)^(-(b_m+k_1))*gamma(b_m+k_1)*gammainc(eta_0*(c_n+c_m),b_m+k_1);
        end
        sum_k_1=sum_k_1+(c_n^k_1/factorial(k_1))*sum_m_1;
    end
    temp_col0_1th_term=temp_col0_1th_term+a_n*c_n^(-b_n)*gamma(b_n)*sum_k_1;
end

% the second term
for n=1:N
    a_n=ch.MG_a(n);
    b_n=ch.MG_b(n);
    c_n=ch.MG_c(n);
    sum_k_1=0;
    for k_1=0:b_n-1
        check_m=0;
        sum_m_1=0;
        for m=1:N
            a_m=ch.MG_a(m);
            b_m=ch.MG_b(m);
            c_m=ch.MG_c(m);
            sum_m_1=sum_m_1+a_m*(c_m)^-(b_m)*gamma(b_m)*gammainc(eta_0*(c_m),b_m);
            check_m=check_m+a_m*(c_m^(-b_m))*gamma(b_m);
        end
        if check_m~=1
        end
        sum_k_1=sum_k_1+((c_n*eta_0)^k_1/factorial(k_1))*sum_m_1;
    end
    temp_col0_2th_term=temp_col0_2th_term+a_n*c_n^(-b_n)*gamma(b_n)*exp(-c_n*eta_0)*sum_k_1;
end

prob=temp_col0_1th_term-temp_col0_2th_term;
end

% Pr{F}
function result = calculatePrF(ch)

N = length(ch.MG_a);
eta_0 = ch.eta_0;

result = 0;

for n = 1:N
    a_n = ch.MG_a(n);
    b_n = ch.MG_b(n);
    c_n = ch.MG_c(n);


    d_n = a_n * (c_n^(-b_n)) * factorial(b_n - 1);

    for m = 1:N
        a_m = ch.MG_a(m);
        b_m = ch.MG_b(m);
        c_m = ch.MG_c(m);

        term1 = exp(-c_n * eta_0) * sumGammaTermsF1(b_n, b_m, c_n, c_m, eta_0);
        term2 = sumGammaTermsF2(b_n, b_m, c_n, c_m, eta_0);


        result = result + a_m * d_n * (term1 - term2);
    end
end
end

function sum1 = sumGammaTermsF1(b_n, b_m, c_n, c_m, eta_0)
sum1 = 0;
for k1 = 0:b_n-1
    innerSum = 0;
    for k2 = 0:k1
        % Upper incomplete gamma function
        Gamma_upper = gamma(b_m + k2) * gammainc((c_m + c_n * eta_0) * eta_0 / (1 - eta_0), b_m + k2, 'upper');
        innerSum = innerSum + nchoosek(k1, k2)*(c_m + c_n * eta_0)^(-(b_m + k2)) * Gamma_upper;
    end
    sum1 = sum1 + ((c_n * eta_0)^k1 / factorial(k1)) * innerSum;
end
end

function sum2 = sumGammaTermsF2(b_n, b_m, c_n, c_m, eta_0)
sum2 = 0;
for k = 0:b_n-1
    % Upper incomplete gamma function
    Gamma_upper = gamma(b_m + k) * gammainc((c_m + c_n) * eta_0 / (1 - eta_0), b_m + k, 'upper');
    sum2 = sum2 + (c_n^k / factorial(k)) * (c_m + c_n)^(-(b_m + k)) * Gamma_upper;
end
end


% Pr{C^\prime}
function result = calculatePrC(ch)
% 引数
N = length(ch.MG_a);
eta_0 = ch.eta_0;

result = 0;

for n = 1:N
    a_n = ch.MG_a(n);
    b_n = ch.MG_b(n);
    c_n = ch.MG_c(n);

    d_n = a_n * (c_n^(-b_n)) * factorial(b_n - 1);

    for m = 1:N
        a_m = ch.MG_a(m);
        b_m = ch.MG_b(m);
        c_m = ch.MG_c(m);

        term1 = sumGammaTermsC1(b_n, b_m, c_n, c_m, eta_0);
        term2 = exp(-c_n * eta_0) * sumGammaTermsC2(b_n, b_m, c_n, c_m, eta_0);


        result = result + a_m * d_n * (term1 - term2);
    end
end
end

function sum1 = sumGammaTermsC1(b_n, b_m, c_n, c_m, eta_0)
sum1 = 0;
for k = 0:b_n-1
    % Lower incomplete gamma function
    Gamma_lower1 = gamma(b_m + k) * gammainc((c_m + c_n) * eta_0 / (1 - eta_0), b_m + k);
    Gamma_lower2 = gamma(b_m + k) * gammainc((c_m + c_n) * eta_0, b_m + k);

    sum1 = sum1 + (c_n^k / factorial(k)) * (c_m + c_n)^(-(b_m + k)) * (Gamma_lower1 - Gamma_lower2);
end
end

function sum2 = sumGammaTermsC2(b_n, b_m, c_n, c_m, eta_0)
sum2 = 0;
for k1 = 0:b_n-1
    innerSum = 0;
    for k2 = 0:k1
        % Lower incomplete gamma function
        Gamma_lower1 = gamma(b_m + k2) * gammainc((c_m + c_n * eta_0) * eta_0 / (1 - eta_0), b_m + k2);
        Gamma_lower2 = gamma(b_m + k2) * gammainc((c_m + c_n * eta_0) * eta_0, b_m + k2);

        innerSum = innerSum + nchoosek(k1, k2) * (c_m + c_n * eta_0)^(-(b_m + k2)) * (Gamma_lower1 - Gamma_lower2);
    end
    sum2 = sum2 + ((c_n * eta_0)^k1 / factorial(k1)) * innerSum;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% P\{E_1^1\}
function prob = Dec11(p,ch,L)
temp_11=0;
N=length(ch.MG_a);
eta_0=ch.eta_0;
for n=1:N
    a_n=ch.MG_a(n);
    b_n=ch.MG_b(n);
    c_n=ch.MG_c(n);
    temp_11=temp_11+a_n*(c_n^(-b_n))*gamma(b_n)*gammainc(c_n*eta_0,b_n,"upper");
end

term = L * log1p(-temp_11);
prob=2*p*(1-p)*( 1 - exp(term));
end

% P\{E_1^2\}
function prob = Dec21(p, ch, L)
A = P_1(ch) + Col_1(ch) + 2 * Col_0(ch);
B = Col_1(ch) + 2 * Col_0(ch);

X = L * log(A);
Y = L * log(B);
M = max(X, Y);

diff_exp = exp(X - M) - exp(Y - M);
diff_exp = max(diff_exp, 0);

prob = 2 * (p^2) * exp(M) * diff_exp;

end

% P\{E_2^2\}
function prob = Dec22(p, ch, L)
A = 2 * P_and(ch) + Col_1(ch) + 2 * P_2(ch) + 2 * P_1(ch);
B = 2 * P_and(ch) + Col_1(ch) + 2 * P_2(ch) + P_1(ch);
C = 2 * P_and(ch) + 2 * P_1(ch);

X = L * log1p(-A);
Y = L * log1p(-B);
Z = L * log1p(-C);

term_X = exp(X);
term_Y = exp(Y);
term_Z = exp(Z);

prob = p^2 * ( 1 + 2*term_X - 2*term_Y - term_Z );

end

function prob = Rec_pot(p,ch,L)
prob = (Dec11(p, ch, L) + Dec21(p, ch, L))/2;
end

function prob = Rec_nop(p,ch,L)
prob = (Dec11(p, ch, L) + Dec21(p, ch, L))/2 + Dec22(p, ch, L);
end