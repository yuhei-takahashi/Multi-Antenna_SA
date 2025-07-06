function R_s = objective_function(x, Antnum,SN,type)
p = x(1);
R = x(2);
ch.snr = SN;

% Rician fading
if type==1
    N=20;
    rice_K=3;
    ch.MG_b=1:N;
    ch.MG_c=((1+rice_K)/(ch.snr)).*ones(1,N);
    ch.MG_a=phi(N,ch.MG_c,ch.MG_b,rice_K,ch.snr);
% Rayleigh fading
elseif type==2
    N=1;
    ch.MG_b=ones(1,N);
    ch.MG_c=(1/ch.snr).*ones(1,N);
    ch.MG_a=(1/ch.snr).*ones(1,N);
% Nakagami-m    
elseif type==3
    N=1;
    m=5;
    ch.MG_b=m;
    ch.MG_c=(m/(ch.snr)).*ones(1,N);
    ch.MG_a=(m^m)/(gamma(m)*(ch.snr)^m);
end

x_0 = 1.0;  % Setting x=1 allows obtaining the conventional state-transition matrix of the Markov process

matrix = K_state_TR_matrix(ch, p, x_0, R, Antnum);

% Treat the 3x3 matrix as a transition matrix to compute the steady-state distribution
[V, D] = eig(matrix');

% Extract the eigenvector corresponding to the eigenvalue closest to 1
[~, idx] = min(abs(diag(D) - 1));
steady_state = V(:, idx);

% Normalize the steady-state distribution
steady_state = steady_state / sum(steady_state);

% Compute throughput based on steady-state distribution
Throughput = K_calculateThroughput(ch, steady_state, p, R, Antnum);

R_s = Throughput * R;

% Negate for minimization (since optimization routine minimizes by default)
R_s = -R_s;
end
