function R_s = objective_function(x, Antnum,SN,type)
p = x(1);
R = x(2);
ch.snr = SN;

% 混合ガンマパラメータ設定
% ライス
if type==1
    N=20;
    rice_K=3;
    ch.MG_b=1:N;
    ch.MG_c=((1+rice_K)/(ch.snr)).*ones(1,N);
    ch.MG_a=phi(N,ch.MG_c,ch.MG_b,rice_K,ch.snr);
elseif type==2
    % レイリー
    N=1;
    ch.MG_b=ones(1,N);
    ch.MG_c=(1/ch.snr).*ones(1,N);
    ch.MG_a=(1/ch.snr).*ones(1,N);
elseif type==3
    % Nakagami-m
    N=1;
    m=5;
    ch.MG_b=m;
    ch.MG_c=(m/(ch.snr)).*ones(1,N);
    ch.MG_a=(m^m)/(gamma(m)*(ch.snr)^m);
end


x_0 = 1.0;  % Setting x=1 allows obtaining the conventional state-transition matrix of the Markov process

matrix = K_state_TR_matrix(ch, p, x_0, R, Antnum);

% 3x3 行列を定常分布計算のための遷移行列とみなし、その固有値と固有ベクトルを計算
[V, D] = eig(matrix');

% 固有値が 1 に最も近い固有ベクトルを取得
[~, idx] = min(abs(diag(D) - 1));
steady_state = V(:, idx);

% 定常分布を正規化
steady_state = steady_state / sum(steady_state);

Throughput = K_calculateThroughput(ch, steady_state, p, R, Antnum);

R_s = Throughput * R;

% R_sを最大化するために符号を反転
R_s = -R_s;
end
