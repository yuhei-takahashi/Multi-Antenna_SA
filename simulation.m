clear;
SNR=25;
power=(10^(SNR/10));
% R=6.129; % 最適化結果 K=1
R=6.8987; % 最適化結果 K=2

eta_0=2^R-1;
K=2; % Number of anntenas
numUsers=2; % This code only supports m=2.
sysnum=1000000; % Number of iterations

user(numUsers)=struct('id',[],'active_flag',[],'decode_flag',[],'snr',[],'sinr',[],'buffer_flag',[],'current_flag',[]);

packet=zeros(1,20);
packet_test=zeros(1,20);

% fprintf('%s\n', 'p Throughput R_s');
fprintf('%s\n', 'p Throughput');

for p=0.0:0.05:1
    correct_packet=0;

    for i = 1:numUsers
        user(i).id = i;
        user(i).active_flag = 0; % 0 => active,1 => inactive
        user(i).decode_flag = 1; % 0 => non-decode,1 => completely decoded
        user(i).snr = zeros(1, K);
        user(i).sinr = zeros(1, K);
        user(i).buffer_flag = 0;
        user(i).current_flag = 0; % 現在のスロットでPDとなり，もう一つのアンテナで復号できてもPDをカウントしない
    end

    for iter = 1:sysnum

        for i = 1:numUsers
            random=rand(1);
            user(i).current_flag=0;
            if random <= p
                user(i).active_flag=1;
                user(i).decode_flag=0;
                user(i).snr=exprnd(power,1,K); % Rayleigh fading
                % user(i).snr=gamrnd(1,power,1,K);
            else
                user(i).active_flag=0;
                user(i).snr=zeros(1, K);
                user(i).sinr=zeros(1, K);
            end
        end

        % decoding
        count=0;
        for c=1:1000 % アンテナ間で繰り返し復号
            count_buf=count;

            for i=1:numUsers % 復号できた現在のパケットはSNRに0を設定
                if user(i).decode_flag == 1 && user(i).active_flag == 1
                    for ant_k = 1:K
                        user(i).snr(ant_k)=0;
                    end
                end
            end

            for ant_k = 1:K
                for i = 1:numUsers
                    other_users_snr = 1;
                    for j = 1:numUsers
                        if j ~= i
                            other_users_snr=other_users_snr+user(j).snr(ant_k);
                        end
                    end
                    user(i).sinr(ant_k)=user(i).snr(ant_k) / other_users_snr;
                end
            end

            for ant_k = 1:K
                for i = 1:numUsers
                    if user(i).decode_flag == 0 && user(i).active_flag == 1 %まだ復号されていない　かつ　アクティブ
                        if user(i).sinr(ant_k) > eta_0                  %i番目のデバイスは復号可能
                            correct_packet=correct_packet+1;
                            count=count+1;
                            user(i).buffer_flag=0;
                            user(i).decode_flag=1;

                            for j = 1:numUsers
                                if j ~= i
                                    if user(j).snr(ant_k)>eta_0 || user(j).buffer_flag==1
                                        if user(j).current_flag==0 && user(j).decode_flag == 0
                                            correct_packet=correct_packet+1;
                                            count=count+1;
                                            user(j).buffer_flag=0;
                                            user(j).decode_flag=1;
                                        end
                                    end
                                end
                            end
                        else
                            if user(i).snr(ant_k)>eta_0 && user(i).decode_flag==0
                                user(i).buffer_flag=1;
                                user(i).current_flag=1;
                            end
                        end
                    end
                end
            end
            if count==count_buf
                break;
            end
        end
    end

    % fprintf('%.2f %.6f\n %.6f\n',p, correct_packet/sysnum, R*correct_packet/sysnum);
    fprintf('%.2f %.6f\n',p, correct_packet/sysnum);
end
