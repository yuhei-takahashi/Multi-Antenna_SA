format long
SNR=20;
R_set=5.8202;
L_set =[2000]; % the number of Antenna

% 1:rician(K=3), 2:Rayleigh, 3:Nakagami-m(m=5)
type=1;

disp("数値計算：混合ガンマ")
for Antnum=L_set
    for snr_set=SNR
        for R=R_set
            tic
            channel.snr = 10^(snr_set/10);

            if type==1
                N=20;
                rice_K=3;
                channel.MG_b=1:N;
                channel.MG_c=((1+rice_K)/(channel.snr)).*ones(1,N);
                channel.MG_a=phi(N,channel.MG_c,channel.MG_b,rice_K,channel.snr);
            elseif type==2
                % レイリー
                N=1;
                channel.MG_b=ones(1,N);
                channel.MG_c=(1/channel.snr).*ones(1,N);
                channel.MG_a=(1/channel.snr).*ones(1,N);
            elseif type==3
                % Nakagami-m
                N=1;
                m=5;
                channel.MG_b=m;
                channel.MG_c=(m/(channel.snr)).*ones(1,N);
                channel.MG_a=(m^m)/(gamma(m)*(channel.snr)^m);
            end

            % ファイル作成
            outputFolder = 'Ana_data_new_SNR';
            if ~exist(outputFolder, 'dir')
                mkdir(outputFolder);
            end
            file_name_Rs=sprintf('Ana_Rs_K=%d_R=%.4f_SNR_%d_MG_rice=%d.txt',Antnum,R,snr_set,rice_K);
            file_path_Rs=fullfile(outputFolder, file_name_Rs);
            filename_Rs=fopen(file_path_Rs,'w');
            fprintf(filename_Rs,'p R_s\n');

            tic
            
            fprintf('Antenna num=%d\n', Antnum);
            disp('p R_s');
            for p =0:0.01:1
                x = 1.0;  % Setting x=1 allows obtaining the conventional state-transition matrix of the Markov process

                matrix = K_state_TR_matrix(channel, p, x, R, Antnum);

                % 3x3 行列を定常分布計算のための遷移行列とみなし、その固有値と固有ベクトルを計算
                [V, D] = eig(matrix');

                % 固有値が 1 に最も近い固有ベクトルを取得
                [~, idx] = min(abs(diag(D) - 1.0));
                steady_state = V(:, idx);

                % 定常分布を正規化
                steady_state = steady_state / sum(steady_state);

                % Throughputを計算
                Throughput = K_calculateThroughput(channel,steady_state, p, R, Antnum);


                % 結果の表示
                R_s = R*Throughput;
                disp([num2str(p), ' ',  num2str(R_s, '%.6f')]);


                % txtファイル出力
                fprintf(filename_Rs,'%f %f\n',p,R_s);
            end
            toc

            % txtファイル出力
            fclose(filename_Rs);


        end
    end
end
