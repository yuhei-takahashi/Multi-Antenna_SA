clear;
disp("Rician (approx. MG) simulation")


SNR_set=20;
R_set=[5.215404];% Code rate
rice_K_set=[3];
slotnum=10^4;% Number of slots
sysnum=10^3;% Number of iterations

for SNR=SNR_set
    for R=R_set
        for K=[1] % # Number of anntenas
            for rice_K=rice_K_set
                tic
                power=(10^(SNR/10));
     
                eta_0=2^R-1;
                numUsers=2;
                
                user(numUsers)=struct('id',[],'active_flag',[],'decode_flag',[],'snr',[],'sinr',[],'buffer_flag',[],'current_flag',[]);
                
                packet=zeros(1,20);
                packet_test=zeros(1,20);
                
                % fprintf('%s\n', 'p Throughput R_s');
                %fprintf('%s\n', 'p Throughput\n');
                
                %混合ガンマパラメータ設定
                %rice_K=3;
                %ライス
                %a=sqrt(2*rice_K*power);
                %A=a^2+2*power;
                %A=a^2/2+power;
                c_n=(1+rice_K)/(power); 
                %c_m=(1+rice_K)/(A);
                N=20;
                
                MG_b=1:N;
                %c_n=(1+K)/(A); 
                MG_c=c_n.*ones(1,N);
                %disp(a)
                MG_a=phi(N,MG_c,MG_b,rice_K,power);
        
                
                %ファイル設定用
                % 出力先のフォルダ名を指定
                outputFolder = 'Sim_data_SNR_new_rice';
                
                % フォルダが存在しない場合は作成
                if ~exist(outputFolder, 'dir')
                    mkdir(outputFolder);
                end
        
                %Ana:K=2_R=0.5_SNR=5_throughput.txt
                %スループット用
                file_name_tp=sprintf('Sim_tp_K=%d_R=%.4f_SNR_%d_MG_rice=%d.txt',K,R,SNR,rice_K);
                file_path_tp=fullfile(outputFolder, file_name_tp);
        
                filename_Throughput=fopen(file_path_tp,'w');
                fprintf(filename_Throughput,'p Throughput\n');
                fprintf('R:%.4f SNR:%d Antnum:%d rice:%d\n',R,SNR,K,rice_K);
        
                %R_s用
                file_name_Rs=sprintf('Sim_Rs_K=%d_R=%.4f_SNR_%d_MG_rice=%d.txt',K,R,SNR,rice_K);
                file_path_Rs=fullfile(outputFolder, file_name_Rs);
        
                filename_Rs=fopen(file_path_Rs,'w');
                fprintf(filename_Rs,'p R_s\n');
                
               
                for p=0.588688:0.1:1
                    
                
                    for i = 1:numUsers
                        user(i).id = i;
                        user(i).active_flag = 0; % 0 => active,1 => inactive
                        user(i).decode_flag = 1; % 0 => non-decode,1 => completely decoded
                        user(i).snr = zeros(1, K);
                        user(i).sinr = zeros(1, K);
                        user(i).buffer_flag = 0;
                        user(i).current_flag = 0; % 現在のスロットでPDとなり，もう一つのアンテナで復号できてもPDをカウントしない
                    end
                    tp_temp=0;
                    for iter = 1:sysnum
                        correct_packet=0;
                    for slot_ind=1:slotnum
                        for i = 1:numUsers
                            random=rand(1);
                            user(i).current_flag=0;
                            if random <= p
                                user(i).active_flag=1;
                                user(i).decode_flag=0;
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %user(i).snr=exprnd(power,1,K); % Rayleigh fading
                                user(i).snr=arrayfun(@(~) MGrnd(MG_a,MG_b,MG_c), 1:K);%MGrnd(MG_a,MG_b,MG_c); % Mixture Gamma
                                % user(i).snr=rice_avg_snr_random(rice_K,SNR,K); %ライス分布
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            else
                                user(i).active_flag=0;
                                user(i).snr=zeros(1, K);
                                user(i).sinr=zeros(1, K);
                            end
                        end
                
                        % decoding
                        count=0;
                        % アンテナ間で繰り返し復号
                        % cの最大が1のときアンテナ間SICを行わない
                        for c=1:1000000000000
                            count_buf=count;
                
                            % 復号できた現在のパケットはSNRに0を設定
                            for i=1:numUsers
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
                                                    if user(j).snr(ant_k)>eta_0
                                                        if user(j).decode_flag == 0
                                                            correct_packet=correct_packet+1;
                                                            count=count+1;
                                                            user(j).buffer_flag=0;
                                                            user(j).decode_flag=1;
                                                        end
                                                    elseif user(j).buffer_flag==1
                                                        if user(j).current_flag==0 && user(j).decode_flag == 0
                                                            correct_packet=correct_packet+1; % Comment out this line for w/o inter-slot SIC
                                                            count=count+1; 
                                                            user(j).buffer_flag=0;
                                                            user(j).decode_flag=1;
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                            if count==count_buf
                                break;
                            end
                        end
                
                        % PDパケットか判断する
                        for ant_k = 1:K
                            for i = 1:numUsers
                                if user(i).snr(ant_k)>eta_0 && user(i).decode_flag==0
                                    if user(i).buffer_flag==0
                                        user(i).buffer_flag=1;
                                        user(i).current_flag=1;
                                    end
                                end
                            end
                        end
                
                    end
                    %correct_packet=correct_packet/slotnum;
                    tp_temp=tp_temp+correct_packet/slotnum;
                    end
                    tp_temp=tp_temp/iter;
                    % 結果の表示
                    % R_s
                    R_s = R*tp_temp;
                
                    % fprintf('%.2f %.6f\n %.6f\n',p, correct_packet/sysnum, R*correct_packet/sysnum);
                    fprintf('p:%.2f Throughput:%.6f Rs:%.6f\n',p, tp_temp,R_s);
                    fprintf(filename_Throughput,'%f %f\n',p,tp_temp);
                    fprintf(filename_Rs,'%f %f\n',p,R_s);
                end
                fclose(filename_Throughput);
                fclose(filename_Rs);
                toc
            end
        end
    end
end
