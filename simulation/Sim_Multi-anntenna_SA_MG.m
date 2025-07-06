clear;
disp("<SIC-Based SA simulation>")

% ========== Parameters ==========
numAntennas=1;
numUsers=2;
SNR=20;
R=5.21540459;
rice_K=3;
slotnum=10^4;
sysnum=10^3;
% --- Show simulation parameters ---
fprintf('=== Simulation Parameters ===\n');
fprintf('# of Antennas : %d\n', numAntennas);
fprintf('# of users    : %d\n', 2);
fprintf('SNR           : %d dB\n', SNR);
fprintf('Code Rate R   : %.4f\n', R);
fprintf('Rician K      : %d\n', rice_K);
fprintf('Slot Number   : %d\n', slotnum);
fprintf('Iterations    : %d\n', sysnum);
fprintf('=============================\n\n');

tic
power=(10^(SNR/10));
eta_0=2^R-1;


% User initialization
user(numUsers)=struct('id',[],'active_flag',[],'decode_flag',[],'snr',[],'sinr',[],'buffer_flag',[],'current_flag',[]);

packet=zeros(1,20);
packet_test=zeros(1,20);

% MG distribution parameters
c_n=(1+rice_K)/(power);
N=20;
MG_b=1:N;
MG_c=c_n.*ones(1,N);
MG_a=phi(N,MG_c,MG_b,rice_K,power);


% ========== Output Setup ==========
outputFolder = 'Sim_data_SNR_new_rice';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

file_name_tp=sprintf('Sim_tp_K=%d_R=%.4f_SNR_%d_MG_rice=%d.txt',numAntennas,R,SNR,rice_K);
file_path_tp=fullfile(outputFolder, file_name_tp);
filename_Throughput=fopen(file_path_tp,'w');
fprintf(filename_Throughput,'p Throughput\n');

file_name_Rs=sprintf('Sim_Rs_K=%d_R=%.4f_SNR_%d_MG_rice=%d.txt',numAntennas,R,SNR,rice_K);
file_path_Rs=fullfile(outputFolder, file_name_Rs);
filename_Rs=fopen(file_path_Rs,'w');
fprintf(filename_Rs,'p R_s\n');

% ========== Simulation Loop ==========
for p=0.1:0.1:1
    for i = 1:numUsers
        user(i).id = i;
        user(i).active_flag = 0;     % 0 => active,1 => inactive
        user(i).decode_flag = 1;     % 0 => non-decode,1 => completely decoded
        user(i).snr = zeros(1, numAntennas);
        user(i).sinr = zeros(1, numAntennas);
        user(i).buffer_flag = 0;     % 0 => no packet in buffer,1 => packet(s) stored in buffer
        user(i).current_flag = 0;    % If a packet becomes PD in the current slot, but is decoded at other antennas, do not count it as PD
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
                    % ========== fading generation ==========
                    %user(i).snr=exprnd(power,1,numAntennas); % Rayleigh fading
                    user(i).snr=arrayfun(@(~) MGrnd(MG_a,MG_b,MG_c), 1:numAntennas);% Mixture Gamma
                    % user(i).snr=rice_avg_snr_random(rice_K,SNR,numAntennas); % Original Rician
                    % =======================================
                else
                    user(i).active_flag=0;
                    user(i).snr=zeros(1, numAntennas);
                    user(i).sinr=zeros(1, numAntennas);
                end
            end

            % ========== Decoding Processing ==========
            count=0;

            % If MAX_ANTENNA_SIC_ITER is set to 1, antenna inter-SIC is disabled
            MAX_ANTENNA_SIC_ITER=1e12;
            for rep = 1:MAX_ANTENNA_SIC_ITER
                count_buf=count;

                % Zero SNR for successfully decoded packets
                for i=1:numUsers
                    if user(i).decode_flag == 1 && user(i).active_flag == 1
                        for ant_k = 1:numAntennas
                            user(i).snr(ant_k)=0;
                        end
                    end
                end

                % SINR computation
                for ant_k = 1:numAntennas
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

                % Decoding check
                for ant_k = 1:numAntennas
                    for i = 1:numUsers
                        if user(i).decode_flag == 0 && user(i).active_flag == 1 % Undecoded and active
                            if user(i).sinr(ant_k) > eta_0                  % The device is decodable
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

            % --- PD packet check ---
            for ant_k = 1:numAntennas
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

        tp_temp=tp_temp+correct_packet/slotnum;
    end
    tp_temp=tp_temp/iter;
    R_s = R*tp_temp;

    fprintf('p:%.2f Throughput:%.6f Rs:%.6f\n',p, tp_temp,R_s);
    fprintf(filename_Throughput,'%f %f\n',p,tp_temp);
    fprintf(filename_Rs,'%f %f\n',p,R_s);
end
fclose(filename_Throughput);
fclose(filename_Rs);
toc
