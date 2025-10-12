SNR = 25; % dB
power = (10^(SNR / 10));
K = 3; % # of antennas
numUsers = 20; % # of users
sysnum = 1000; % Number of slots
eta_0_func = @(R) 2^R - 1; % decoding threshold
p_values = 0.095:0.005:0.105; % trans prob
R_values = 7.40466875:0.01:7.40466875; % coding rate
iteration = 10^2;


% Mixture Gamma parameter for Rician fading
global rice_K;
rice_K=3;
global N;
N=20;

% initialize a 3D array
Rs_matrix = zeros(length(R_values), length(p_values));

fprintf('# of Ant = %d, # of users = %d\n', K, numUsers);
fprintf('p       R       R_s\n');
for j = 1:length(p_values)
    before_Rs_matrix=0;
    for i = 1:length(R_values)
        R = R_values(i);
        p = p_values(j);
        total_Rs = 0;

        % average the results over multiple iterations
        for t = 1:iteration
            total_Rs = total_Rs + optimize_func([R, p], power, eta_0_func, K, numUsers, sysnum);
        end
        Rs_matrix(i, j) = total_Rs / iteration;
        fprintf('%.4f\t%.4f\t%.4f\n', p, R, Rs_matrix(i, j));

    end
end

% 3D surface plot
[P_mesh, R_mesh] = meshgrid(p_values, R_values);
surf(P_mesh, R_mesh, Rs_matrix)
xlabel('p');
ylabel('R');
zlabel('Rs');
title('Throughput Rs vs Transmission Probability p and Code Rate R');
% combined_matrix = [P_mesh(:), R_mesh(:), Rs_matrix(:)];
% disp(combined_matrix);

% objective function
function Rs = optimize_func(x, power, eta_0_func, K, numUsers, sysnum)

% Mixture Gamma parameter for Rician fading
global N;
global rice_K;
c_n=(1+rice_K)/(power);
MG_b=1:N;
MG_c=c_n.*ones(1,N);
MG_a=phi(N,MG_c,MG_b,rice_K,power);

R = x(1); % coding rate
p = x(2); % transmission probability
eta_0 = eta_0_func(R); % SINR threshold

correct_packet = 0;
user_decode_flag = zeros(1, numUsers);

% initialize the slot
slot(sysnum) = struct('user_snr', [], 'user_sinr', [], 'decode_flag', []);
for s = 1:sysnum
    for i = 1:numUsers
        slot(s).user_snr(i).ant = zeros(1, K);
        slot(s).user_sinr(i).ant = zeros(1, K);
    end
end

% set the SNR for each antenna of the user
for s = 1:sysnum
    for i = 1:numUsers
        random = rand(1);
        if random <= p
            %%%%% Rayleigh %%%%
             slot(s).user_snr(i).ant = exprnd(power, 1, K);
            %%%%% Mixture Gamma %%%%
            % slot(s).user_snr(i).ant = arrayfun(@(~) MGrnd(MG_a,MG_b,MG_c), 1:K);
        else
            slot(s).user_snr(i).ant = zeros(1, K);
        end
        slot(s).user_sinr(i).ant = zeros(1, K);
    end
end


slot_index = 0;
sic_times = 0;
for c = 0:10000000
    count_buf = sic_times;

    % calculate the SINR for each user
    for s = 1:sysnum
        for ant_num = 1:K
            for i = 1:numUsers
                other_users_snr = 1;
                for j = 1:numUsers
                    if j ~= i
                        other_users_snr = other_users_snr + slot(s).user_snr(j).ant(ant_num);
                    end
                end
                slot(s).user_sinr(i).ant(ant_num) = slot(s).user_snr(i).ant(ant_num) / other_users_snr;
            end
        end
    end

    % SIC-based decoding
    flag = false;
    current_slot = 0;
    for s = 1:sysnum
        current_slot = current_slot + 1;

        for ant_num = 1:K
            for i = 1:numUsers
                if slot(s).user_sinr(i).ant(ant_num) > eta_0
                    sic_times = sic_times + 1;
                    for j = 1:1:current_slot
                        slot(j).user_snr(i).ant = zeros(1, K);
                    end
                    if user_decode_flag(i) == 0
                        correct_packet = correct_packet + 1;
                        user_decode_flag(i) = 1;
                    end
                    flag = true;
                    break;
                end
            end
            if flag
                break;
            end
        end

        buffer_slot_index = slot_index;
        slot_index = max(current_slot, buffer_slot_index);

        if flag
            break;
        end
    end

    for i = 1:numUsers
        if user_decode_flag(i) == 1
            for s = 1:1:slot_index
                slot(s).user_snr(i).ant = zeros(1, K);
            end
        end
    end

    if buffer_slot_index <= slot_index
        user_decode_flag = zeros(1, numUsers);
    end

    if count_buf == sic_times
        break;
    end
end

% calculation of the sum rate
Rs = R * (correct_packet / sysnum);
end
