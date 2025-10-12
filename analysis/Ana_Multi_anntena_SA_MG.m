format long
SNR=20;
R_set=5;
L_set =[2]; % the number of Antenna

% 1:rician(K=3), 2:Rayleigh, 3:Nakagami-m(m=5)
type=1;

disp("< MG analysis>")
for Antnum=L_set
    for snr_set=SNR
        for R=R_set
            tic
            channel.snr = 10^(snr_set/10);

            if type==1
                % rician
                N=20;
                rice_K=3;
                channel.MG_b=1:N;
                channel.MG_c=((1+rice_K)/(channel.snr)).*ones(1,N);
                channel.MG_a=phi(N,channel.MG_c,channel.MG_b,rice_K,channel.snr);
            elseif type==2
                % Rayleigh
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

            % output to a file
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

                % Consider the 3x3 matrix as a transition matrix for computing the steady-state distribution
                [V, D] = eig(matrix');

                % Select the eigenvector corresponding to the eigenvalue closest to 1
                [~, idx] = min(abs(diag(D) - 1.0));
                steady_state = V(:, idx);

                % Normalize the steady-state distribution
                steady_state = steady_state / sum(steady_state);

                % Compute the throughput
                Throughput = K_calculateThroughput(channel,steady_state, p, R, Antnum);

                % Display the result
                R_s = R*Throughput;
                disp([num2str(p), ' ',  num2str(R_s, '%.6f')]);

                % Write to a text file
                fprintf(filename_Rs,'%f %f\n',p,R_s);
            end
            toc

            fclose(filename_Rs);
        end
    end
end
