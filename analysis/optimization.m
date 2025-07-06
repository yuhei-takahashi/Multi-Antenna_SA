tic
L = [2]; % Number of antennas

SNR_dB=15;
snr = 10^(SNR_dB/10);

% 1:rician(K=3), 2:Rayleigh, 3:Nakagami-m(m=5)
type=1;

for Antnum=L

    % initial settings [p,R]
    x0 = [0.6, 4]; % Please adjust these settings to avoid local optima!!!!!

    % define the range [p,R]
    lb = [0.000001,0.000001];  % lower bound
    ub = [1, Inf];  % upper bound

    % fmincon options
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');

    % optimization
    [x_opt, fval_opt] = fmincon(@(x) objective_function(x, Antnum,snr,type), x0, [], [], [], [], lb, ub, @nonlcon, options);

    % Display of optimization results
    p_opt = x_opt(1);
    R_opt = x_opt(2);
    R_s_opt = -fval_opt;

    disp(['L: ', num2str(Antnum)]);
    disp(['Optimal p: ', num2str(p_opt, '%.4f')]);
    disp(['Optimal R: ', num2str(R_opt, '%.8f')]);
    disp(['Optimal R_s: ', num2str(R_s_opt, '%.8f')]);

    % output to a file
    save_optimization_results(Antnum, SNR_dB, p_opt, R_opt, R_s_opt);

    % generate a 3D plot
    generate_3D_plot_and_save(Antnum, SNR_dB, snr, type, p_opt, R_opt, R_s_opt);

    toc
end

function save_optimization_results(Antnum, SNR_dB, p_opt, R_opt, R_s_opt)
%   Save optimization results to text files in a specified folder

outputFolder = 'Ana_data_new_SNR';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

file_name_Rs = sprintf('Ana_Rs_L=%d_SNR%d.txt', Antnum, SNR_dB);
file_path_Rs = fullfile(outputFolder, file_name_Rs);
filename_Rs = fopen(file_path_Rs, 'w');

file_name_opt = sprintf('Opt_Ana_Rs_L=%d_SNR%d.txt', Antnum, SNR_dB);
file_path_opt = fullfile(outputFolder, file_name_opt);
filename_opt = fopen(file_path_opt, 'w');

fprintf(filename_opt, 'Optimal p: %f\n', p_opt);
fprintf(filename_opt, 'Optimal R: %f\n', R_opt);
fprintf(filename_opt, 'Optimal R_s: %f\n', R_s_opt);

fclose(filename_Rs);
fclose(filename_opt);
end

function generate_3D_plot_and_save(Antnum, SNR_dB, snr, type, p_opt, R_opt, R_s_opt)
%   Generate a 3D surface plot, save data and optimal point info to files.

% ---- Generate p-R grid ----
p_values = linspace(0, 1, 30);
R_values = linspace(0, 10, 30);
[P, R] = meshgrid(p_values, R_values);
R_s_values = arrayfun(@(p, R) -objective_function([p, R], Antnum, snr, type), P, R);

% ---- Generate 3D surface plot ----
figure;
surf(P, R, R_s_values, 'EdgeColor', 'none');
xlabel('$p$', 'Interpreter', 'latex');
ylabel('$R$', 'Interpreter', 'latex');
zlabel('$R_s$', 'Interpreter', 'latex');

hold on;
scatter3(p_opt, R_opt, R_s_opt, 50, 'black', 'filled');

z_offset = 0.1;
text_str = sprintf('$p$: %.2f\n$R$: %.2f\n$R_s$: %.2f', p_opt, R_opt, R_s_opt);
text(p_opt, R_opt, R_s_opt + z_offset, text_str, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10, 'Interpreter', 'latex');


% ---- Save data to file ----
outputFolder = 'Ana_data_new_SNR';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

file_name_Rs = sprintf('Ana_Rs_L=%d_SNR%d.txt', Antnum, SNR_dB);
file_path_Rs = fullfile(outputFolder, file_name_Rs);
filename_Rs = fopen(file_path_Rs, 'w');
if filename_Rs == -1
    error('Failed to open Rs file for writing.');
end
fprintf(filename_Rs, 'p R R_s\n');
fprintf(filename_Rs, '%f %f %f\n', [P(:), R(:), R_s_values(:)]');
fclose(filename_Rs);

file_name_opt = sprintf('Opt_Ana_Rs_L=%d_SNR%d.txt', Antnum, SNR_dB);
file_path_opt = fullfile(outputFolder, file_name_opt);
filename_opt = fopen(file_path_opt, 'w');
if filename_opt == -1
    error('Failed to open opt file for writing.');
end
fprintf(filename_opt, 'Optimal p: %f\n', p_opt);
fprintf(filename_opt, 'Optimal R: %f\n', R_opt);
fprintf(filename_opt, 'Optimal R_s: %f\n', R_s_opt);
fclose(filename_opt);

hold off;
end
