clc; clear all; close all;




R = 8;
ss = 1;     % snapshot 1 2 3 4 5




% Read atom data and store chain information
N = 200;



marker = '.';
line = 2;



ch_id_directory = "/projects/p31790/arman/KG_PGN/R"+ R + "_rho05_analysis/Ri/" + N + "";

data_directory = "/projects/p31790/arman/KG_PGN/PGN_R" + R + "_rho05_N" + N + "/KG/KG_R" + R + "_N" + N + "/equ_bond_swap/simulations/npt_simulation/trajectory/";
num_atoms = readmatrix(data_directory + "dump.lammpsdata", FileType="text",Range=[4 1 4 1]);
chain_id = readmatrix(ch_id_directory + "/R" + R + "_N" + N + "_chain_id.txt", FileType="text");



data = readmatrix(data_directory + "dump.lammpsdata", FileType="text");
    
dim0 = readmatrix(data_directory + "dump.lammpsdata", FileType="text",Range=[6 1 8 2]);
dim1 = readmatrix(data_directory + "dump.lammpsdata", FileType="text",Range=[ss*(9+num_atoms)+6 1 9+ ss*(9+num_atoms)-1 2]);

data0 = data(1:num_atoms, :);
data1 = data(ss*(9+num_atoms)+1:ss*(9+num_atoms)+num_atoms, :);


% Sort data_traj_2 by the first column in ascending order
data0 = sortrows(data0, 1);
data1 = sortrows(data1, 1);





dim0(:,3) = dim0(:,2) - dim0(:,1);
dim1(:,3) = dim1(:,2) - dim1(:,1);


% Unwrapping coordinate, finding which bead has image flag
data0(:,4) = (data0(:,4) + data0(:,7))*dim0(1,3);
data0(:,5) = (data0(:,5) + data0(:,8))*dim0(2,3);
data0(:,6) = (data0(:,6) + data0(:,9))*dim0(3,3);

data1(:,4) = (data1(:,4) + data1(:,7))*dim1(1,3);
data1(:,5) = (data1(:,5) + data1(:,8))*dim1(2,3);
data1(:,6) = (data1(:,6) + data1(:,9))*dim1(3,3);

disp(:,1:3) = data0(:,4:6) - data1(:,4:6); 
msd(:,1) = (disp(:,1).^2+disp(:,2).^2+disp(:,3).^2);
    

%collecting au dwf
dwf = zeros(size(chain_id,2),1);
dwf_matrix = zeros(size(chain_id, 1), size(chain_id, 2));

for j = 1:size(chain_id,1)
    for k = 1:size(chain_id,2)

        dwf(k,1) = dwf(k,1) + msd(chain_id(j,k),1);
        dwf_matrix(j, k) = msd(chain_id(j, k), 1);
        
    end 
end 

std_dwf = std(dwf_matrix,1);

dwf_val = dwf/(size(chain_id,1));
 




%errorbar(1:N,dwf,std_dwf)


% Step 3: Plot the averaged data

plot(1:N, dwf_val, 'Color', 'r', 'Marker', marker, 'LineWidth', line,'DisplayName', num2str(N)); % Example plot style; customize as needed

hold on 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars -except R marker line ss;
% Read atom data and store chain information
N = 150;




ch_id_directory = "/projects/p31790/arman/KG_PGN/R"+ R + "_rho05_analysis/Ri/" + N + "";


data_directory = "/projects/p31790/arman/KG_PGN/PGN_R" + R + "_rho05_N" + N + "/KG/KG_R" + R + "_N" + N + "/equ_bond_swap/simulations/npt_simulation/trajectory/";
num_atoms = readmatrix(data_directory + "dump.lammpsdata", FileType="text",Range=[4 1 4 1]);
chain_id = readmatrix(ch_id_directory + "/R" + R + "_N" + N + "_chain_id.txt", FileType="text");



data = readmatrix(data_directory + "dump.lammpsdata", FileType="text");
    
dim0 = readmatrix(data_directory + "dump.lammpsdata", FileType="text",Range=[6 1 8 2]);
dim1 = readmatrix(data_directory + "dump.lammpsdata", FileType="text",Range=[ss*(9+num_atoms)+6 1 9+ ss*(9+num_atoms)-1 2]);

data0 = data(1:num_atoms, :);
data1 = data(ss*(9+num_atoms)+1:ss*(9+num_atoms)+num_atoms, :);


% Sort data_traj_2 by the first column in ascending order
data0 = sortrows(data0, 1);
data1 = sortrows(data1, 1);





dim0(:,3) = dim0(:,2) - dim0(:,1);
dim1(:,3) = dim1(:,2) - dim1(:,1);


% Unwrapping coordinate, finding which bead has image flag
data0(:,4) = (data0(:,4) + data0(:,7))*dim0(1,3);
data0(:,5) = (data0(:,5) + data0(:,8))*dim0(2,3);
data0(:,6) = (data0(:,6) + data0(:,9))*dim0(3,3);

data1(:,4) = (data1(:,4) + data1(:,7))*dim1(1,3);
data1(:,5) = (data1(:,5) + data1(:,8))*dim1(2,3);
data1(:,6) = (data1(:,6) + data1(:,9))*dim1(3,3);

disp(:,1:3) = data0(:,4:6) - data1(:,4:6); 
msd(:,1) = (disp(:,1).^2+disp(:,2).^2+disp(:,3).^2);
    

%collecting au dwf
dwf = zeros(size(chain_id,2),1);
dwf_matrix = zeros(size(chain_id, 1), size(chain_id, 2));

for j = 1:size(chain_id,1)
    for k = 1:size(chain_id,2)

        dwf(k,1) = dwf(k,1) + msd(chain_id(j,k),1);
        dwf_matrix(j, k) = msd(chain_id(j, k), 1);
        
    end 
end 

std_dwf = std(dwf_matrix,1);

dwf_val = dwf/(size(chain_id,1));
 




%errorbar(1:N,dwf,std_dwf)


% Step 3: Plot the averaged data
plot(1:N, dwf_val,  'Marker', marker, 'Color', [0.8500, 0.3250, 0.0980],  'LineWidth', line,'DisplayName', num2str(N)); % Example plot style; customize as needed

hold on 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars -except R marker line ss;
% Read atom data and store chain information
N = 100;



ch_id_directory = "/projects/p31790/arman/KG_PGN/R"+ R + "_rho05_analysis/Ri/" + N + "";


data_directory = "/projects/p31790/arman/KG_PGN/PGN_R" + R + "_rho05_N" + N + "/KG/KG_R" + R + "_N" + N + "/equ_bond_swap/simulations/npt_simulation/trajectory/";
num_atoms = readmatrix(data_directory + "dump.lammpsdata", FileType="text",Range=[4 1 4 1]);
chain_id = readmatrix(ch_id_directory + "/R" + R + "_N" + N + "_chain_id.txt", FileType="text");



data = readmatrix(data_directory + "dump.lammpsdata", FileType="text");
    
dim0 = readmatrix(data_directory + "dump.lammpsdata", FileType="text",Range=[6 1 8 2]);
dim1 = readmatrix(data_directory + "dump.lammpsdata", FileType="text",Range=[ss*(9+num_atoms)+6 1 9+ ss*(9+num_atoms)-1 2]);

data0 = data(1:num_atoms, :);
data1 = data(ss*(9+num_atoms)+1:ss*(9+num_atoms)+num_atoms, :);


% Sort data_traj_2 by the first column in ascending order
data0 = sortrows(data0, 1);
data1 = sortrows(data1, 1);





dim0(:,3) = dim0(:,2) - dim0(:,1);
dim1(:,3) = dim1(:,2) - dim1(:,1);


% Unwrapping coordinate, finding which bead has image flag
data0(:,4) = (data0(:,4) + data0(:,7))*dim0(1,3);
data0(:,5) = (data0(:,5) + data0(:,8))*dim0(2,3);
data0(:,6) = (data0(:,6) + data0(:,9))*dim0(3,3);

data1(:,4) = (data1(:,4) + data1(:,7))*dim1(1,3);
data1(:,5) = (data1(:,5) + data1(:,8))*dim1(2,3);
data1(:,6) = (data1(:,6) + data1(:,9))*dim1(3,3);

disp(:,1:3) = data0(:,4:6) - data1(:,4:6); 
msd(:,1) = (disp(:,1).^2+disp(:,2).^2+disp(:,3).^2);
    

%collecting au dwf
dwf = zeros(size(chain_id,2),1);
dwf_matrix = zeros(size(chain_id, 1), size(chain_id, 2));

for j = 1:size(chain_id,1)
    for k = 1:size(chain_id,2)

        dwf(k,1) = dwf(k,1) + msd(chain_id(j,k),1);
        dwf_matrix(j, k) = msd(chain_id(j, k), 1);
        
    end 
end 

std_dwf = std(dwf_matrix,1);

dwf_val = dwf/(size(chain_id,1));
 




%errorbar(1:N,dwf,std_dwf)


% Step 3: Plot the averaged data
plot(1:N, dwf_val, 'Marker', marker, 'Color', [0.4940 0.1840 0.5560], 'LineWidth', line,'DisplayName', num2str(N)); % Example plot style; customize as needed

hold on 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars -except R marker line ss;
% Read atom data and store chain information
N = 50;



ch_id_directory = "/projects/p31790/arman/KG_PGN/R"+ R + "_rho05_analysis/Ri/" + N + "";


data_directory = "/projects/p31790/arman/KG_PGN/PGN_R" + R + "_rho05_N" + N + "/KG/KG_R" + R + "_N" + N + "/equ_bond_swap/simulations/npt_simulation/trajectory/";
num_atoms = readmatrix(data_directory + "dump.lammpsdata", FileType="text",Range=[4 1 4 1]);
chain_id = readmatrix(ch_id_directory + "/R" + R + "_N" + N + "_chain_id.txt", FileType="text");



data = readmatrix(data_directory + "dump.lammpsdata", FileType="text");
    
dim0 = readmatrix(data_directory + "dump.lammpsdata", FileType="text",Range=[6 1 8 2]);
dim1 = readmatrix(data_directory + "dump.lammpsdata", FileType="text",Range=[ss*(9+num_atoms)+6 1 9+ ss*(9+num_atoms)-1 2]);

data0 = data(1:num_atoms, :);
data1 = data(ss*(9+num_atoms)+1:ss*(9+num_atoms)+num_atoms, :);


% Sort data_traj_2 by the first column in ascending order
data0 = sortrows(data0, 1);
data1 = sortrows(data1, 1);





dim0(:,3) = dim0(:,2) - dim0(:,1);
dim1(:,3) = dim1(:,2) - dim1(:,1);


% Unwrapping coordinate, finding which bead has image flag
data0(:,4) = (data0(:,4) + data0(:,7))*dim0(1,3);
data0(:,5) = (data0(:,5) + data0(:,8))*dim0(2,3);
data0(:,6) = (data0(:,6) + data0(:,9))*dim0(3,3);

data1(:,4) = (data1(:,4) + data1(:,7))*dim1(1,3);
data1(:,5) = (data1(:,5) + data1(:,8))*dim1(2,3);
data1(:,6) = (data1(:,6) + data1(:,9))*dim1(3,3);

disp(:,1:3) = data0(:,4:6) - data1(:,4:6); 
msd(:,1) = (disp(:,1).^2+disp(:,2).^2+disp(:,3).^2);
    

%collecting au dwf
dwf = zeros(size(chain_id,2),1);
dwf_matrix = zeros(size(chain_id, 1), size(chain_id, 2));

for j = 1:size(chain_id,1)
    for k = 1:size(chain_id,2)

        dwf(k,1) = dwf(k,1) + msd(chain_id(j,k),1);
        dwf_matrix(j, k) = msd(chain_id(j, k), 1);
        
    end 
end 

std_dwf = std(dwf_matrix,1);

dwf_val = dwf/(size(chain_id,1));
 




%errorbar(1:N,dwf,std_dwf)


% Step 3: Plot the averaged data
plot(1:N, dwf_val, 'Marker', marker,'Color', [0.4660 0.6740 0.1880],  'LineWidth', line,'DisplayName', num2str(N)); % Example plot style; customize as needed

hold on 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clearvars -except R marker line ss;
% Read atom data and store chain information
N = 20;



ch_id_directory = "/projects/p31790/arman/KG_PGN/R"+ R + "_rho05_analysis/Ri/" + N + "";


data_directory = "/projects/p31790/arman/KG_PGN/PGN_R" + R + "_rho05_N" + N + "/KG/KG_R" + R + "_N" + N + "/equ_bond_swap/simulations/npt_simulation/trajectory/";
num_atoms = readmatrix(data_directory + "dump.lammpsdata", FileType="text",Range=[4 1 4 1]);
chain_id = readmatrix(ch_id_directory + "/R" + R + "_N" + N + "_chain_id.txt", FileType="text");



data = readmatrix(data_directory + "dump.lammpsdata", FileType="text");
    
dim0 = readmatrix(data_directory + "dump.lammpsdata", FileType="text",Range=[6 1 8 2]);
dim1 = readmatrix(data_directory + "dump.lammpsdata", FileType="text",Range=[ss*(9+num_atoms)+6 1 9+ ss*(9+num_atoms)-1 2]);

data0 = data(1:num_atoms, :);
data1 = data(ss*(9+num_atoms)+1:ss*(9+num_atoms)+num_atoms, :);


% Sort data_traj_2 by the first column in ascending order
data0 = sortrows(data0, 1);
data1 = sortrows(data1, 1);





dim0(:,3) = dim0(:,2) - dim0(:,1);
dim1(:,3) = dim1(:,2) - dim1(:,1);


% Unwrapping coordinate, finding which bead has image flag
data0(:,4) = (data0(:,4) + data0(:,7))*dim0(1,3);
data0(:,5) = (data0(:,5) + data0(:,8))*dim0(2,3);
data0(:,6) = (data0(:,6) + data0(:,9))*dim0(3,3);

data1(:,4) = (data1(:,4) + data1(:,7))*dim1(1,3);
data1(:,5) = (data1(:,5) + data1(:,8))*dim1(2,3);
data1(:,6) = (data1(:,6) + data1(:,9))*dim1(3,3);

disp(:,1:3) = data0(:,4:6) - data1(:,4:6); 
msd(:,1) = (disp(:,1).^2+disp(:,2).^2+disp(:,3).^2);
    

%collecting au dwf
dwf = zeros(size(chain_id,2),1);
dwf_matrix = zeros(size(chain_id, 1), size(chain_id, 2));

for j = 1:size(chain_id,1)
    for k = 1:size(chain_id,2)

        dwf(k,1) = dwf(k,1) + msd(chain_id(j,k),1);
        dwf_matrix(j, k) = msd(chain_id(j, k), 1);
        
    end 
end 

std_dwf = std(dwf_matrix,1);

dwf_val = dwf/(size(chain_id,1));
 
plot(1:N, dwf_val, 'Marker', marker, 'Color', 'b', 'LineWidth', line,'DisplayName', num2str(N)); % Example plot style; customize as needed








xlabel('N_i', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('<u^2_i>', 'FontSize', 16, 'FontWeight', 'bold');
set(gca, 'TickDir', 'out');  % 'out' places the ticks on the opposite side
set(gca, 'FontSize', 16);  % Set the font size for the ticks
set(gca, 'FontWeight', 'bold');  % Set the font weight for the ticks
set(gca, 'LineWidth', 2);  % Set the line width for the plot box
ylim([0.02,0.04]);
xlim([-1,201])
legend('Location', 'southeast', 'FontSize', 16, 'FontWeight', 'bold');
titleStr = ['NP Radius (σ): ' num2str(R)];
bbox = annotation('rectangle', [0.58, 0.875, 0.325, 0.05], 'Color', 'black', 'FaceColor', 'white', 'FaceAlpha', 0, 'LineWidth', 2);
titleText = text('Units', 'normalized', 'Position', [0.6, 0.97], 'String', titleStr, 'FontSize', 16, 'FontWeight', 'bold');

output_figure_file = "/projects/p31790/arman/KG_PGN/R_all_analysis/u2_v_n/R" + R + "_ss" + ss + "_averaged_dwf_plot.png";
saveas(gcf, output_figure_file);













%     module load matlab/r2022a

 %           matlab -nodisplay -nosplash -nodesktop -r "run('/projects/p31790/arman/KG_PGN/R_all_analysis/u2_v_n/combined_debye_waller.m')"



