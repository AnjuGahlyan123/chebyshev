tic
clc
clear;
n = 28;                  % Number of Chebyshev Functions taken
k = 0;                   % a^2 = k^2 + m^2 ; k and m being the x and y wavenumber
Rvv = 100;      % Vertical Rayleigh Number
syms a

%% Values of Parameters

rho = 1180;               % Density
g = 9.8;                  % Gravity 
alpha = 5.2e-4;           % Temperature volume expansion coefficient
d = 0.001;                % Thickness
Kt = 0.59;                % Thermal conductivity
Cvh = 3545.76;            % Heat Capacity
mhu = 0.007;              % Dynamic viscosity
kappa = Kt / (rho * Cvh); % Thermal Diffusivity
Kv = mhu / rho;           % Kinematic viscosity

%% Define Rh (Horizontal Rayleigh Number) value 

Rh = @(beta_initial) (rho * alpha * g * beta_initial * (d^3) * Kv) / ((kappa^2) * mhu);

%% Eigenvalue Calculation Function

function eigval1 = eigenvalue(beta, rho, g, alpha, d, mhu, kappa, Kv, Rv, a, n, k)
    Rh = @(beta_initial) (rho * alpha * g * beta_initial * (d^3) * Kv) / ((kappa^2) * mhu);    
    Rh_val = Rh(beta);  
    [A, B] = define_matrix(n, k, Rh_val, Rv, a);
    eigenvalues = eig(A, B);
    eigenvalues = eigenvalues(isfinite(eigenvalues) );
    [~, sortIdx] = sort(real(eigenvalues), 'descend');  
    sorted_eigenvalues = eigenvalues(sortIdx);
    eigval1 = real(sorted_eigenvalues(1));
end


%% Finding Critical Beta using fzero

function beta = betacri(n, rho, g, alpha, d, mhu, kappa, Kv, Rv, a, k)
    betaa = logspace(-20, 20, 120);
    Eig = zeros(1, length(betaa));  
    x0 = [1.0e-20, 1.0e20];  
    for idx = 1:length(betaa)
        beta1 = betaa(idx);
        eigval1 = eigenvalue(beta1, rho, g, alpha, d, mhu, kappa, Kv, Rv, a, n, k);  
        Eig(idx) = eigval1;
        if (idx >= 2 && Eig(idx - 1) * Eig(idx) < 0)
            x0(1) = betaa(idx - 1);
            x0(2) = betaa(idx);
        end
    end
    beta = fzero(@(beta) eigenvalue(beta, rho, g, alpha, d, mhu, kappa, Kv, Rv, a, n, k), x0);
end

%% Optimization for a

format longG
Critical_values=zeros(length(Rvv),5);
for idx=1:length(Rvv)
    Rv=Rvv(idx);
    a_lower = 0.1;  
    a_upper = 30;  
    [a_critical, min_sigma] = fminbnd(@(a) betacri(n, rho, g, alpha, d, mhu, kappa, Kv, Rv, a, k), a_lower, a_upper);         
    critical_beta = betacri(n, rho, g, alpha, d, mhu, kappa, Kv, Rv, a_critical, k);
    critical_Rh = Rh(critical_beta);
    [Eig, EigVec] =eigenvaluee(n,k,a_critical,critical_Rh,Rv)
    result=[Rv,critical_Rh,a_critical];
    Critical_values(idx,:)=[result,real(Eig),imag(Eig)];
end

colNames = {'Rv_c','Rh_c', 'a_c', 'Real_part_Eig','Ima_part_Eig'};
criticalTable = array2table(Critical_values, 'VariableNames', colNames);
filename = 'critical_matrix1.xlsx'; 
writetable(criticalTable, filename);
loadedTable = readtable(filename);
disp(loadedTable)

% figure
% plot(Critical_values(:,1),Critical_values(:,2))
% title("Rh and Rv plot for k=0")
% xlabel("Rv")
% ylabel("Rh")

toc




function [A,B]=define_matrix(n,k,Rh,Rv,a)
D = zeros(n+2, n+2);
for i = 1:n+2
    for j = 1:n+2
        if i == 1
            D(i, 2*j) = 2*j - 1;
        else
            D(i, i+2*j-1) = 4*j + 2*i - 4;
        end
    end
end
D = D(1:n+2, 1:n+2);
D2 = D^2;
P = Pz(n);
M = Mz(n);

% Boundary conditions
if mod(n+2,2) == 1
    for i = 1:2:n
        D(:,i) = D(:,i) - D(:,n+2);
        D2(:,i) = D2(:,i) - D2(:,n+2);
    end
    for i = 2:2:n
        D(:,i) = D(:,i) - D(:,n+1);
        D2(:,i) = D2(:,i) - D2(:,n+1);
    end
else
    for i = 1:2:n
        D(:,i) = D(:,i) - D(:,n+1);
        D2(:,i) = D2(:,i) - D2(:,n+1);
    end
    for i = 2:2:n
        D(:,i) = D(:,i) - D(:,n+2);
        D2(:,i) = D2(:,i) - D2(:,n+2);
    end
end
D = D(1:n, 1:n);
D2 = D2(1:n, 1:n);

% Construct matrices
I = eye(n);
O = zeros(n);
A11 = 4*D2 - a^2 * I; 
A12 = a^2 * I; 
A21 = (1i * 2 * k * Rh / a^2) * D + (Rv - ((Rh^2) / 24)) * I + ((Rh^2) / 8) * P;
A22 = 4*D2 - a^2 * I - (1i * k * Rh / 2) * M;
A = [A11, A12; A21, A22];
B = [O, O; O, I];    
end

%% Helper function for P

function P = Pz(N)
P = zeros(N, N);
P(1, 1) = 1/2;
P(2, 2) = 3/4;
P(3, 1) = 1/2;
for i = 3:N
    P(i, i) = 1/2;
end
for i = 1:N-2
    P(i, i+2) = 1/4;
end
for i = 2:N-2
    P(i+2, i) = 1/4;
end
end

%% Helper Function for M

function M = Mz(N)
M = zeros(N, N);
for i = 1:N-1
    M(i, i+1) = 1/2;
end
for i = 3:N
    M(i, i-1) = 1/2;
end
M(2, 1) = 1;
end


function [Eig,eigenvectors ] = eigenvaluee(n, k, a, Rh, Rv)

    [A, B] = define_matrix(n, k, Rh, Rv, a);

    %% Solve the generalized eigenvalue problem by QZ Algorithm

    [S, T, Q, Z] = qz(A, B); % QZ decomposition
    alpha = diag(S); 
    beta = diag(T);
    values = alpha ./ beta;
    id = isfinite(values);
    id = logical(id == 1);
    eigenvalues = values(id);
    [~, sortIdx] = sort(real(eigenvalues), 'descend'); 
    sorted_eigenvalues = eigenvalues(sortIdx);
    Eig = sorted_eigenvalues(1);


    eigenvectors = Z'*null(double(A-Eig*B));
    for k = 1:size(eigenvectors, 2)
        max_mag = max(abs(eigenvectors(:, k)));
        eigenvectors(:, k) = eigenvectors(:, k) / max_mag;
    end
    

    half_length = length(eigenvectors) / 2;
    W = eigenvectors(1:half_length);
    Q = eigenvectors(half_length+1:end);
    x = linspace(0, 1, half_length);
    
    figure;
    subplot(2, 1, 1);
    plot(x, real(W), 'LineWidth', 2);
    title('W components of the eigenvector');
    grid on;
    
    subplot(2, 1, 2);
    plot(x, real(Q), 'LineWidth', 2);
    title('Q components of the eigenvector');
    grid on;
end
