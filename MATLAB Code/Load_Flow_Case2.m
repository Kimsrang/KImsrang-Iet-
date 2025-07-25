clc; clear all; 

% Twisted Cable Specification
cablespec = [
 % No   CSA   Dt     R      Ir @ 40degree
    1   16   17.5   1.91   38;
    2   25   20     1.20   51;
    3   35   23     0.868  63;
    4   50   27     0.641  89;
    5   70   31.5   0.443 145;
    6   95   38     0.320/3 160;
    7  120   41.5   0.253 185;
    8  150   44.5   0.206 214
];

f = 50;       % Frequency [Hz]
rho = 100;    % Soil resistivity [Ohm·m]

% Cable indices for each feeder (you can modify these)
n1 = 1; % Feeder 1st cable number
n2 = 4; % Feeder 2nd cable number
n3 = 6; % Feeder 3rd cable number (3x1x4C-95mm2 Al/XLPE)
n4 = 1; % Feeder 4th cable number

feeder_labels = {'1st', '2nd', '3rd', '4th'};
cable_indices = [n1 n2 n3 n4];

for i = 1:4
    cableNo = cable_indices(i);
    CSA = cablespec(cableNo,2);
    Dt = cablespec(cableNo,3);
    R = cablespec(cableNo,4);

    r = (Dt / 4) / 1000; % radius in meters
    GMR = r * exp(-1/4);
    d_mn = sqrt((2*r)^2 + (2*r)^2); % distance between conductors in meters

    % Calculate impedance matrix Z
    Z = complex(zeros(3,3));
    for m = 1:3
        for n = 1:3
            if m == n
                Z(m,n) = R + 0.00158836*f + 1i*0.00202237*f*(log(1/GMR) + 7.6786 + 0.5*log(rho/f));
            else
                Z(m,n) = 0.00158836*f + 1i*0.00202237*f*(log(1/d_mn) + 7.6786 + 0.5*log(rho/f));
            end
        end
    end

    if i == 1
        Z1 = Z; % cable for feeder 1
    elseif i == 2
        Z2 = Z; % cable for feeder 2
    elseif i == 3
        Z3 = Z; % cable for feeder 3
    elseif i == 4
        Z4 = Z; % cable for feeder 4
    end
    % fprintf('\nCable for feeder %s:\n', feeder_labels{i});
    % fprintf('  Cross-sectional Area (CSA): %d mm^2\n', CSA);
    % fprintf('  GMR: %.3f mm\n', GMR*1000);
    % fprintf('  Impedance Matrix Z (Ohm/km):\n');
    % disp(Z);
end


%%% Prepare BusData
data = readtable('BusBranchData.xlsx', 'Sheet',  'Clustering Strategy');

% Type of Bus
Slack = 1; % Slack bus
PQ = 2; % Load bus
Bus = data.Bus; 
Type = data.Type; % slack bus = 1 and PQ bus = 2
Pa = data.P_A_kW_; % power for each bus on phase A [kW]
Pb = data.P_B_kW_; % power for each bus on phase B [kW]
Pc = data.P_C_kW_; % power for each bus on phase C [kW]

Databus = table(Bus(1:14), Type(1:14), Pa(1:14), Pb(1:14), Pc(1:14), ...
    'VariableNames', {'Bus', 'Type', 'Pa_kW', 'Pb_kW', 'Pc_kW'});

%disp(DataBus(1:47,:));

%%% Prepare Branch Data
fb = data.FromBus; % frombus
tb = data.ToBus;   % tobus
L = data.Length_m_; % distance between buses
feeder = data.ithFeeder; % ith feeder in network

% Prepare branch data table first
BranchData= table(fb(1:13), tb(1:13), L(1:13), feeder(1:13), ...
    'VariableNames', {'FromBus', 'ToBus', 'L[m]', 'ith Feeder'});

% Preallocate a cell array to hold the total impedance matrices for each branch
Z_total = cell(height(BranchData), 1);

% Loop through each row to multiply length with correct impedance matrix
for i = 1:height(BranchData)
    length_km = BranchData.('L[m]')(i) / 1000;
    feederID = BranchData.('ith Feeder')(i);

    switch feederID
        case 1
            Z_total{i} = Z1 * length_km;
        case 2
            Z_total{i} = Z2 * length_km;
        case 3
            Z_total{i} = Z3 * length_km;
        case 4
            Z_total{i} = Z4 * length_km;
    end
end

% Add the computed impedance matrices back to the BranchData table
BranchData.Z = Z_total(1:13); % size 3*3 for each branch [0hm]

%disp(BranchData(1:46,:))

% Impedance for each branch
nBranches = height(BranchData);
Z = zeros(3,3,nBranches);
for k = 1:nBranches
    if isempty(BranchData.Z{k})
        error('Branch %d has empty impedance matrix', k);
    else
        Z(:,:,k) = BranchData.Z{k};
    end
end

%% Backward/Forward Sweep Method 
% --- Base values ---
Sbase = 560e3; % [VA]
Vbase = 400 / sqrt(3); % [V]
Zbase = Vbase^2 / Sbase; % [ohm]

nBranches = height(BranchData);
nBus = height(Databus);

% --- Prepare Impedance Matrix ---
Z = zeros(3,3,nBranches);
for k = 1:nBranches
    if isempty(BranchData.Z{k})
        error('Branch %d has empty impedance matrix', k);
    else
        Z(:,:,k) = BranchData.Z{k};
    end
end

% --- Load Powers (in pu) ---
pf = 0.95; theta = acos(pf);
Sload = zeros(3, nBus);
for i = 1:nBus
    Sload(1,i) = Databus.Pa_kW(i) * 1000 / Sbase * (1 + 1j*tan(theta)); % [pu]
    Sload(2,i) = Databus.Pb_kW(i) * 1000 / Sbase * (1 + 1j*tan(theta)); % [pu]
    Sload(3,i) = Databus.Pc_kW(i) * 1000 / Sbase * (1 + 1j*tan(theta)); % [pu]
end

% --- Initialization ---
Vbus = ones(3, nBus);
Vbus(:,1) = [1; exp(-1j*2*pi/3); exp(1j*2*pi/3)];  % Slack bus at bus 1

max_iter = 100;
tol = 1e-6;

for iter = 1:max_iter
    V_prev = Vbus;

    %% === BACKWARD SWEEP ===
    k1 = 0.96; % correction factor @ 35 degree celcius
    Iinj = conj(Sload ./ Vbus)./k1;  % 3 x nBus complex current injections

    Ibranch = zeros(3, nBranches);  % Branch currents

    for k = nBranches:-1:1
        toBus = BranchData.ToBus(k);
        fromBus = BranchData.FromBus(k);

        I = Iinj(:,toBus); % [pu]

        % Add downstream branch currents
        downstreamBranches = find(BranchData.FromBus == toBus);
        for j = 1:length(downstreamBranches)
            I = I + Ibranch(:,downstreamBranches(j));
        end

        Ibranch(:,k) = I;
    end

    %% === FORWARD SWEEP ===
    Vbus(:,1) = [1; exp(-1j*2*pi/3); exp(1j*2*pi/3)];  % Reset slack

    for k = 1:nBranches
        fromBus = BranchData.FromBus(k);
        toBus = BranchData.ToBus(k);

        Zpu = Z(:,:,k) / Zbase;
        Vbus(:,toBus) = Vbus(:,fromBus) - Zpu * Ibranch(:,k); % [pu]
    end

    %% === Check convergence ===
    deltaV = max(max(abs(abs(Vbus) - abs(V_prev))));
    if deltaV < tol
        fprintf('Converged in %d iterations.\n', iter);
        break;
    end
end


%% ---Calculation Power Losse by I^2 (3x3) .* R (3x3) per branch ---
nBranches    = size(Z,3);
Ploss_W       = zeros(3,3,nBranches);   % element-wise I^2 .* R (in W)
P_loss_phase = zeros(3,nBranches);     % per-phase loss (sum rows)
P_loss_total = zeros(1,nBranches);     % per-branch total loss

for k = 1:nBranches
    Rk  = real(Z(:,:,k));        % 3x3 (pu)
    Iph = Ibranch(:,k);          % 3x1 (pu)

    % Element-wise I^2 .* R, then to Watts
    Ploss_W(:,:,k) = abs(Iph.^2 .* Rk * Sbase);

    % Sum each row -> per-phase loss (A,B,C) including your I^2*R row-scaling
    P_loss_phase(:,k) = sum(Ploss_W(:,:,k), 2);

    % Total branch loss
    P_loss_total(k)  = sum(Ploss_W(:,:,k), 'all');
end

%% --- Table (per branch) ---
BranchFrom = BranchData.FromBus(:);
BranchTo   = BranchData.ToBus(:);

Lineloss = table(BranchFrom, BranchTo, ...
          P_loss_phase(1,:)', P_loss_phase(2,:)', P_loss_phase(3,:)', ...
          'VariableNames', {'FromBus','ToBus', ...
                            'P_loss_A_W','P_loss_B_W','P_loss_C_W'});
disp(Lineloss);

%% --- Totals ---
total_loss_A = sum(P_loss_phase(1,:));
total_loss_B = sum(P_loss_phase(2,:));
total_loss_C = sum(P_loss_phase(3,:));
total_loss   = sum(P_loss_total);

TotalLosses = table( ...
    {'Total Loss Phase A'; 'Total Loss Phase B'; 'Total Loss Phase C'; 'Total System Loss'}, ...
    [total_loss_A; total_loss_B; total_loss_C; total_loss], ...
    'VariableNames', {'Description','Loss_W'});
disp(TotalLosses);

%%% Extract to Excel 
% writetable(Lineloss, 'Result_case2.xlsx','Sheet', 'LineLoss');


% Voltage Profile 
nBus = size(Vbus, 2); % number of poles (buses)

% Preallocate cell arrays for formatted polar strings
Va_str = cell(nBus,1);
Vb_str = cell(nBus,1);
Vc_str = cell(nBus,1);

% Fill strings with "magnitude ∠ angle°"
for i = 1:nBus
    Va = Vbus(1,i);
    Vb = Vbus(2,i);
    Vc = Vbus(3,i);

    Va_str{i} = sprintf('%.3f ∠ %6.2f°', abs(Va), angle(Va)*180/pi);
    Vb_str{i} = sprintf('%.3f ∠ %6.2f°', abs(Vb), angle(Vb)*180/pi);
    Vc_str{i} = sprintf('%.3f ∠ %6.2f°', abs(Vc), angle(Vc)*180/pi);
end

% Create the table
Pole_No = (1:nBus)';
VoltageProfile = table(Pole_No, Va_str, Vb_str, Vc_str, ...
    'VariableNames', {'Pole_No', 'Va_deg', 'Vb_deg', 'Vc_deg'});

% Display the table
disp(VoltageProfile)

% writetable(VoltageProfile, 'Result_Case2.xlsx', 'Sheet', 'VoltageProfile');


% --- Plot Voltages ---
figure;
bus_indices = 1:nBus;
plot(bus_indices, abs(Vbus(1,:)), 'r-o', 'LineWidth', 1.5); hold on;
plot(bus_indices, abs(Vbus(2,:)), 'g--s', 'LineWidth', 1.5);
plot(bus_indices, abs(Vbus(3,:)), 'b-.*', 'LineWidth', 1.5);
yline(0.9, '--k', '0.9 p.u. limit', 'LineWidth', 1.2);
xlabel('Bus Number');
ylabel('Voltage Magnitude (p.u.)');
title('Voltage Profile for Each Bus');
legend('Phase A', 'Phase B', 'Phase C');
grid on;
