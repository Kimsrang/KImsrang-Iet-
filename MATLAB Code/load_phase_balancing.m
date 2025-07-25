function [phases, phaseCurrents, deviations, UC] = load_phase_balancing(filename)

    % Constants
    pf = 0.95;     % Power factor
    V1 = 230;      % Single-phase voltage (V)
    V3 = 400;      % Three-phase line voltage (V)

    % Load data
    data = readtable(filename);
    data = data(1:min(129, height(data)), :);

    % Extract power values and phase types
    loadValues = data.P_kW__Pf_0_95;  % kW
    phaseType = data.Phase;           % 1 = single-phase, 3 = three-phase
    x = data.X;
    y = data.Y;
    p = [x, y];

    numLoads = length(loadValues);
    phases = zeros(1, numLoads);      
    phaseCurrents = zeros(1, 3);      % A, B, C
    currentList = zeros(numLoads,1);  % Current per load

    % Assign loads
    for i = 1:numLoads
        P = loadValues(i) * 1000; % kW -> W
        if phaseType(i) == 1  % Single-phase
            I = P / (V1 * pf);
            [~, minIdx] = min(phaseCurrents);
            phases(i) = minIdx;  % 1 = A, 2 = B, 3 = C
            phaseCurrents(minIdx) = phaseCurrents(minIdx) + I;
        elseif phaseType(i) == 3  % Three-phase
            I = P / (sqrt(3) * V3 * pf);
            phaseCurrents = phaseCurrents + I;
            phases(i) = 0;  % Three-phase
        else
            error('Unknown phase type at row %d: %d', i, phaseType(i));
        end
        currentList(i) = I;
    end

    % Check balance
    I_avg = mean(phaseCurrents);
    deviations = abs(phaseCurrents - I_avg);
    UC = (1/3) * sum((phaseCurrents / I_avg).^2);

    % Output
    fprintf('\n⚡ Phase Currents Including Three-Phase Loads:\n');
    fprintf('Phase A: %.3f A\n', phaseCurrents(1));
    fprintf('Phase B: %.3f A\n', phaseCurrents(2));
    fprintf('Phase C: %.3f A\n', phaseCurrents(3));
    fprintf('Average Current: %.3f A\n', I_avg);
    fprintf('Deviation from average: [%.3f, %.3f, %.3f] A\n', deviations);
    fprintf('Unbalance Coefficient (UC): %.4f\n', UC);

    if UC < 1.1
        fprintf('✅ Load is balanced (UC < 1.1).\n');
    else
        fprintf('⚠️ Load is unbalanced (UC ≥ 1.1).\n');
    end

    % Prepare table for export
    phaseStr = strings(numLoads,1);
    for i = 1:numLoads
        if phases(i) == 0
            phaseStr(i) = "Three-Phase";
        elseif phases(i) == 1
            phaseStr(i) = "Phase A";
        elseif phases(i) == 2
            phaseStr(i) = "Phase B";
        elseif phases(i) == 3
            phaseStr(i) = "Phase C";
        end
    end

    resultTable = table((1:numLoads).', p,  loadValues, currentList, phaseType, phases.', phaseStr, ...
        'VariableNames', {'LoadIndex', 'Coordinate', 'Power_kW', 'Current_A', 'PhaseCode', 'AssignedPhaseNum', 'AssignedPhase'});

    % Write results to Excel
    % writetable(resultTable, 'Phase_Assignments.xlsx', 'Sheet', 1);

    % Add summary at the end
    summaryTable = table( ...
        phaseCurrents(1), phaseCurrents(2), phaseCurrents(3), I_avg, UC, ...
        'VariableNames', {'PhaseA_A', 'PhaseB_A', 'PhaseC_A', 'Average_A', 'UC'});
    % writetable(summaryTable, 'Phase_Assignments.xlsx', 'Sheet', 2);

    % Detailed Load Assignment (printed)
    fprintf('\n--- Load Assignments ---\n');
    for i = 1:numLoads
        fprintf('Load %3d: %.3f kW (%.2f A) -> %s\n', ...
            i, loadValues(i), currentList(i), phaseStr(i));
    end
    % fprintf('\nResults saved to: Phase_Assignments.xlsx\n');
end
