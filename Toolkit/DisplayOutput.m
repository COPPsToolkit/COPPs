function [] = DisplayOutput(params, exitflag, Lopt)

disp('------------------------------------------------------------------')
if exitflag == 1
    disp('                          COPPs converged!                        ')
    disp('------------------------------------------------------------------')
elseif exitflag == 0
    warning('Maximum number of iterations exceeded')
elseif exitflag == -2
    warning('No feasible point found')
elseif exitflag == -3
    warning('Problem is unbounded')
else
    warning('Algorithm not converged')
end

fprintf('\nModel: %s\n\n', params.ModelFilename)

fprintf('Loss function: Weight x variables\n')
for LossVars_index = 1:length(params.LossFunctionVariables)
    fprintf('               %6.2f x %s\n', params.LossFunctionWeights(LossVars_index), params.LossFunctionVariables{LossVars_index})
end

fprintf('\nPolicy variables: name \t shock\n')
for Instr_index = 1:size(params.PolicyInstrumentsAndShocks,1)
    fprintf('                  %s \t %s\n', char(params.PolicyInstrumentsAndShocks(Instr_index)), char(params.PolicyInstrumentsAndShocks{Instr_index,2}))
end

if strcmp(params.PolicyType,'COM') && strcmp(params.Constraint,'OFF')
    fprintf('\nType of policy: Optimal commitment\n')
elseif strcmp(params.PolicyType,'DIS') && strcmp(params.Constraint,'OFF')
    fprintf('\nType of policy: Optimal discretion\n')
elseif strcmp(params.PolicyType,'COM') && strcmp(params.Constraint,'ON')
    fprintf('\nType of policy: Optimal constrained commitment\n')
elseif strcmp(params.PolicyType,'DIS') && strcmp(params.Constraint,'ON')
    fprintf('\nType of policy: Optimal constrained discretion\n')
end

if strcmp(params.Constraint,'ON')
    fprintf('\nConstraints\n')
    fprintf('\t ELB:    %6.2f\n', params.RateMin)
    fprintf('\t QE Min: %6.2f\n', params.QeMin)
    fprintf('\t QE Max: %6.2f\n', params.QeMax)
    if isfield(params.UserProvided,'RateDevMin')
            fprintf(['\n\t RateDevMin: ', num2str(params.UserProvided.RateDevMin)])
    end
    if isfield(params.UserProvided,'RateDevMax')
            fprintf(['\n\t RateDevMax: ', num2str(params.UserProvided.RateDevMax)])
    end
    if isfield(params.UserProvided,'QeDevMin')
            fprintf(['\n\t QeDevMin: ', num2str(params.UserProvided.QeDevMin)])
    end
    if isfield(params.UserProvided,'QeDevMin')
            fprintf(['\ns\t QeDevMin: ', num2str(params.UserProvided.QeDevMin)])
    end
    
end

disp(' ')
disp('-------------------------------------------------')
disp('Optional parameters')
disp(params.UserProvided)
disp(' ')

disp('-------------------------------------------------')
disp('Optimal loss relative to baseline')
fprintf('\t \t \t %.2f%%\n',100*Lopt);
disp(' ')
disp('------------------------------------------------------------------')
