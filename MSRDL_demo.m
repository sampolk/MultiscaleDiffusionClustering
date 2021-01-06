%% Compare M-SRDL against M-LUND on Salinas A HSI data

profile off;
profile on;

%% Choose whether to save results

prompt = 'Should we save everything? \n 1) Yes \n 2) No\n';
SaveSelected = input(prompt);

if SaveSelected == 1

    save_on = 1;
    
elseif SaveSelected == 2

    save_on = 0;    
    
else
    disp('Incorrect prompt input. Please enter one of [1,2].')
end

%% Choose whether to plot results

prompt = 'Should we plot everything? \n 1) Yes \n 2) No\n';
PlotSelected = input(prompt);

if PlotSelected == 1

    plot_on = 1;
    
elseif PlotSelected == 2

    plot_on = 0;    
    
else
    disp('Incorrect prompt input. Please enter one of [1,2].')
end

%% Run M-LUND & M-SRDL comparisons

[X,Y] = extract_salinasA(); %load data
data_name = 'SalinasA';
load('salinasA-SRDL-HP.mat')

p = KDE(X,Hyperparameters);
G_SpatialRegularization = extract_graph(X, Hyperparameters);
Clusterings_SRDL = M_SRDL(X, Hyperparameters, G_SpatialRegularization, p);

disp('M-SRDL run complete.')

G_NoSpatialRegularization = extract_graph(X, rmfield(Hyperparameters, 'SpatialParams'));
Clusterings_LUND = M_LUND(X, rmfield(Hyperparameters, 'SpatialParams'), G_NoSpatialRegularization, p);
 
disp('M-LUND run complete.')

if save_on 
    save(strcat('M_LUND_Results_', data_name, '.mat'), 'Clusterings_SRDL','Clusterings_LUND', 'X', 'Y', 'Hyperparameters', 'data_name')
end

if plot_on
    results_plot = plot_results_MSRDL(Clusterings_SRDL, Clusterings_LUND);
end


%% 

[~,t] = min(Clusterings_LUND.TotalVI);
[NMI_MLUND] = nmi(Clusterings_LUND.Labels(:,t), Y);

disp('The optimal M-LUND clustering returned an NMI of')
disp(NMI_MLUND)


[~,t] = min(Clusterings_SRDL.TotalVI);
[NMI_MSRDL] = nmi(Clusterings_SRDL.Labels(:,t), Y);

disp('The optimal M-SRDL clustering returned an NMI of')
disp(NMI_MSRDL)

