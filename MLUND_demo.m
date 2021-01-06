%%

profile off;
profile on;

%% Choose the dataset

prompt = 'Which dataset? \n 1) 3D Gaussians \n 2) 2D Nonlinear \n 3) 2D Bottleneck \n 4) Salinas A HSI';
DataSelected = input(prompt);

if DataSelected == 1
    
    [X,Y] = gaussian_sample(1,3,1000,1000);
    data_name = 'Gaussians';
    load('gaussian-HP.mat')

elseif DataSelected ==2
    
    [X,Y] = nonlinear_sample(0.2, 1.2, 4, 180, 1200, 4000);
    data_name = 'Nonlinear';
    load('nonlinear-HP.mat')
    
elseif DataSelected == 3 
    
    [X,Y] = bottleneck_sample(1500, 1500, 0.35,  1, 5, 12);
    data_name = 'Bottleneck';
    load('bottleneck-HP.mat')
    
elseif DataSelected == 4
    
    [X,Y] = extract_salinasA();
    data_name = 'SalinasA';
    load('salinasA-HP.mat')
    
else
    disp('Incorrect prompt input. Please enter one of [1,2,3,4].')
end

%% Choose whether to save results

prompt = 'Should we save everything? \n 1) Yes \n 2) No';
SaveSelected = input(prompt);

if SaveSelected == 1

    save_on = 1;
    
elseif SaveSelected == 2

    save_on = 0;    
    
else
    disp('Incorrect prompt input. Please enter one of [1,2].')
end

%% Choose whether to plot results

prompt = 'Should we plot everything? \n 1) Yes \n 2) No';
PlotSelected = input(prompt);

if PlotSelected == 1

    plot_on = 1;
    
    % Choose whether to plot stochastic complements
    prompt = 'Should we plot the intervals? These are computationally expensive. \n 1) Yes \n 2) No';
    PlotSelected = input(prompt);

    if PlotSelected == 1

        sc_on = 1;

    elseif PlotSelected == 2

        sc_on = 0;    

    else
        disp('Incorrect prompt input. Please enter one of [1,2].')
    end
    
    
    
elseif PlotSelected == 2

    plot_on = 0;    
    
else
    disp('Incorrect prompt input. Please enter one of [1,2].')
end





%% Run M-LUND

Clusterings = M_LUND(X, Hyperparameters);

if save_on 
    save(strcat('M_LUND_Results_', data_name, '.mat'), 'Clusterings', 'X', 'Y', 'Hyperparameters', 'data_name')
end

if plot_on
    results_plot = plot_results(X, Clusterings, data_name, sc_on);
end