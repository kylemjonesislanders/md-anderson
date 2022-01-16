%% Clear Workspace
clear all; clc; close all;

%% Load Image Data
load('ImageData.mat');
% CEST Data
images_preinjection = ROIavg_PreInjectionSimulatedSpectrums;
images_postinjection = ROIavg_PostInjectionCESTSpectrums;
SatFrqList = SatFrqList*2*pi;

%% Remove Dummy Scans from PostInjection Images
num_dummyfreq=length(SatFrqList(SatFrqList<-10000*2*pi)); % Determines which CEST datapoints are dummy scans
images_postinjection(:,:,1:num_dummyfreq)=[]; % Removes dummy scans from postinjection CEST spectrums
SatFrqList(1:num_dummyfreq)=[]; % Removes Sat Frequencies that correspond to dummy scans

%% Define number of samples to analyze and initialize vector to hold pH results
num_phantoms = 24;
pH_results = zeros(num_phantoms, 1);

for phantom = 1:num_phantoms
    
    %% Isolate single Z spectrum
    phantom_num = phantom;
    ZSignal_preinjection = squeeze(images_preinjection(1, phantom_num, :));
    ZSignal_postinjection = squeeze(images_postinjection(1, phantom_num, :));

    %% Normalize Z spectra by the setting the maximum to 1
    ZSignal_preinjection=ZSignal_preinjection./max(ZSignal_preinjection);
    ZSignal_postinjection=ZSignal_postinjection./max(ZSignal_postinjection);

    %% Load Chemical Data file
    load('chemical_data.mat');
    
    %% Feed in experimental Z signal and Saturation Frequency List into chemical data file
    chemical_data.experiment_list(1).Z_signal = ZSignal_preinjection; 
    chemical_data.experiment_list(2).Z_signal = ZSignal_postinjection; 
    chemical_data.experiment_list(1).Hz_offset = SatFrqList;
    chemical_data.experiment_list(2).Hz_offset = SatFrqList;
    
    %% Convert the saturation pulse power to nutation rate and save value in the waveform list of the chemical_data file
    chemical_data.waveform_list(1,2) = 2*pi*SatPower*42.58; % Saturation Power in Nutation Rate

    %% Make sure to update chemical_data file with parameter map info
    % Add R1 and R2 values (for water pool) to chemical_data file if you have T1 and T2 maps
    if (chemical_data.spinsystem_list(1).fit_R1 == false)
        T1 = (T1map.avg_final(phantom_num))./1000;
        R1 = 1/T1;
        chemical_data.spinsystem_list(1).R1 = R1;
    end
    if (chemical_data.spinsystem_list(1).fit_R2 == false)
       T2 = (T2map.avg_final(phantom_num))./1000;
       R2 = 1/T2;
       chemical_data.spinsystem_list(1).R2 = R2;
    end
    
    % Add B1 correction to chemical data file if you have a B1 Map
    if (chemical_data.fit_B1_correction == false)
        B1_correction = ROIavg_B1(phantom_num)/100;
        chemical_data.B1_correction = B1_correction;
    end

    % Add B0 shift information into chemical_data file if you have a B0 map. If not, estimate B0 shift with the offset of the water signal.
    if (exist('B0map', 'var')==1)
        B0_shift = B0map(phantom_num)*MHz_value*2*pi;
        chemical_data.B0_shift = B0_shift;
    else
        water_minimum_signal_index_number=find(ZSignal_preinjection==min(ZSignal_preinjection));
        water_minimum_signal_index_number=water_minimum_signal_index_number(1);
        chemical_data.B0_shift = SatFrqList(water_minimum_signal_index_number);
    end

    %% Create initial guesses, lowerbound guesses, and upperbound guesses for parameter that will be fit
    [parameters, lower_bound, upper_bound, typical_parameters, chemical_data] = initialize_parameters(chemical_data);

    %% Perform Bloch Fitting
    options_lsqnonlin = optimoptions('lsqnonlin', 'MaxFunEvals', 300000, 'MaxIter', 3000, 'TypicalX', typical_parameters, 'TolX', 1e-7, 'TolFun', 1e-7);%TLi changed the MaxFunEvals, MaxIter, TolX and TolFun on 20200607
    [final_parameters, resnorm, R, exitflag, output, lambda, J] = lsqnonlin(@(parameters) lsqnonlin_contrast_renorm(chemical_data, parameters), parameters, lower_bound, upper_bound, options_lsqnonlin);

    %% Update chemical data file with final parameters
    [chemical_data_final] = bloch_spectra_model_function(final_parameters, chemical_data);
    fitted_pH_val = chemical_data_final.experiment_list(2).pH;
    pH_results(phantom, 1) = fitted_pH_val;

    %% Plot fitting results
    % Fitted Difference Spectrum
    figure(1)
    raw_difference_spectrum = (chemical_data_final.experiment_list(1).Z_signal-chemical_data_final.experiment_list(2).Z_signal)';
    plot(chemical_data_final.experiment_list(2).Hz_offset, raw_difference_spectrum, 'o') 
    hold on;
    simulated_difference_spectrum = chemical_data_final.experiment_list(1).simulated_Z_signal-chemical_data_final.experiment_list(2).simulated_Z_signal;
    plot(chemical_data_final.experiment_list(2).Hz_offset, simulated_difference_spectrum)
    x_coord = min(chemical_data_final.experiment_list(2).Hz_offset);
    y_coord = min(simulated_difference_spectrum)+ (max(simulated_difference_spectrum) - min(simulated_difference_spectrum))/2;
    text(x_coord, y_coord, append("Fitted pH: ", num2str(fitted_pH_val)));
    text(x_coord, y_coord+.025, append("Phantom number: ", num2str(phantom)));
    title("Fitted Difference Spectrum")
    hold off;
    pause(5)
end