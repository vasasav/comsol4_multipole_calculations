% Written by Vassili Savinov on 10/07/2019
% extract multipoles for mir-ir MM

% on 22/07/19 added the code to compute the correct wavelengts as given by
% giorgio
% how about displacing the orign?

%clear;
close all;

% phyiscal constants
light_speed = 299792458;

lam_range = dlmread('g_wavelengths_SI.csv');%[7.5, 8.0, 8.5]*1e-6;

% get the comsol model
cmsl_model = Get_Bi2Te3_VS_loop_z_offset_good_model(lam_range);

% will assign later on
data_mat = [];

% the geometry is only half along the y-axis this will need to be taken
% into account, but lukily the y=0 is on the symmetry plane so we can
% simpoly fold y back and add.

% the ambient index above the material is 1 (air)

% extract at each wavelength
for i_lam=1:length(lam_range)
    tic;
    
    % specifics of volume integration for this case
    int_callback = @(int_str)(mphint2(cmsl_model, int_str, 'volume', 'selection', [2,5], 'outersolnum', [i_lam], 'dataset', 'dset3'));
    
    % get the datatowe
    light_freq = light_speed/lam_range(i_lam);
    data_row = CSL4v4_ArbMed_ComputeAllMultipoleData_Callback_July19(int_callback, light_freq);
    
    % initialize
    if isempty(data_mat)
        data_mat = zeros(length(lam_range), length(data_row) + 1 + 1 + 1);% wavelength, transmission, reflection
    end
    
    % assign
    data_mat(i_lam, 1) = lam_range(i_lam);
    %
    data_mat(i_lam, 2) = mphglobal(cmsl_model, 'emw.S21', 'outersolnum', [i_lam], 'dataset','dset3');% transmission amplitude
    data_mat(i_lam, 3) = mphglobal(cmsl_model, 'emw.S11', 'outersolnum', [i_lam], 'dataset','dset3');% reflectionamplitude
    % 
    data_mat(i_lam, 4:end) = data_row; % multipoles
    
    dlmwrite('Bi2Te3_Mult_fixed_offset.csv', data_mat);
    
    fprintf('%d/%d (%.3f min)\n', i_lam, length(lam_range), toc/60);
end