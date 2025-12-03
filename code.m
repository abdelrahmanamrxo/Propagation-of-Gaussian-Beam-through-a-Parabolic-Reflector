%% Electromagnetic Waves Project

%% Define Integral Parameters

% Base Constants
lambda = 3.2e-3;
Wo= 10e-3;
k = (2*pi)/lambda;
Zo = (pi*(Wo^2))/lambda;

% Spacial Step Width
dx_min = (sqrt(2)*pi)/k;
dx = dx_min + 0.2e-3;

% Spacial Mesh
L = 40*Wo;
N = floor(L/dx);
x_dist = ( -N/2 : ((N/2)-1) )*dx;
[x, y] = meshgrid(x_dist);

% Spacial Frequency Mesh
dk = (2*pi)/L;
k_dist= ( -N/2 : ((N/2)-1) )*dk;
[kx, ky] = meshgrid(k_dist);

% Define Propagation Constant Kz

kz = k - (kx.^2 + ky.^2) / (2*k);

%% Define U(r) & Plot Initial Density @ Z=0

Ao = 1;
p_sqrd = x.^2 + y.^2;
U_r = Ao*exp(-(p_sqrd)/(Wo^2));

plot_intensity(x_dist, x_dist, U_r, 'Figure 1: Plots @ Z=0');
%% Find Intensity @ Z=0.5Zo & Zo

% First at 0.5Zo
U_1 = propagate_wave(U_r, 0.5*Zo, kz);
plot_intensity(x_dist, x_dist, U_1, 'Figure 2a: Plots @ Z=0.5Zo');

% Next, at Zo
U_2 = propagate_wave(U_r, Zo, kz);
plot_intensity(x_dist, x_dist, U_2, 'Figure 2b: Plots @ Z=Zo');
%% Define the Reflector & Input Field

% Define the reflector in the spacial domain
f_magnitude = 4*Zo;
f = -f_magnitude;
M = exp(-1j*k*p_sqrd/(2*f));
%% Case 1: Reflector @ 3Zo

reflector_position = 3*Zo;
% Find the input field at the reflector position
U_in = propagate_wave(U_r, reflector_position, kz);

% Find the output field from the reflector in the spacial domain
U_out_reflector = U_in .* M;

% Propagate the output field from the reflector for the required cases:

% At Zo
U_z0 = propagate_wave(U_out_reflector, Zo, kz);
plot_intensity(x_dist, x_dist, U_z0, 'Reflector @ 3Zo, Plots @ Zo from Reflector');

% At 4Zo
U_4z0 = propagate_wave(U_out_reflector, 4*Zo, kz);
plot_intensity(x_dist, x_dist, U_4z0, 'Reflector @ 3Zo, Plots @ 4Zo from Reflector');

% At 6Zo
U_6z0 = propagate_wave(U_out_reflector, 6*Zo, kz);
plot_intensity(x_dist, x_dist, U_6z0, 'Reflector @ 3Zo, Plots @ 6Zo from Reflector');
%% Case 2: Reflector @ 4Zo

reflector_position = 4*Zo;
% Find the input field at the reflector position
U_in = propagate_wave(U_r, reflector_position, kz);

% Find the output field from the reflector in the spacial domain
U_out_reflector = U_in .* M;

% Propagate the output field from the reflector for the required cases:

% At Zo
U_z0 = propagate_wave(U_out_reflector, Zo, kz);
plot_intensity(x_dist, x_dist, U_z0, 'Reflector @ 4Zo, Plots @ Zo from Reflector');

% At 4Zo
U_4z0 = propagate_wave(U_out_reflector, 4*Zo, kz);
plot_intensity(x_dist, x_dist, U_4z0, 'Reflector @ 4Zo, Plots @ 4Zo from Reflector');

% At 6Zo
U_6z0 = propagate_wave(U_out_reflector, 6*Zo, kz);
plot_intensity(x_dist, x_dist, U_6z0, 'Reflector @ 4Zo, Plots @ 6Zo from Reflector');
%% Case 3: Reflector @ 5Zo

reflector_position = 5*Zo;
% Find the input field at the reflector position
U_in = propagate_wave(U_r, reflector_position, kz);

% Find the output field from the reflector in the spacial domain
U_out_reflector = U_in .* M;

% Propagate the output field from the reflector for the required cases:

% At Zo
U_z0 = propagate_wave(U_out_reflector, Zo, kz);
plot_intensity(x_dist, x_dist, U_z0, 'Reflector @ 5Zo, Plots @ Zo from Reflector');

% At 4Zo
U_4z0 = propagate_wave(U_out_reflector, 4*Zo, kz);
plot_intensity(x_dist, x_dist, U_4z0, 'Reflector @ 5Zo, Plots @ 4Zo from Reflector');

% At 6Zo
U_6z0 = propagate_wave(U_out_reflector, 6*Zo, kz);
plot_intensity(x_dist, x_dist, U_6z0, 'Reflector @ 5Zo, Plots @ 6Zo from Reflector');
%% Helping Functions

function U_out = propagate_wave(U_in, z_prop, kz)
    U_k = fftshift(fft2(U_in));
    H = exp(-1j * kz * z_prop);
    U_out_k = U_k .* H;
    U_out = ifft2(ifftshift(U_out_k));
end

function plot_intensity(X, Y, U, title_str)
    figure('Name', title_str);
    I = abs(U).^2;

    % Plot the 2D Intensity Heatmap
    subplot(2,1,1);
    imagesc(X*1000, Y*1000, I)
    colormap('jet'); 
    colorbar;
    title(title_str)
    axis image;
    xlabel('x (mm)');
    ylabel('y (mm)');
    xlim([-100, 100]);
    ylim([-100, 100]);

    % Plot the Intensity slice at y=0
    subplot(2,1,2);
    [~, center_y_idx] = min(abs(X - 0)); 
    slice_data = I(center_y_idx, :);
    plot(X*1000, slice_data, 'LineWidth', 2);
    ylabel('Intensity Distribution @ y=0');
    xlabel('x (mm)');
    grid on;
    xlim([-100, 100]);
end