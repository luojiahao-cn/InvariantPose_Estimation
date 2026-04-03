function params = get_params()
% GET_PARAMS Get simplified experiment parameter configuration

%% Magnet parameters (3 magnets)
params.magnet.m_hat = [
    [-sin(1/12*pi);0;cos(1/12*pi)], ...
    [sin(1/12*pi);0;cos(1/12*pi)],...
    [0; sin(1/12*pi); cos(1/12*pi)]
    ];
params.magnet.m_pos = [
    [-0.1; 2e-4; 0.3],...
    [0.1; 2e-4; 0.3],...
    [0; -0.1; 0.3]
    ];
params.magnet.m_norm = [500, 500, 500];
params.magnet.m_hat = params.magnet.m_hat ./ vecnorm(params.magnet.m_hat);

%% Sensor parameters (7 sensors in cross configuration)
params.sensor.d_list = [
    [0; 0; 0], ...
    [1e-3; 0; 0], ...
    [-1e-3; 0; 0], ...
    [0; 1e-3; 0], ...
    [0; -1e-3; 0], ...
    [0; 0; 1e-3], ...
    [0; 0; -1e-3]
    ];

%% Uncertainty parameters
params.uncertainty.noise_strength = 0e-4;

%% Optimization parameters
params.optimization.options = optimoptions('lsqnonlin', ...
    'Algorithm', 'levenberg-marquardt', ...
    'Display', 'off');
params.optimization.mu = 1e2;
params.optimization.beta = 1e2;
params.optimization.W = diag([1e-4, 1, 1]);

%% Workspace constraints
params.workspace.center = [0; 0; 0];
params.workspace.radius = 0.15;

end
