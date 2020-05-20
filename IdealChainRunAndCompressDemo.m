

% Load .Net engines
NET.addAssembly([cd '\Avinery.Sim.dll']);
NET.addAssembly([cd '\Avinery.Compression.dll']);

%% Setup

num_monomers = 100;
random_seed = 1;

R = 0;
sampling_interval = 1000;
num_of_samples = 1024;
steps = num_of_samples * sampling_interval;
number_of_burned_steps = 10*sampling_interval;

%% Simulation parameters
minimal_move_space = 3;
maximal_move_space = num_monomers / 2;
probability_of_moving_edges = 0.0; % 0 for stretched polymer

mc = Avinery.Sim.IdealChainMC();
mc.RandomSeed = random_seed;
mc.ShouldSampleConfigurations = true;

mc.MinimalMoveInterval = minimal_move_space;
mc.MaximalMoveInterval = maximal_move_space;
mc.N = num_monomers;
mc.ProbabilityOfMovingEdges = probability_of_moving_edges;
mc.InitializeConfiguration(R);

if (0)
    %%
    coordinates = double(mc.ConfigurationMatrix);
    
    figure(1);
    plot(coordinates(:, 1), coordinates(:, 2), '-*k');
    %xlim(1.5*[-R R])
    %ylim(xlim());
    hold on;
    plot([0 R], [0 0], '--r');
    hold off;
    axis equal;
end

mc.Burn(number_of_burned_steps);

%%
icLog = mc.Run(steps, sampling_interval);
icLog.WriteBinary('ideal_chain_configurations.binary');

%%
f = fopen('ideal_chain_configurations.binary');
configurations = fread(f, '*double');
fclose(f);

configurations = reshape(configurations, [2 num_monomers size(configurations, 1)/(2*num_monomers)]);
d = diff(configurations, 1, 2);
angles = squeeze(atan2(d(2, :, :), d(1, :, :)));

coarse_graining_level = 11;

discretized_angles = floor((angles + pi) * (coarse_graining_level / (2*pi)));
% just in case (due to numerical innacuracy)
discretized_angles(discretized_angles < 0) = 0;
discretized_angles(discretized_angles >= coarse_graining_level) = coarse_graining_level - 1;

compressed_size = double(Avinery.Compression.LZMA.CompressedSize(discretized_angles(:)));
compression_ratio = compressed_size / numel(discretized_angles);

%% scale
zero_compressed_size = double(Avinery.Compression.LZMA.CompressedSize(zeros([numel(discretized_angles) 1])));
zero_compressed_ratio = zero_compressed_size / numel(discretized_angles);

random_compressed_size_list = ...
    arrayfun(@(i)double(Avinery.Compression.LZMA.CompressedSize(randi(coarse_graining_level, [numel(discretized_angles) 1]))), ...
    1:10);

random_compressed_size = mean(random_compressed_size_list);

random_compression_ratio = random_compressed_size / numel(discretized_angles);

adjusted_incompressibility = (compression_ratio - zero_compressed_ratio) / (random_compression_ratio - zero_compressed_ratio);
entropy_estimate = size(discretized_angles, 1) * log(coarse_graining_level) * (adjusted_incompressibility - 1); % in k_B
entropy_estimate2 = size(discretized_angles, 1) * log2(coarse_graining_level) * (adjusted_incompressibility - 1); % in bits
