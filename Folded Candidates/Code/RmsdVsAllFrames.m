function [rmsd] = RmsdVsAllFrames(state, coordinates)

if (numel(state.Trajectories) < 27)
    state.Trajectories{27} = [];
end

for i = 1:27
    if (isempty(state.Trajectories{i}))
        fprintf('Loading trajectory #%d \n', i);
        state.Trajectories{i} = state.Loaders.LoadTrajectoryFileCAlpha(i);
    end
end

if (0)
    display('Calculating RMSD for all frames');
    
    all_frames = arrayfun(@(i)state.Trajectories{i}.frames, 1:numel(state.Trajectories), 'UniformOutput', false);
    all_frames = vertcat(all_frames{:});
    
    rmsd2 = arrayfun(@(i)...
        CalculateRmsd(representative_frame, squeeze(all_frames(i, :, :))), ...
        1:size(all_frames, 1));
else
    rmsd2 = [];
end
end
