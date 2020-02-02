classdef TrajectoriesClass < handle
    properties
        Data
        SelectedIndex
        Trajectories = {}
        TrajectoriesRama = {}
    end

    properties (SetAccess = private)
        Loaders
        WindowLength
        WindowStride
    end
    
    methods
        function [obj] = TrajectoriesClass(loaders, windowLength, windowStride)
            obj.Loaders = loaders;
            obj.WindowLength = windowLength;
            obj.WindowStride = windowStride;
        end
        
        function [frames] = GetCAFrames(obj, simNum, frame_indexes)
            if (numel(obj.Trajectories) < (simNum + 1) || isempty(obj.Trajectories{simNum + 1}))
                obj.Trajectories = {}; % clear up memory
                obj.Trajectories{simNum + 1} = obj.Loaders.LoadTrajectoryFileCAlpha(simNum);
            end

            frames = obj.Trajectories{simNum + 1}.frames;
            frames = frames(frame_indexes, :, :);
        end
        
        function [frames] = GetWindowCAFrames(obj, simNum, startFrame, subset)
            if (numel(obj.Trajectories) < (simNum + 1) || isempty(obj.Trajectories{simNum + 1}))
                obj.Trajectories{simNum + 1} = obj.Loaders.LoadTrajectoryFileCAlpha(simNum);
            end

            frames = obj.Trajectories{simNum + 1}.frames;
            framesIndexes = (startFrame - 1) + [1:obj.WindowLength];
            
            if (nargin >= 4 && ~isempty(subset))
                frames = frames(framesIndexes(subset), :, :);
            else
                frames = frames(framesIndexes, :, :);
            end
        end
    end
end
