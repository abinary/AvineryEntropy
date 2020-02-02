classdef StateClass < handle
    properties
        WindowLength
        WindowStride
        Data
        SelectedIndex
        Trajectories = {}
        TrajectoriesRama = {}
        Loaders
    end
end
