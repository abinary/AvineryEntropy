classdef AvinerySlidingEntropyGUIStateClass < handle

    properties
        CurrentFolder = [];
        FileList = {};
        FilesAreSeparateSimulations = true;
        UseCacheFiles = true;
        WindowLengths = [1000];
        NumOfRandomWindowsForCoarseGrainingOpt = 5;
        CoarseGrainingLevels = 11;
        WindowLength = 1000;
        WindowStridePercent = 10;
        NumOfRandomWindowsForLengthOpt = 1e2;
        MaxWindowLengthForOpt = 2000;
        SlidingCompressionOutputFileName = '*.running_entropy.mat';
        
        CoarseGrainingLevelsToTest = [];
        
        UseWorkerPool = false;
    end

    methods
        function [obj] = AvinerySlidingEntropyGUIStateClass(obj)
        end
        
        function [] = Save(state)
            my_file_path = mfilename('fullpath');
            [folder, ~, ~] = fileparts(my_file_path);
            
            save([folder filesep 'AvinerySlidingEntropyGUI.state.mat'], 'state');
        end
        
        function [] = Load(obj)
            my_file_path = mfilename('fullpath');
            [folder, ~, ~] = fileparts(my_file_path);
            
            state_file_path = [folder filesep 'AvinerySlidingEntropyGUI.state.mat'];
            
            if (exist(state_file_path, 'file'))
                state = load(state_file_path);
                state = state.state;
                
                % copy fields
                fn = fieldnames(state);
                for i = 1:numel(fn)
                    if (isprop(obj, fn{i}))
                        obj.(fn{i}) = state.(fn{i});
                    end
                end
            end
        end
    end
end
