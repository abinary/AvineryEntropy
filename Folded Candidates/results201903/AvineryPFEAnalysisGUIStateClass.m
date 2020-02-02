classdef AvineryPFEAnalysisGUIStateClass < handle

    properties
        CurrentFolder = [];
        
        EntropyFileList = {};
        MultipleSimulations = true;
        ChooseCandidateWindowsEvenly = true;
        WindowLength = 1000;
        CoarseGrainingLevels = 11;
        WindowStridePercent = 10;
        
        PercentWindowsInAnalysis = 5;
        NumFinalCandidates = 10;
        FoldedCandidateAnalysisOrdering = 3; % 1 - enthalpy, 2 - entropy, 3 - pseudo-free-energy
        
        EntropyData = [];
        
        EnthalpyFilePaths = [];
        EnthalpyData = [];
        ShouldSubtractEnthalpyMean = false;
        
        TrajectoryFileList = {};
        
        UpdateInterfaceFuncs = {};
    end

    methods
        function [obj] = AvineryPFEAnalysisGUIStateClass(obj)
        end
        
        function [] = Save(state)
            my_file_path = mfilename('fullpath');
            [folder, ~, ~] = fileparts(my_file_path);
            
            save([folder filesep 'AvineryPFEAnalysisGUI.state.mat'], 'state');
        end
        
        function [] = Load(obj)
            my_file_path = mfilename('fullpath');
            [folder, ~, ~] = fileparts(my_file_path);
            
            state_file_path = [folder filesep 'AvineryPFEAnalysisGUI.state.mat'];
            
            if (exist(state_file_path, 'file'))
                state = load(state_file_path);
                state = state.state;
                
                % copy fields
                fn = fieldnames(state);
                for i = 1:numel(fn)
                    if (strcmp(fn, 'UpdateInterfaceFuncs')) % skip this one
                        continue;
                    end
                    
                    if (isprop(obj, fn{i}))
                        obj.(fn{i}) = state.(fn{i});
                    end
                end
            end
        end
    end
end
