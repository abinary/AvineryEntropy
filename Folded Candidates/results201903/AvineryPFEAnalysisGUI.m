function varargout = AvineryPFEAnalysisGUI(varargin)
% AVINERYPFEANALYSISGUI MATLAB code for AvineryPFEAnalysisGUI.fig
%      AVINERYPFEANALYSISGUI, by itself, creates a new AVINERYPFEANALYSISGUI or raises the existing
%      singleton*.
%
%      H = AVINERYPFEANALYSISGUI returns the handle to a new AVINERYPFEANALYSISGUI or the handle to
%      the existing singleton*.
%
%      AVINERYPFEANALYSISGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AVINERYPFEANALYSISGUI.M with the given input arguments.
%
%      AVINERYPFEANALYSISGUI('Property','Value',...) creates a new AVINERYPFEANALYSISGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AvineryPFEAnalysisGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AvineryPFEAnalysisGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AvineryPFEAnalysisGUI

% Last Modified by GUIDE v2.5 09-Nov-2019 19:42:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AvineryPFEAnalysisGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @AvineryPFEAnalysisGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before AvineryPFEAnalysisGUI is made visible.
function AvineryPFEAnalysisGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AvineryPFEAnalysisGUI (see VARARGIN)

% Choose default command line output for AvineryPFEAnalysisGUI
handles.output = hObject;
handles.state = AvineryPFEAnalysisGUIStateClass();
handles.state.Load();

% Update handles structure
guidata(hObject, handles);

UpdateInterface(hObject, handles);

% UIWAIT makes AvineryPFEAnalysisGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function [] = GenerateFileListMenu(m, file_list, should_update)

if (nargin < 3)
    should_update = false;
end

if (isempty(file_list))
    return;
end

[prefix, suffix] = FindCommonPreSuf(file_list);
%display([prefix '*' suffix]);

% A misuse of the "fileparts" function, but oh well...
[folder, filename, ext] = fileparts(prefix);
full_prefix = prefix;

expected_num_of_items = numel(file_list) + (~isempty(folder));

if (should_update)
    should_update = (expected_num_of_items == numel(m.Children));
end

if (~isempty(folder))
    prefix = prefix(length(folder) + 2:end);
    folder_text = sprintf('(%s)', folder);

    if (should_update)
        set(m.Children(1), 'Text', folder_text);
    else
        uimenu(m, 'Text', folder_text);
    end
    %set(handles.chosenFilesEditBox, 'String', sprintf('%s*%s (%s)', prefix, suffix, folder));
else
end

varying_part = cellfun(@(c)c(length(full_prefix)+1:end-length(suffix)), file_list, 'UniformOutput', false);
varying_part_num = cellfun(@(c)str2double(c), varying_part);

% all good numbers?
if (all(~isnan(varying_part_num)))
    [~, order] = sort(varying_part_num);
    file_list = file_list(order);
end

for i = 1:numel(file_list)
    f = file_list{i};
    f = f(length(folder) + 2:end);
    
    if (should_update)
        set(m.Children(i + (~isempty(folder))), 'Text', f);
    else
        if (i == 1 && ~isempty(folder))
            uimenu(m, 'Text', f, 'Separator', 'on');
        else
            uimenu(m, 'Text', f);
        end
    end
end

function [valid] = ValidateData(handles, data)
valid = 0;

sliding_data = handles.state.EntropyData;

if (isempty(sliding_data))
    return;
end

window_length = handles.state.WindowLength;
window_stride = round(handles.state.WindowStridePercent * handles.state.WindowLength / 100);
x = cellfun(@(c)numel(c), sliding_data);

data = NormalizeData(handles, data);
ExpectedMeanValues = @(x)floor((numel(x) - window_length + window_stride) / window_stride);

valid = (numel(data) == numel(x)) && all(x == cellfun(@(c)ExpectedMeanValues(c), data));

1;

function [n] = NumWindows(handles)
s = handles.state.EntropyData;

if (isnumeric(s))
    n = numel(s);
else
    n = sum(cellfun(@(c)numel(c), s));
end

function [] = UpdateInterface(hObject, handles)
if (nargin < 2 || isempty(handles))
    handles = guidata(hObject);
end

try
    set(handles.separateSimulationsCheckbox, 'Value', handles.state.MultipleSimulations);
catch
end

try
    set(handles.chooseEvenlyCheckbox, 'Value', handles.state.MultipleSimulations && ...
        handles.state.ChooseCandidateWindowsEvenly);
    
    if (handles.state.MultipleSimulations)
        set(handles.chooseEvenlyCheckbox, 'Enable', 'on');
    else
        set(handles.chooseEvenlyCheckbox, 'Enable', 'off');
    end
catch
end

set(handles.subtractMeanEnthalpyCheckbox, 'Value', handles.state.ShouldSubtractEnthalpyMean);

set(handles.coarseGrainingLevelsEditbox, 'String', num2str(handles.state.CoarseGrainingLevels));
set(handles.windowLengthEditbox, 'String', num2str(handles.state.WindowLength));
set(handles.strideEditBox, 'String', num2str(handles.state.WindowStridePercent));

updateEnthalpyGenericDataDisplayFunc(hObject, handles, ...
    handles.state.EnthalpyFilePaths, ...
    handles.enthalpyDataEditbox, ...
    handles.state.EnthalpyData);

if (~isempty(handles.state.EnthalpyData) && ValidateData(handles, handles.state.EnthalpyData))
    set(handles.enthalpyDataEditbox, 'BackgroundColor', [0.7 1 0.7]);
else
    set(handles.enthalpyDataEditbox, 'BackgroundColor', [1 0.7 0.7]);
end

if (~isempty(handles.state.EnthalpyData))
    set(handles.loadEnthalpyDataBtn, 'String', 'Clear enthalpy data');
else
    set(handles.loadEnthalpyDataBtn, 'String', 'Load enthalpy data ...');
end


updateEnthalpyGenericDataDisplayFunc(hObject, handles, ...
    handles.state.TrajectoryFileList, ...
    handles.cartesianTrajectoriesEditbox, []);

file_list = handles.state.EntropyFileList;

[prefix, suffix] = FindCommonPreSuf(file_list);

[folder, filename, ext] = fileparts(prefix);
full_prefix = prefix;

if (~isempty(folder))
    prefix = prefix(length(folder) + 2:end);
end

if (length(suffix) > 0)
    set(handles.chosenFilesEditBox, 'String', sprintf('%s*%s', prefix, suffix));
else
    set(handles.chosenFilesEditBox, 'String', prefix);
end


m = uicontextmenu();
GenerateFileListMenu(m, file_list);

set(handles.chosenFilesEditBox, 'UIContextMenu', m);

if (numel(handles.state.TrajectoryFileList) == numel(handles.state.EntropyData))
    set(handles.cartesianTrajectoriesEditbox, 'BackgroundColor', [0.7 1 0.7]);
else
    set(handles.cartesianTrajectoriesEditbox, 'BackgroundColor', [1 0.7 0.7]);
end

set(handles.percentWindowsEdit, 'String', num2str(handles.state.PercentWindowsInAnalysis));
set(handles.foldedCandidatesLowestPopup, 'Value', handles.state.FoldedCandidateAnalysisOrdering);

if (~isempty(handles.state.EntropyData))
    n = NumWindows(handles);

    num_of_candidates = round(n * 1e-2 * handles.state.PercentWindowsInAnalysis);
    if (1)
        if (handles.state.ChooseCandidateWindowsEvenly && iscell(handles.state.EntropyData))
            candidates_per_sim = ceil(num_of_candidates / numel(handles.state.EntropyData));
            num_of_candidates = candidates_per_sim * numel(handles.state.EntropyData);
        end
    end
        
    set(handles.numWindowsInCandidateAnalysisEdit, 'String', num2str(num_of_candidates));
end

set(handles.numFinalCandidatesEdit, 'String', num2str(handles.state.NumFinalCandidates));

1;  

function [file_name] = ExtractFileName(file_path)
[~, file_name, ~] = fileparts(file_path);

function [prefix, suffix] = FindCommonPreSuf(string_list)

if (isempty(string_list))
    prefix = '';
    suffix = '';
    return;
end

if (numel(string_list) == 1)
    prefix = string_list{1};
    suffix = '';
    return;
end

lengths = cellfun(@(c)length(c), string_list);
min_length = min(lengths);
max_length = max(lengths);

prefix = string_list{1};
suffix = prefix(end:-1:1);

if (min_length < length(prefix))
    prefix = prefix(1:min_length);
    suffix = suffix(1:min_length);
end

for i = 2:numel(string_list)
    s = string_list{i};
    
    if (length(prefix) > length(s))
        prefix = prefix(1:length(s));
    end
    
    if (length(suffix) > length(s))
        suffix = suffix(1:length(s));
    end
    
    whichNotEqual = ~(lower(prefix) == lower(s(1:numel(prefix))));
    firstNotEqual = find(whichNotEqual, 1);
    
    if (~isempty(firstNotEqual))
        prefix = prefix(1:firstNotEqual-1);
    end

    s = s(end:-1:1); % reverse string for suffix comparison
    whichNotEqual = ~(lower(suffix) == lower(s(1:numel(suffix))));
    firstNotEqual = find(whichNotEqual, 1);
    
    if (~isempty(firstNotEqual))
        suffix = suffix(1:firstNotEqual-1);
    end
end

suffix = suffix(end:-1:1); % reverse suffix back

% This condition means all strings are the same
if (length(suffix) == max_length)
    suffix = '';
end

% --- Outputs from this function are returned to the command line.
function varargout = AvineryPFEAnalysisGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
1;

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function windowLengthEditbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to windowLengthEditbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function strideEditBox_Callback(hObject, eventdata, handles)
% hObject    handle to strideEditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.state.WindowStridePercent = str2double(get(hObject,'String'));
handles.state.Save();
UpdateInterface(hObject, handles);

% --- Executes during object creation, after setting all properties.
function strideEditBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to strideEditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in calculateCoarseGrainingPlotsButton.
function calculateCoarseGrainingPlotsButton_Callback(hObject, eventdata, handles)
% hObject    handle to calculateCoarseGrainingPlotsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
f = gcf();
commandwindow();
file_list = handles.state.FileList;
GenerateCoarseGrainingPlotsOnRandomWindows(2:51, file_list, handles.state.NumOfRandomWindowsForCoarseGrainingOpt, handles.state.WindowLengths);

f_generated = gcf();
figure(f);
figure(f_generated);


function numOfRandomWindowsEditbox_Callback(hObject, eventdata, handles)
% hObject    handle to numOfRandomWindowsEditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.state.NumOfRandomWindowsForCoarseGrainingOpt = str2double(get(hObject,'String'));
handles.state.Save();
UpdateInterface(hObject, handles);

% --- Executes during object creation, after setting all properties.
function numOfRandomWindowsEditbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numOfRandomWindowsEditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function windowLengthsEditbox_Callback(hObject, eventdata, handles)
% hObject    handle to windowLengthsEditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of windowLengthsEditbox as text
%        str2double(get(hObject,'String')) returns contents of windowLengthsEditbox as a double

window_lengths_str = get(hObject,'String');
window_lengths = eval(['[' window_lengths_str ']']);

handles.state.WindowLengths = window_lengths;
handles.state.Save();
UpdateInterface(hObject, handles);


% --- Executes during object creation, after setting all properties.
function windowLengthsEditbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to windowLengthsEditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function chosenFilesEditBox_Callback(hObject, eventdata, handles)
% hObject    handle to chosenFilesEditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of chosenFilesEditBox as text
%        str2double(get(hObject,'String')) returns contents of chosenFilesEditBox as a double


% --- Executes during object creation, after setting all properties.
function chosenFilesEditBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to chosenFilesEditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [fn] = ChooseMultipleFiles(handles, title, filter)
if (isempty(handles.state.CurrentFolder))
    handles.state.CurrentFolder = cd();
end

[fn, pathn] = uigetfile(filter, title, handles.state.CurrentFolder, 'MultiSelect', 'on');

% pressed cancel?
if (pathn == 0)
    fn = {};
    return;
end

handles.state.CurrentFolder = pathn;

if (~iscell(fn))
    fn = {fn};
end

fn = cellfun(@(c)[pathn c], fn, 'UniformOutput', false);


% --- Executes on button press in chooseFilesBtn.
function chooseFilesBtn_Callback(hObject, eventdata, handles)
% hObject    handle to chooseFilesBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fn = ChooseMultipleFiles(handles, 'Choose sliding entropy files...', {...
   'entropy for *.mat;*.entropy.*','Common entropy files (entropy for *.mat,*.entropy.*)'; ...
   '*.mat','MAT-files (*.mat)'; ...
   '*.*',  'All Files (*.*)'});

if (isempty(fn))
    return;
end

LoadEntropyFiles(handles, fn);

handles.state.Save();
UpdateInterface(hObject, handles);

1;

function [] = LoadEntropyFiles(handles, file_list)

if (isempty(file_list))
    handles.state.EntropyFileList = {};
    handles.state.WindowLength = 0;
    handles.state.WindowStridePercent = 0;
    handles.state.CoarseGrainingLevels = 0;
    handles.state.MultipleSimulations = 0;
end

handles.state.EntropyFileList = file_list;

all_data = {};

for file_index = 1:numel(file_list)
    data = load(file_list{file_index});
    all_data{file_index} = data;
    1;
end

handles.state.MultipleSimulations = double(numel(all_data) > 1);
handles.state.WindowLength = all_data{1}.parameters.window_length;
handles.state.WindowStridePercent = all_data{1}.parameters.stride_percent;
handles.state.CoarseGrainingLevels = all_data{1}.parameters.coarse_graining_levels;

handles.state.EntropyData = cellfun(@(c)c.entropy_estimates, all_data, 'UniformOutput', false)
1;

% --- Executes on button press in separateSimulationsCheckbox.
function separateSimulationsCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to separateSimulationsCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.state.FilesAreSeparateSimulations = get(hObject, 'Value');
handles.state.Save();

UpdateInterface(hObject, handles);

% --- Executes on button press in useCacheFilesCheckbox.
function useCacheFilesCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to useCacheFilesCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.state.UseCacheFiles = get(hObject, 'Value');
handles.state.Save();


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over chosenFilesEditBox.
function chosenFilesEditBox_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to chosenFilesEditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
1;

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function coarseGrainingLevelsEditbox_Callback(hObject, eventdata, handles)
% hObject    handle to coarseGrainingLevelsEditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.state.CoarseGrainingLevels = str2double(get(hObject,'String'));
handles.state.Save();

UpdateInterface(hObject, handles);


% --- Executes during object creation, after setting all properties.
function coarseGrainingLevelsEditbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to coarseGrainingLevelsEditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function windowLengthEditbox_Callback(hObject, eventdata, handles)
% hObject    handle to windowLengthEditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.state.WindowLength = str2double(get(hObject,'String'));
handles.state.Save();
UpdateInterface(hObject, handles);

% --- Executes during object creation, after setting all properties.
function windowLengthEditbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to windowLengthEditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function coarseGrainingLevelsEditbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to coarseGrainingLevelsEditbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function windowLengthOptNumWindowsEditbox_Callback(hObject, eventdata, handles)
% hObject    handle to windowLengthOptNumWindowsEditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.state.NumOfRandomWindowsForLengthOpt = str2double(get(hObject,'String'));
handles.state.Save();
UpdateInterface(hObject, handles);


% --- Executes during object creation, after setting all properties.
function windowLengthOptNumWindowsEditbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to windowLengthOptNumWindowsEditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in calculateWindowLenOptPlotsBtn.
function calculateWindowLenOptPlotsBtn_Callback(hObject, eventdata, handles)
% hObject    handle to calculateWindowLenOptPlotsBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file_list = handles.state.FileList;

max_window_len = handles.state.MaxWindowLengthForOpt;
window_lengths = ceil(linspace(50, max_window_len, 31));

f = gcf();
commandwindow();
GenerateWindowLengthOptimizationPlot(handles.state.CoarseGrainingLevels, ...
    file_list, handles.state.NumOfRandomWindowsForLengthOpt, ...
    window_lengths);
f_generated = gcf();
figure(f);
figure(f_generated);



function maxWindowLengthEditbox_Callback(hObject, eventdata, handles)
% hObject    handle to maxWindowLengthEditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.state.MaxWindowLengthForOpt = str2double(get(hObject,'String'));
handles.state.Save();
UpdateInterface(hObject, handles);

% --- Executes during object creation, after setting all properties.
function maxWindowLengthEditbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxWindowLengthEditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [output_file_list] = GenerateOutputFileList(handles)
file_list = handles.state.FileList;
output_file_list = {};

[output_folder, output_name, output_ext] = fileparts(handles.state.SlidingCompressionOutputFileName);

% as a default, set the output folder to the folder of the input
if (isempty(output_folder))
    [output_folder, ~, ~] = fileparts(file_list{1});
end

if (handles.state.FilesAreSeparateSimulations)
    suffix = [output_name output_ext];
    
    [list_prefix, list_suffix] = FindCommonPreSuf(cellfun(@(f)ExtractFileName(f), file_list, 'UniformOutput', false));

    for file_index = 1:numel(file_list)
        filename = ExtractFileName(file_list{file_index});
        
        if (contains(suffix, {'*', '+'}))
            middle = filename(length(list_prefix) + 1:end - length(list_suffix));
            filename = strrep(suffix, '*', filename);
            filename = strrep(filename, '+', middle);
            outputfile = [output_folder filesep filename];
        else
            outputfile = [output_folder filesep filename suffix];
        end
        
        output_file_list{file_index} = outputfile;
    end
else
    if (list_prefix(end) == '.')
        list_prefix(end) = [];
    end
    
    output_name = [output_name output_ext];
    if (contains(output_name, '*'))
        [~, filename, ~] = fileparts(file_list{1});
        outputfile = [output_folder filesep strrep(output_name, '*', prefix)];
    elseif (contains(output_name, '+'))
        [~, filename, ~] = fileparts(file_list{1});
        outputfile = [output_folder filesep strrep(output_name, '+', prefix)];
    else
        outputfile = [output_folder filesep output_name output_ext];
    end
    
    outputfile = {outputfile};
end

% --- Executes on button press in runSlidingCompressionBtn.
function runSlidingCompressionBtn_Callback(hObject, eventdata, handles)
% hObject    handle to runSlidingCompressionBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file_list = handles.state.FileList;

stride = round(handles.state.WindowStridePercent * handles.state.WindowLength / 100);
parameters = struct();
parameters.stride = stride;
parameters.stride_percent = handles.state.WindowStridePercent;
parameters.window_length = handles.state.WindowLength;

f = gcf();
commandwindow(); % go to watch output

[entropy_estimates, sliding_entropy_results] = GenerateSlidingWindowEntropyEstimate(file_list, ...
    false, handles.state.CoarseGrainingLevels, ...
    handles.state.WindowLength, ...
    stride, handles.state.UseWorkerPool);

output_file_list = GenerateOutputFileList(handles);

if (handles.state.FilesAreSeparateSimulations)
    entropy_estimates_all = entropy_estimates;
    sliding_entropy_results_all = sliding_entropy_results;

    for file_index = 1:numel(file_list)
        output_filename = output_file_list{file_index};
                
        entropy_estimates = entropy_estimates_all{file_index};
        sliding_entropy_results = sliding_entropy_results_all{file_index};
        save(output_filename, 'parameters', 'entropy_estimates', 'sliding_entropy_results');
    end
else
    output_filename = output_file_list{1};
    save(output_filename, 'parameters', 'entropy_estimates', 'sliding_entropy_results');
end

% return to app
figure(f);

1;

function [] = updateSlidingCompressionOutputFileEditboxMenu(handles)
%m = get(handles.slidingCompressionOutputFileEditbox, 'UIContextMenu');

%if (isempty(m))
    m = uicontextmenu();
    set(handles.slidingCompressionOutputFileEditbox, 'UIContextMenu', m);
    m.Callback = @slidingCompressionOutputFileEditboxMenu_Callback;
%end

% children = m.Children;
% for c = children
%     c.delete();
% end

output_file_list = GenerateOutputFileList(handles);

GenerateFileListMenu(m, output_file_list, true);

function slidingCompressionOutputFileEditbox_Callback(hObject, eventdata, handles)
% hObject    handle to slidingCompressionOutputFileEditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.state.SlidingCompressionOutputFileName = get(hObject,'String');
handles.state.Save();
UpdateInterface(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slidingCompressionOutputFileEditbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slidingCompressionOutputFileEditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function slidingCompressionOutputFileEditboxMenu_Callback(hObject, eventdata)
handles = guidata(hObject);
slidingCompressionOutputFileEditbox_Callback(handles.slidingCompressionOutputFileEditbox, [], handles);
1;


% --- Executes on button press in useWorkerPoolCheckbox.
function useWorkerPoolCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to useWorkerPoolCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.state.UseWorkerPool = get(hObject, 'Value');
handles.state.Save();
UpdateInterface(hObject, handles);



function edit27_Callback(hObject, eventdata, handles)
% hObject    handle to windowLengthEditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of windowLengthEditbox as text
%        str2double(get(hObject,'String')) returns contents of windowLengthEditbox as a double


% --- Executes during object creation, after setting all properties.
function edit27_CreateFcn(hObject, eventdata, handles)
% hObject    handle to windowLengthEditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit28_Callback(hObject, eventdata, handles)
% hObject    handle to strideEditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of strideEditBox as text
%        str2double(get(hObject,'String')) returns contents of strideEditBox as a double


% --- Executes during object creation, after setting all properties.
function edit28_CreateFcn(hObject, eventdata, handles)
% hObject    handle to strideEditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit29_Callback(hObject, eventdata, handles)
% hObject    handle to coarseGrainingLevelsEditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of coarseGrainingLevelsEditbox as text
%        str2double(get(hObject,'String')) returns contents of coarseGrainingLevelsEditbox as a double


% --- Executes during object creation, after setting all properties.
function edit29_CreateFcn(hObject, eventdata, handles)
% hObject    handle to coarseGrainingLevelsEditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [data] = NormalizeData(handles, data)
% TODO: match given data to sliding entropy data (i.e., multiple vs. single
% simulation)
window_length = handles.state.WindowLength;
window_stride = round(handles.state.WindowStridePercent * handles.state.WindowLength / 100);
sliding_data = handles.state.EntropyData;

%original_lengths = cellfun(@(c)numel(c) * window_stride + (window_length - window_stride), sliding_data);

% Make all data column vectors
data = cellfun(@(c)c(:), data, 'UniformOutput', false);

% concatenate to match the sliding entropy
if (numel(sliding_data) == 1) % a continuous simulation
    data = { vertcat(data{:}) };
end
1;

function [data] = SlidingMean(data, window_length, window_stride)
sliding_avg_kernel = ones(1, window_length) * (1/window_length);

if (isnumeric(data))
    data = conv(data, sliding_avg_kernel, 'valid');
    data = data(1:window_stride:end);
else
    data = cellfun(@(c)conv(c, sliding_avg_kernel, 'valid'), data, 'UniformOutput', false);
    data = cellfun(@(c)c(1:window_stride:end), data, 'UniformOutput', false);
end

function [window_length, window_stride] = GetWindowLengthAndStride(handles)
window_length = handles.state.WindowLength;
window_stride = round(handles.state.WindowStridePercent * handles.state.WindowLength / 100);

function [data] = FetchData(handles, data_name, subtract_mean)
entropy = handles.state.EntropyData;
window_length = handles.state.WindowLength;
window_stride = round(handles.state.WindowStridePercent * handles.state.WindowLength / 100);

if (nargin < 3)
    subtract_mean = false;
end

switch (data_name)
    case 'Entropy'
        data = entropy;
        
    case 'Enthalpy'
        data = handles.state.EnthalpyData;
        
        if (~isempty(data))
            data = NormalizeData(handles, data);
            data = SlidingMean(data, window_length, window_stride);
        end

    case 'PFE'
        h = FetchData(handles, 'Enthalpy', false);
        s = FetchData(handles, 'Entropy', false);
        
        if (isnumeric(h))
            data = h - s;
        else
            data = arrayfun(@(i)h{i} - s{i}, 1:numel(h), 'UniformOutput', false);
        end
        
    case 'Random'
        data = entropy;
        
        if (isnumeric(data))
            data = rand(size(data));
        else
            data = cellfun(@(x)rand(size(x)), data, 'UniformOutput', false);
        end

    case 'Sim'
        window_length = handles.state.WindowLength;
        window_stride = round(handles.state.WindowStridePercent * handles.state.WindowLength / 100);
        
        if (isnumeric(entropy))
            data = ones(size(PFE));
        else
            data = arrayfun(@(i)i * ones(size(entropy{i})), 1:numel(entropy), 'UniformOutput', false);
        end
        
    case 'WindowStart'
        window_length = handles.state.WindowLength;
        window_stride = round(handles.state.WindowStridePercent * handles.state.WindowLength / 100);
        
        if (isnumeric(entropy))
            data = 1 + window_stride * reshape([1:numel(PFE)] - 1, size(PFE));
        else
            data = cellfun(@(c)1 + ([1:numel(c)] - 1) * window_stride, entropy, 'UniformOutput', false);
        end

    otherwise
        data = [];
end

if (subtract_mean)
    if (iscell(data))
        data = cellfun(@(c)c - mean(c), data, 'UniformOutput', false);
    elseif (isnumeric(data))
        data = data - mean(data);
    end
end
1;

% --- Executes on button press in plot_S_H_Scatter_Btn.
function plot_S_H_Scatter_Btn_Callback(hObject, eventdata, handles)
% hObject    handle to plot_S_H_Scatter_Btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
x = FetchData(handles, 'Enthalpy', handles.state.ShouldSubtractEnthalpyMean);
y = FetchData(handles, 'Entropy', false);

figure(1);
clf;
[plot_handles] = PlotData(x, y, 'H (k_B \cdot T)', 'S (k_B)');
hold off;
axis equal

if (1)
    for i = 1:numel(plot_handles)
        s = plot_handles(i);
        set(s, 'UserData', handles.state.EntropyFileList{i});
    end
    
    datacursormode on;
    dcm_obj = datacursormode(gcf());
    dcm_obj.UpdateFcn = @PolyChargedDataCursorFunc;
end

function enthalpyDataEditbox_Callback(hObject, eventdata, handles)
% hObject    handle to enthalpyDataEditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of enthalpyDataEditbox as text
%        str2double(get(hObject,'String')) returns contents of enthalpyDataEditbox as a double


% --- Executes during object creation, after setting all properties.
function enthalpyDataEditbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to enthalpyDataEditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit33_Callback(hObject, eventdata, handles)
% hObject    handle to edit33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit33 as text
%        str2double(get(hObject,'String')) returns contents of edit33 as a double


% --- Executes during object creation, after setting all properties.
function edit33_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox10.
function checkbox10_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox10



function edit37_Callback(hObject, eventdata, handles)
% hObject    handle to edit37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit37 as text
%        str2double(get(hObject,'String')) returns contents of edit37 as a double


% --- Executes during object creation, after setting all properties.
function edit37_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit36_Callback(hObject, eventdata, handles)
% hObject    handle to edit36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit36 as text
%        str2double(get(hObject,'String')) returns contents of edit36 as a double


% --- Executes during object creation, after setting all properties.
function edit36_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit35_Callback(hObject, eventdata, handles)
% hObject    handle to edit35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit35 as text
%        str2double(get(hObject,'String')) returns contents of edit35 as a double


% --- Executes during object creation, after setting all properties.
function edit35_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cartesianTrajectoriesEditbox_Callback(hObject, eventdata, handles)
% hObject    handle to cartesianTrajectoriesEditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cartesianTrajectoriesEditbox as text
%        str2double(get(hObject,'String')) returns contents of cartesianTrajectoriesEditbox as a double


% --- Executes during object creation, after setting all properties.
function cartesianTrajectoriesEditbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cartesianTrajectoriesEditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function numWindowsInCandidateAnalysisEdit_Callback(hObject, eventdata, handles)
% hObject    handle to numWindowsInCandidateAnalysisEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
n = NumWindows(handles);
num_windows_in_analysis = str2double(get(hObject,'String'));
handles.state.PercentWindowsInAnalysis = 100 * num_windows_in_analysis / n;
handles.state.Save();
UpdateInterface(hObject, handles);


% --- Executes during object creation, after setting all properties.
function numWindowsInCandidateAnalysisEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numWindowsInCandidateAnalysisEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function [] = updateEnthalpyGenericDataDisplayFunc(hObject, handles, file_list, editbox, data)
[prefix, suffix] = FindCommonPreSuf(file_list);

[folder, filename, ext] = fileparts(prefix);
full_prefix = prefix;

if (~isempty(folder))
    prefix = prefix(length(folder) + 2:end);
end

if (length(suffix) > 0)
    display_text = sprintf('%s*%s', prefix, suffix);
else
    display_text = prefix;
end

if (~isempty(data))
    if (iscell(data))
        num_values_read = cellfun(@(c)numel(c), data);
        num_values_read = sum(num_values_read);
    else
        num_values_read = numel(data);
    end
    
    display_text = sprintf('%s (%d values)', display_text, num_values_read);
end

set(editbox, 'String', display_text);

m = uicontextmenu();
GenerateFileListMenu(m, file_list);
set(editbox, 'UIContextMenu', m);
1;

function [] = updateEnthalpyDisplayFunc(hObject, handles)
updateEnthalpyGenericDataDisplayFunc(hObject, handles, ...
    handles.state.EnthalpyFilePaths, ...
    handles.enthalpyDataEditbox, ...
    handles.state.EnthalpyData);
1;

% --- Executes on button press in loadEnthalpyDataBtn.
function loadEnthalpyDataBtn_Callback(hObject, eventdata, handles)
% hObject    handle to loadEnthalpyDataBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (~isempty(handles.state.EnthalpyData))
    handles.state.EnthalpyData = [];
    handles.state.EnthalpyFilePaths = [];
else
    fn = ChooseMultipleFiles(handles, 'Choose sliding entropy files...', {...
   '*.energy.*;*.pe.*','Common enthalpy files (*.energy.*,*.pe.*)'; ...
   '*.mat','MAT-files (*.mat)'; ...
   '*.*',  'All Files (*.*)'});
    if (isempty(fn))
        return;
    end
    
    handles.state.EnthalpyFilePaths = fn;
    
    data = cellfun(@(f)LoadMatrixFromFile(f), fn, 'UniformOutput', false);
    handles.state.EnthalpyData = data;
end

handles.state.Save();
UpdateInterface(hObject, handles);


% --- Executes during object creation, after setting all properties.
function loadEnthalpyDataBtn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loadEnthalpyDataBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function additionalEntropyEditbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to additionalEntropyEditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes during object creation, after setting all properties.
function additionalEnthalpyEditbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to additionalEnthalpyEditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function additionalEnthalpyEditbox_Callback(hObject, eventdata, handles)
% hObject    handle to additionalEnthalpyEditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% hObject    handle to loadEnthalpyDataBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in additionalEntropyBtn.
function additionalEntropyBtn_Callback(hObject, eventdata, handles)
% hObject    handle to additionalEntropyBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (~isempty(handles.state.AdditionalEntropyData))
    handles.state.AdditionalEntropyData = [];
    handles.state.AdditionalEntropyFilePaths = [];
else
    fn = ChooseMultipleFiles(handles, 'Choose additional entropy data files...', '*.*');
    if (isempty(fn))
        return;
    end
    
    handles.state.AdditionalEntropyFilePaths = fn;
    
    data = cellfun(@(f)LoadMatrixFromFile(f), fn, 'UniformOutput', false);
    handles.state.AdditionalEntropyData = data;
end

handles.state.Save();
UpdateInterface(hObject, handles);


% --- Executes on button press in additionalEnthalpyBtn.
function additionalEnthalpyBtn_Callback(hObject, eventdata, handles)
% hObject    handle to additionalEnthalpyBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (~isempty(handles.state.AdditionalEnthalpyData))
    handles.state.AdditionalEnthalpyData = [];
    handles.state.AdditionalEnthalpyFilePaths = [];
else
    fn = ChooseMultipleFiles(handles, 'Choose additional enthalpy data files...', '*.*');
    if (isempty(fn))
        return;
    end
    
    handles.state.AdditionalEnthalpyFilePaths = fn;
    
    data = cellfun(@(f)LoadMatrixFromFile(f), fn, 'UniformOutput', false);
    handles.state.AdditionalEnthalpyData = data;
    
end

handles.state.Save();
UpdateInterface(hObject, handles);



function additionalEntropyEditbox_Callback(hObject, eventdata, handles)
% hObject    handle to additionalEntropyEditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of additionalEntropyEditbox as text
%        str2double(get(hObject,'String')) returns contents of additionalEntropyEditbox as a double
% --- Executes on button press in selectCartesianTrajBtn.


function selectCartesianTrajBtn_Callback(hObject, eventdata, handles)
% hObject    handle to selectCartesianTrajBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fn = ChooseMultipleFiles(handles, 'Select cartesian trajectory files ...', '*.*');
if (isempty(fn))
    return;
end

handles.state.TrajectoryFileList = [];

display('Trying to load first trajectory file... ');
try
    data = LoadCartesianTrajectoryFile(fn{1}); % try loading one file
    handles.state.TrajectoryFileList = fn;
    display('Success!');
catch
    display('FAILED LOADING');
    fn = {};
end

handles.state.Save();
UpdateInterface(hObject, handles);


function edit41_Callback(hObject, eventdata, handles)
% hObject    handle to edit41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit41 as text
%        str2double(get(hObject,'String')) returns contents of edit41 as a double


% --- Executes during object creation, after setting all properties.
function edit41_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit42_Callback(hObject, eventdata, handles)
% hObject    handle to edit42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit42 as text
%        str2double(get(hObject,'String')) returns contents of edit42 as a double


% --- Executes during object creation, after setting all properties.
function edit42_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [x] = ColumnizeCells(x)
x = cellfun(@(c)c(:), x, 'UniformOutput', false);

function [x] = VCatCells(x)
x = vertcat(x{:});

% --- Executes on button press in saveDataBtn.
function saveDataBtn_Callback(hObject, eventdata, handles)
% hObject    handle to saveDataBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, path, filterindex] = uiputfile(...
    {'*.txt', 'Tab separated text file (*.txt)'; '*.mat', 'Matlab file (*.mat)'}, ...
    'Choose file path to save to', handles.state.CurrentFolder);

if (filename)
    handles.state.CurrentFolder = path;
    handles.state.Save();
    
    path = [path filename];
    
    entropy = FetchData(handles, 'Entropy', false);
    enthalpy = FetchData(handles, 'Enthalpy', handles.state.ShouldSubtractEnthalpyMean);
    
    window_length = handles.state.WindowLength;
    window_stride = round(handles.state.WindowStridePercent * handles.state.WindowLength / 100);

    if (isnumeric(enthalpy))
        %PFE =  enthalpy - entropy;
        sim = ones(size(PFE));
        windows_start = 1 + window_stride * reshape([1:numel(PFE)] - 1, size(PFE));
    else
        assert(iscell(enthalpy));
        %PFE = arrayfun(@(i)enthalpy{i} - entropy{i}, 1:numel(enthalpy), 'UniformOutput', false);
        sim = arrayfun(@(i)i * ones(size(enthalpy{i})), 1:numel(enthalpy), 'UniformOutput', false);
        windows_start = cellfun(@(c)1 + ([1:numel(c)] - 1) * window_stride, sim, 'UniformOutput', false);
    end
        
    if (~isnumeric(enthalpy))
        sim = VCatCells(ColumnizeCells(sim));
        windows_start = VCatCells(ColumnizeCells(windows_start));
        enthalpy = VCatCells(ColumnizeCells(enthalpy));
        entropy = VCatCells(ColumnizeCells(entropy));
        %PFE = VCatCells(ColumnizeCells(PFE));
    end
    
    PFE =  enthalpy - entropy;
    
    switch (filterindex)
        case 1
            m = [sim, windows_start, enthalpy, entropy, PFE];
            
            f = fopen(path, 'w');
            
            if (f > 0)
                fprintf(f, 'Sim\tWindowStart\tEnthalpy\tEntropy\tPFE\n');
                %fprintf(f, '%d\t%d\t%g\t%g\t%g\n', m');
                fprintf(f, '%d\t%d\t%e\t%e\t%e\n', m');
                fclose(f);
            end
            
        case 2
            save(path, 'sim', 'windows_start', 'entropy', 'enthalpy', 'PFE');
            
    end
end


% --- Executes on button press in plot_S_H_Density_Btn.
function plot_S_H_Density_Btn_Callback(hObject, eventdata, handles)
% hObject    handle to plot_S_H_Density_Btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

x = FetchData(handles, 'Enthalpy', handles.state.ShouldSubtractEnthalpyMean);
y = FetchData(handles, 'Entropy', false);

ColumnizeCells = @(x)cellfun(@(c)c(:), x, 'UniformOutput', false);
VCatCells = @(x)vertcat(x{:});

if (~isnumeric(x))
    x = VCatCells(ColumnizeCells(x));
    y = VCatCells(ColumnizeCells(y));
end

[n, x_edges, y_edges] = histcounts2(x, y, 51);

x_centers = 0.5 * (x_edges(1:end-1) + x_edges(2:end));
y_centers = 0.5 * (y_edges(1:end-1) + y_edges(2:end));

figure(1);
imagesc(x_centers, y_centers, n');
axis xy equal



function numFinalCandidatesEdit_Callback(hObject, eventdata, handles)
% hObject    handle to numFinalCandidatesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.state.NumFinalCandidates = str2double(get(hObject,'String'));
handles.state.Save();
UpdateInterface(hObject, handles);


% --- Executes during object creation, after setting all properties.
function numFinalCandidatesEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numFinalCandidatesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chooseEvenlyCheckbox.
function chooseEvenlyCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to chooseEvenlyCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.state.ChooseCandidateWindowsEvenly = get(hObject, 'Value');
handles.state.Save();
UpdateInterface(hObject, handles);

% --- Executes on button press in subtractMeanEnthalpyCheckbox.
function subtractMeanEnthalpyCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to subtractMeanEnthalpyCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.state.ShouldSubtractEnthalpyMean = get(hObject, 'Value');
handles.state.Save();
UpdateInterface(hObject, handles);


function percentWindowsEdit_Callback(hObject, eventdata, handles)
% hObject    handle to percentWindowsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.state.PercentWindowsInAnalysis = str2double(get(hObject,'String'));
handles.state.Save();
UpdateInterface(hObject, handles);


% --- Executes during object creation, after setting all properties.
function percentWindowsEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to percentWindowsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in foldedCandidatesLowestPopup.
function foldedCandidatesLowestPopup_Callback(hObject, eventdata, handles)
% hObject    handle to foldedCandidatesLowestPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.state.FoldedCandidateAnalysisOrdering = get(hObject,'Value');
handles.state.Save();
UpdateInterface(hObject, handles);




% --- Executes during object creation, after setting all properties.
function foldedCandidatesLowestPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to foldedCandidatesLowestPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in goFoldedCandidateAnalysisBtn.
function goFoldedCandidateAnalysisBtn_Callback(hObject, eventdata, handles)
% hObject    handle to goFoldedCandidateAnalysisBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[window_length, window_stride] = GetWindowLengthAndStride(handles);

ordering = handles.state.FoldedCandidateAnalysisOrdering;

enthalpy = FetchData(handles, 'Enthalpy', handles.state.ShouldSubtractEnthalpyMean);
entropy = FetchData(handles, 'Entropy', false);

if (isempty(enthalpy))
    assert(ordering == 2, 'Enthalpy data is missing!');
    enthalpy = cellfun(@(c)c .* 0, entropy, 'UniformOutput', false);
end

pfe = FetchData(handles, 'PFE');

switch(ordering)
    case 1
        ordering_key = enthalpy;
        
    case 2
        ordering_key = entropy;
        
    case 3
        ordering_key = pfe;
        
    case 4
        ordering_key = FetchData(handles, 'Random');
end

window_start = FetchData(handles, 'WindowStart');
sim = FetchData(handles, 'Sim');

n = NumWindows(handles);
num_of_candidates = round(n * 1e-2 * handles.state.PercentWindowsInAnalysis);
candidate_windows = [];

if (iscell(ordering_key))
    if (handles.state.ChooseCandidateWindowsEvenly)
        candidates_per_sim = ceil(num_of_candidates / numel(ordering_key));
        num_of_candidates = candidates_per_sim * numel(ordering_key);
        
        for s = 1:numel(ordering_key)
            [~, order] = sort(ordering_key{s});
            w = order(1:candidates_per_sim);
            
            candidate_windows = [candidate_windows; ...
                s + zeros(candidates_per_sim, 1), window_start{s}(w)', ...
                enthalpy{s}(w)', entropy{s}(w)', pfe{s}(w)'];
        end
    else
        ordering_key = VCatCells(ColumnizeCells(ordering_key));
        sim = VCatCells(ColumnizeCells(sim));
        window_start = VCatCells(ColumnizeCells(window_start));

        pfe = VCatCells(ColumnizeCells(pfe));
        entropy = VCatCells(ColumnizeCells(entropy));
        enthalpy = VCatCells(ColumnizeCells(enthalpy));
    end
end

if (isnumeric(ordering_key))
    [~, order] = sort(ordering_key);
    w = order(1:num_of_candidates);
    %candidate_windows = [sim(w) window_start(w)];
    candidate_windows = [sim(w) window_start(w) enthalpy(w) entropy(w) pfe(w)];
end

options = FoldedCandidateAnalysis();
options.NumOfFinalists = handles.state.NumFinalCandidates;

time_stamp_str = datestr(now, 'yyyymmdd-hhMMss');

results = FoldedCandidateAnalysis(handles.state.TrajectoryFileList, ...
    candidate_windows, window_length, options);

save(sprintf('folded-analysis %s.results.mat', time_stamp_str), 'results');
