function varargout = AvinerySlidingEntropyGUI(varargin)
% AVINERYSLIDINGENTROPYGUI MATLAB code for AvinerySlidingEntropyGUI.fig
%      AVINERYSLIDINGENTROPYGUI, by itself, creates a new AVINERYSLIDINGENTROPYGUI or raises the existing
%      singleton*.
%
%      H = AVINERYSLIDINGENTROPYGUI returns the handle to a new AVINERYSLIDINGENTROPYGUI or the handle to
%      the existing singleton*.
%
%      AVINERYSLIDINGENTROPYGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AVINERYSLIDINGENTROPYGUI.M with the given input arguments.
%
%      AVINERYSLIDINGENTROPYGUI('Property','Value',...) creates a new AVINERYSLIDINGENTROPYGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AvinerySlidingEntropyGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AvinerySlidingEntropyGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AvinerySlidingEntropyGUI

% Last Modified by GUIDE v2.5 19-Nov-2019 10:39:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AvinerySlidingEntropyGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @AvinerySlidingEntropyGUI_OutputFcn, ...
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


% --- Executes just before AvinerySlidingEntropyGUI is made visible.
function AvinerySlidingEntropyGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AvinerySlidingEntropyGUI (see VARARGIN)

% Choose default command line output for AvinerySlidingEntropyGUI
handles.output = hObject;
handles.state = AvinerySlidingEntropyGUIStateClass();
try
    handles.state.Load();
catch err
    display(err.message());
end

% Update handles structure
guidata(hObject, handles);

% Link related edit boxes
set(handles.coarseGrainingLevelsEditbox2, 'Callback', get(handles.coarseGrainingLevelsEditbox, 'Callback'));
set(handles.windowLengthEditbox2, 'Callback', get(handles.windowLengthEditbox, 'Callback'));

UpdateInterface(hObject, handles);

% UIWAIT makes AvinerySlidingEntropyGUI wait for user response (see UIRESUME)
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
 
function [] = UpdateInterface(hObject, handles)
if (nargin < 2 || isempty(handles))
    handles = guidata(hObject);
end

set(handles.separateSimulationsCheckbox, 'Value', handles.state.FilesAreSeparateSimulations);
set(handles.useCacheFilesCheckbox, 'Value', handles.state.UseCacheFiles);

if (handles.state.FilesAreSeparateSimulations)
    set(handles.slidingCompressionOutputFileLabel, 'String', 'Output files suffix:');
else
    set(handles.slidingCompressionOutputFileLabel, 'String', 'Output file:');
end

window_lengths = arrayfun(@(x)num2str(x), handles.state.WindowLengths, 'UniformOutput', false);

window_lengths = [ '[' strjoin(window_lengths, ',') ']' ];
set(handles.windowLengthsEditbox, 'String', window_lengths);

set(handles.numOfRandomWindowsEditbox, 'String', num2str(handles.state.NumOfRandomWindowsForCoarseGrainingOpt));

updateLevelsToTestEditbox(hObject, handles);

set(handles.coarseGrainingLevelsEditbox, 'String', num2str(handles.state.CoarseGrainingLevels));
set(handles.coarseGrainingLevelsEditbox2, 'String', num2str(handles.state.CoarseGrainingLevels));

set(handles.windowLengthEditbox, 'String', num2str(handles.state.WindowLength));
set(handles.windowLengthEditbox2, 'String', num2str(handles.state.WindowLength));
set(handles.strideEditBox, 'String', num2str(handles.state.WindowStridePercent));

set(handles.slidingCompressionOutputFileEditbox, 'String', handles.state.SlidingCompressionOutputFileName);
updateSlidingCompressionOutputFileEditboxMenu(handles);

set(handles.useWorkerPoolCheckbox, 'Value', handles.state.UseWorkerPool);

set(handles.windowLengthOptNumWindowsEditbox, 'String', num2str(handles.state.NumOfRandomWindowsForLengthOpt));
set(handles.maxWindowLengthEditbox, 'String', num2str(handles.state.MaxWindowLengthForOpt));

file_list = handles.state.FileList;

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
function varargout = AvinerySlidingEntropyGUI_OutputFcn(hObject, eventdata, handles) 
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
GenerateCoarseGrainingPlotsOnRandomWindows(handles.state.CoarseGrainingLevelsToTest, ...
    file_list, handles.state.NumOfRandomWindowsForCoarseGrainingOpt, handles.state.WindowLengths, 1234); % fix the lottery
%    file_list, handles.state.NumOfRandomWindowsForCoarseGrainingOpt, handles.state.WindowLengths);

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
   '*.rama.*;*.angles.*','Common dihedral angle files (*.rama.*,*.angles.*)'; ...
   '*.mat','MAT-files (*.mat)'; ...
   '*.*',  'All Files (*.*)'});
if (isempty(fn))
    return;
end

handles.state.FileList = fn;
handles.state.Save();

UpdateInterface(hObject, handles);

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
    window_lengths, 1234); % fix the lottery
%    window_lengths);
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

if (isempty(file_list))
    return;
end

[output_folder, output_name, output_ext] = fileparts(handles.state.SlidingCompressionOutputFileName);

% as a default, set the output folder to the folder of the input
if (isempty(output_folder))
    [output_folder, ~, ~] = fileparts(file_list{1});
end

[list_prefix, list_suffix] = FindCommonPreSuf(cellfun(@(f)ExtractFileName(f), file_list, 'UniformOutput', false));

if (handles.state.FilesAreSeparateSimulations)
    suffix = [output_name output_ext];
    
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
        outputfile = [output_folder filesep strrep(output_name, '*', list_prefix)];
    elseif (contains(output_name, '+'))
        [~, filename, ~] = fileparts(file_list{1});
        outputfile = [output_folder filesep strrep(output_name, '+', list_prefix)];
    else
        outputfile = [output_folder filesep output_name output_ext];
    end
    
    output_file_list = {outputfile};
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
parameters.coarse_graining_levels = handles.state.CoarseGrainingLevels;

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

function updateLevelsToTestEditbox(hObject, handles)
l = handles.state.CoarseGrainingLevelsToTest;

if (~isempty(l))
    dl = mean(diff(l));
    
    l1 = min(l):dl:max(l);
    
    if (all(l1 == l))
        if (dl == 1)
            s = sprintf('%d:%d', min(l), max(l));
        else
            s = sprintf('%d:%d:%d', min(l), dl, max(l));
        end
    else
        s = arrayfun(@(x)num2str(x), l, 'UniformOutput', false);
        s = strjoin(s, ', ');
    end
    
    set(handles.levelsToTestEditbox, 'String', s);
end

function levelsToTestEditbox_Callback(hObject, eventdata, handles)
% hObject    handle to levelsToTestEditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
levels_str = get(hObject,'String');
levels = eval(['[' levels_str ']']);

handles.state.CoarseGrainingLevelsToTest = levels;
handles.state.Save();
UpdateInterface(hObject, handles);




% --- Executes during object creation, after setting all properties.
function levelsToTestEditbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to levelsToTestEditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
