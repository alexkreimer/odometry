function varargout = odometry_dialog(varargin)
% ODOMETRY_DIALOG MATLAB code for odometry_dialog.fig
%      ODOMETRY_DIALOG, by itself, creates a new ODOMETRY_DIALOG or raises the existing
%      singleton*.
%
%      H = ODOMETRY_DIALOG returns the handle to a new ODOMETRY_DIALOG or the handle to
%      the existing singleton*.
%
%      ODOMETRY_DIALOG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ODOMETRY_DIALOG.M with the given input arguments.
%
%      ODOMETRY_DIALOG('Property','Value',...) creates a new ODOMETRY_DIALOG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before odometry_dialog_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to odometry_dialog_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help odometry_dialog

% Last Modified by GUIDE v2.5 16-Feb-2016 13:54:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @odometry_dialog_OpeningFcn, ...
                   'gui_OutputFcn',  @odometry_dialog_OutputFcn, ...
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

% --- Executes just before odometry_dialog is made visible.
function odometry_dialog_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to odometry_dialog (see VARARGIN)

% Choose default command line output for odometry_dialog
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using odometry_dialog.
if strcmp(get(hObject,'Visible'),'off')
    plot(rand(5));
end

% UIWAIT makes odometry_dialog wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = odometry_dialog_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
cla;

popup_sel_index = get(handles.popupmenu1, 'Value');
switch popup_sel_index
    case 1
        plot(rand(5));
    case 2
        plot(sin(1:0.01:25.99));
    case 3
        bar(1:.5:10);
    case 4
        plot(membrane);
    case 5
        surf(peaks);
end


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
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

set(hObject, 'String', {'plot(rand(5))', 'plot(sin(1:0.01:25))', 'bar(1:.5:10)', 'plot(membrane)', 'surf(peaks)'});


% --- Executes on selection change in alglist.
function alglist_Callback(hObject, eventdata, handles)
% hObject    handle to alglist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns alglist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from alglist
handle = findobj('Tag', 'listbox2');
contents = cellstr(get(hObject,'String'));
selected = {contents{get(hObject,'Value')}};
D = dir(fullfile('/home','kreimer','KITTI','results',selected{1},'data'));
for i=length(D):-1:1
    if D(i).isdir || D(i).name(1) == '.' 
        D(i) = [];
    end
end

val = cell(length(D),1);
for i=1:length(D)
    [pathstr,name,ext] = fileparts(D(i).name);  
    val{i} = name;
end
set(handle,'String',val);


% --- Executes during object creation, after setting all properties.
function alglist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alglist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

D = dir(fullfile('/home','kreimer','KITTI','results'));
for i=length(D):-1:1
    if ~D(i).isdir || D(i).name(1) == '.' 
        D(i) = [];
    end
end

set(hObject,'String',{D.name});

% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = cellstr(get(handles.alglist,'String'));
algs = contents(get(handles.alglist,'Value'));

contents = cellstr(get(handles.listbox2,'String'));
sequence = contents{get(handles.listbox2,'Value')};

frame = get(handles.edit2,'String');
frame = str2double(frame);
DATA_ROOT  = '/home/kreimer/KITTI/';
KITTI_HOME = fullfile(DATA_ROOT, 'dataset');
RESULT_DIR = fullfile(DATA_ROOT, 'results');
DBG_DIR    = fullfile(DATA_ROOT, 'debug');
image_dir  = fullfile(KITTI_HOME, 'sequences', sequence);

load(fullfile(DATA_ROOT,'data',sequence,sprintf('%06d.mat',frame)));

if length(algs)==1
    images(sequence,frame,algs{1});
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(handles.alglist,'String'));
algs = contents(get(handles.alglist,'Value'));

contents = cellstr(get(handles.listbox2,'String'));
sequence = contents{get(handles.listbox2,'Value')};

frame = get(handles.edit2,'String');
increment = str2double(get(handles.edit3,'String'));
frame = str2double(frame);
frame = frame-increment;
set(handles.edit2,'String',num2str(frame));
plot_data(frame,sequence,algs,handles)

function plot_data(frame,sequence,algs,handles)

DATA_ROOT  = '/home/kreimer/KITTI/';
KITTI_HOME = fullfile(DATA_ROOT, 'dataset');
RESULT_DIR = fullfile(DATA_ROOT, 'results');
DBG_DIR    = fullfile(DATA_ROOT, 'debug');
image_dir  = fullfile(KITTI_HOME, 'sequences', sequence);

load(fullfile(DATA_ROOT,'data',sequence,sprintf('%06d.mat',frame)));
cla(handles.axes3);
hold(handles.axes3,'on');
num_alg = length(algs);

info.gt.poses = util.read_poses(fullfile(DATA_ROOT, 'dataset', 'poses', [sequence, '.txt']));    

for i = 1:num_alg
    alg = algs{i};
    info.(alg).poses = util.read_poses(fullfile(DATA_ROOT, 'results', alg, 'data', [sequence, '.txt']));
    T = [info.(alg).poses(:,:,frame-1);0 0 0 1]\[info.(alg).poses(:,:,frame); 0 0 0 1];
    X0 = [0 0 0 1]';
    X1 = util.h2e(T*X0);
    plot([X0(1);X1(1)],[X0(3);X1(3)],'DisplayName',alg,'parent',handles.axes3);
end

legend(handles.axes3,'show');
axis(handles.axes3,'equal');

num_frames = size(info.(algs{1}).poses,3);
cla(handles.axes4);cla(handles.axes4);
hold(handles.axes4,'on');
poses = info.gt.poses(:,:,1:num_frames);
x = poses(1,4,:);
z = poses(3,4,:);
plot(x(:),z(:),'DisplayName','gt','parent',handles.axes4);
for i=1:num_alg
    field = algs{i};
    poses = info.(field).poses(:,:,1:num_frames);
    x = poses(1,4,:);
    z = poses(3,4,:);
    plot(x(:),z(:),'DisplayName',field, 'parent',handles.axes4);
end
legend(handles.axes4,'show'); axis(handles.axes4,'equal');

for i=1:num_alg
    field = algs{i};
    poses = info.(field).poses(:,:,1:num_frames);
    x = permute(poses(1,4,:),[3 2 1]);
    z = permute(poses(3,4,:),[3 2 1]);
    plot(x(frame-1:frame),z(frame-1:frame),'--gs',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5],...
    'parent',handles.axes4);
end

if num_alg>2
    return;
end

alg1 = algs{1};
alg2 = algs{2};
rot_error = zeros(num_frames,1);
trans_error = zeros(num_frames,1);
abs_rot_error = zeros(num_frames,1);
abs_trans_error = zeros(num_frames,1);

for i = 2:num_frames
    T1 = [info.(alg1).poses(:,:,i-1);0 0 0 1]\[info.(alg1).poses(:,:,i); 0 0 0 1];
    T2 = [info.(alg2).poses(:,:,i-1);0 0 0 1]\[info.(alg2).poses(:,:,i); 0 0 0 1];
    pose_delta = T1\T2;
    rot_error(i) = util.rot_error(pose_delta);
    trans_error(i) = util.trans_error(pose_delta);
    
    pose_error = [info.(alg1).poses(:,:,i);0 0 0 1]\[info.(alg2).poses(:,:,i); 0 0 0 1];
    abs_rot_error(i) = util.rot_error(pose_error);
    abs_trans_error(i) = util.trans_error(pose_error);
end

cla(handles.axes6);
hold(handles.axes6,'on');
plot(rot_error','DisplayName','rotation error','parent',handles.axes6);
legend(handles.axes6,'show');
hold(handles.axes6,'off');

cla(handles.axes7);
hold(handles.axes7,'on');
plot(trans_error','DisplayName','translation error','parent',handles.axes7);
legend(handles.axes7,'show');
hold(handles.axes7,'off');

cla(handles.axes8);
hold(handles.axes8,'on');
plot(abs_rot_error','DisplayName','rotation error','parent',handles.axes8);
legend(handles.axes8,'show');
hold(handles.axes8,'off');

cla(handles.axes9);
hold(handles.axes9,'on');
plot(abs_trans_error','DisplayName','translation error','parent',handles.axes9);
legend(handles.axes9,'show');
hold(handles.axes9,'off');


function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(handles.alglist,'String'));
algs = contents(get(handles.alglist,'Value'));

contents = cellstr(get(handles.listbox2,'String'));
sequence = contents{get(handles.listbox2,'Value')};

frame = get(handles.edit2,'String');
increment = str2double(get(handles.edit3,'String'));
frame = str2double(frame);
frame = frame+increment;
set(handles.edit2,'String',num2str(frame));
plot_data(frame,sequence,algs,handles)
