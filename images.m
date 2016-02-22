function varargout = images(varargin)
% IMAGES MATLAB code for images.fig
%      IMAGES, by itself, creates a new IMAGES or raises the existing
%      singleton*.
%
%      H = IMAGES returns the handle to a new IMAGES or the handle to
%      the existing singleton*.
%
%      IMAGES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMAGES.M with the given input arguments.
%
%      IMAGES('Property','Value',...) creates a new IMAGES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before images_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to images_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help images

% Last Modified by GUIDE v2.5 09-Feb-2016 11:50:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @images_OpeningFcn, ...
                   'gui_OutputFcn',  @images_OutputFcn, ...
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


% --- Executes just before images is made visible.
function images_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to images (see VARARGIN)

% Choose default command line output for images
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
sequence = varargin{1};
DATA_ROOT  = '/home/kreimer/KITTI/';
KITTI_HOME = fullfile(DATA_ROOT, 'dataset');
image_dir  = fullfile(KITTI_HOME, 'sequences', sequence);
frame = varargin{2};
load(fullfile(DATA_ROOT,'data',sequence,sprintf('%06d.mat',frame)));
[i1, i2] = read_images(image_dir, frame);
[i1p,i2p]= read_images(image_dir, frame-1);
image = [i1 i2; i1p i2p];
imshow(image,'parent',handles.axes1);

x2(1,:)  = size(i1,2) + x2(1,:);
x1p(2,:) = size(i1,1) + x1p(2,:);
x2p(1,:) = size(i1,2) + x2p(1,:);
x2p(2,:) = size(i1,1) + x2p(2,:);

alg = varargin{3};
inliers = st.(alg).inliers{1};
hold on;
plot([x1(1,inliers); x2(1,inliers); x2p(1,inliers); x1p(1,inliers); x1(1,inliers)],...
    [x1(2,inliers); x2(2,inliers); x2p(2,inliers); x1p(2,inliers); x1(2,inliers)]);


% UIWAIT makes images wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = images_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
