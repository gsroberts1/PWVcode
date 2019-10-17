function varargout = Anatomical_Images(varargin)
% ANATOMICAL_IMAGES MATLAB code for Anatomical_Images.fig
%      ANATOMICAL_IMAGES, by itself, creates a new ANATOMICAL_IMAGES or raises the existing
%      singleton*.
%
%      H = ANATOMICAL_IMAGES returns the handle to a new ANATOMICAL_IMAGES or the handle to
%      the existing singleton*.
%
%      ANATOMICAL_IMAGES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ANATOMICAL_IMAGES.M with the given input arguments.
%
%      ANATOMICAL_IMAGES('Property','Value',...) creates a new ANATOMICAL_IMAGES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Anatomical_Images_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Anatomical_Images_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Anatomical_Images

% Last Modified by GUIDE v2.5 11-Oct-2019 17:53:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Anatomical_Images_OpeningFcn, ...
                   'gui_OutputFcn',  @Anatomical_Images_OutputFcn, ...
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


% --- Executes just before Anatomical_Images is made visible.
function Anatomical_Images_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Anatomical_Images (see VARARGIN)

% Choose default command line output for Anatomical_Images
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Anatomical_Images wait for user response (see UIRESUME)
% uiwait(handles.figure1);

handles.anatDatasets = varargin{1};
set(handles.DatasetListbox,'String',{handles.anatDatasets.Names});
guidata(hObject, handles);
updateImages(handles)
global spanUD spanLR
spanUD = 0;
spanLR = 0;



% --- Outputs from this function are returned to the command line.
function varargout = Anatomical_Images_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on selection change in DatasetListbox.
function DatasetListbox_Callback(hObject, eventdata, handles)
% hObject    handle to DatasetListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns DatasetListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from DatasetListbox
global spanUD spanLR
spanUD = 0;
spanLR = 0;
updateImages(handles);


% --- Executes during object creation, after setting all properties.
function DatasetListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DatasetListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ZoomButton.
function ZoomButton_Callback(hObject, eventdata, handles)
% hObject    handle to ZoomButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global zoomed spanUD spanLR
disp('Draw an ROI to zoom in');
rect = drawrectangle;
positions = round(rect.Position);
spanUD = spanUD(1)+positions(2):(spanUD(1)+positions(2)+positions(4));
spanLR = spanLR(1)+positions(1):(spanLR(1)+positions(1)+positions(3));

zoomed = 1;
updateImages(handles,spanUD,spanLR);


% --- Executes on button press in UnzoomButton.
function UnzoomButton_Callback(hObject, eventdata, handles)
% hObject    handle to UnzoomButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global zoomed
zoomed = 0;
updateImages(handles);


% --- Executes on slider movement.
function ImageSlider_Callback(hObject, eventdata, handles)
% hObject    handle to ImageSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
updateImages(handles);


% --- Executes during object creation, after setting all properties.
function ImageSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImageSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function updateImages(varargin)

handles = varargin{1};
datasetNum = get(handles.DatasetListbox,'Value');
dataset = handles.anatDatasets;
imageSet = dataset(datasetNum).Data;

dim3size = size(imageSet,3);
steps = [1/(dim3size-1) 10/(dim3size-1)];
set(handles.ImageSlider, 'SliderStep', steps);
sliceNum = 1+round( get(handles.ImageSlider,'Value').*(dim3size-1) );
anatSlice = imageSet(:,:,sliceNum);

if nargin>2
    spanUD = varargin{2};
    spanLR = varargin{3};
    anatSlice = anatSlice(spanUD,spanLR);
end 

axes(handles.ImagePlot); %show image
imshow(anatSlice,[]);
set(gca,'XTickLabel','') %remove tick marks
set(gca,'YTickLabel','')
daspect([1 1 1]);
drawnow;




