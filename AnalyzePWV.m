function varargout = AnalyzePWV(varargin)
% ANALYZEPWV MATLAB code for AnalyzePWV.fig
%      ANALYZEPWV, by itself, creates a new ANALYZEPWV or raises the existing
%      singleton*.
%
%      H = ANALYZEPWV returns the handle to a new ANALYZEPWV or the handle to
%      the existing singleton*.
%
%      ANALYZEPWV('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ANALYZEPWV.M with the given input arguments.
%
%      ANALYZEPWV('Property','Value',...) creates a new ANALYZEPWV or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AnalyzePWV_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AnalyzePWV_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AnalyzePWV

% Last Modified by GUIDE v2.5 09-Oct-2019 21:05:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AnalyzePWV_OpeningFcn, ...
                   'gui_OutputFcn',  @AnalyzePWV_OutputFcn, ...
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


% --- Executes just before AnalyzePWV is made visible.
function AnalyzePWV_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AnalyzePWV (see VARARGIN)

% Choose default command line output for AnalyzePWV
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AnalyzePWV wait for user response (see UIRESUME)
% uiwait(handles.figure1);
handles.pcDatasets = varargin{1};
handles.magDatasets = varargin{2};
handles.pcDatasets(1).ROI = [];
handles.magDatasets(1).ROI = [];
set(handles.PlanePopup,'String',{handles.pcDatasets.Names});
set(handles.DatasetPopup,'String',fieldnames(handles.pcDatasets(1).Data));
guidata(hObject, handles);
updateImages(handles);
if numel(handles.magDatasets)==0
    set(handles.MAGRadio,'Enable','off')
end 


% --- Outputs from this function are returned to the command line.
function varargout = AnalyzePWV_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


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
updateImages(handles);


% --- Executes on button press in UnzoomButton.
function UnzoomButton_Callback(hObject, eventdata, handles)
% hObject    handle to UnzoomButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global zoomed spanUD spanLR
zoomed = 0; spanUD = 0; spanLR = 0;
updateImages(handles);


% --- Executes on button press in DrawROIbutton.
function DrawROIbutton_Callback(hObject, eventdata, handles)
% hObject    handle to DrawROIbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

planeNum = get(handles.PlanePopup,'Value');
if isempty(handles.pcDatasets(planeNum).ROI)
    mydlg = warndlg('Press enter when the ROI is set');
    waitfor(mydlg);
    circle = drawcircle('FaceAlpha',0.1,'Color','g','LineWidth',1,'Deletable',0);
    
while true
    w = waitforbuttonpress; 
    switch w 
        case 1 % (keyboard press) 
          key = get(gcf,'currentcharacter'); 
              switch key
                  case 27 % 27 is the escape key
                      disp('User pressed the escape key. Deleting ROI.')
                      DeleteROIbutton_Callback(hObject, eventdata, handles)
                      break % break out of the while loop
                  case 13 % 13 is the return key 
                      disp('ROI selected');
                      circle.InteractionsAllowed = 'none';
                      break
                  otherwise 
                      % Wait for a different command. 
              end
   end
end

    if get(handles.PCRadio,'Value')
        handles.pcDatasets(planeNum).ROI = circle;
    else 
        handles.magDatasets(planeNum).ROI = circle;
    end 
    guidata(hObject,handles);
    updateImages(handles);
    set(handles.DrawROIbutton,'Enable','off');
    set(handles.PlanePopup,'Enable','off');
else
    fprintf('An ROI has already been placed!\n');
end 


% --- Executes on button press in DeleteROIbutton.
function DeleteROIbutton_Callback(hObject, eventdata, handles)
% hObject    handle to DeleteROIbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hold off
planeNum = get(handles.PlanePopup,'Value');
if get(handles.PCRadio,'Value')
    handles.pcDatasets(planeNum).ROI = [];
else 
    handles.magDatasets(planeNum).ROI = [];
end 
guidata(hObject,handles);
updateImages(handles);
set(handles.DrawROIbutton,'Enable','on');
set(handles.PlanePopup,'Enable','on');

% --- Executes on button press in AnalyzeROIbutton.
function AnalyzeROIbutton_Callback(hObject, eventdata, handles)
% hObject    handle to AnalyzeROIbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in PlanePopup.
function PlanePopup_Callback(hObject, eventdata, handles)
% hObject    handle to PlanePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PlanePopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PlanePopup
global zoomed spanUD spanLR
zoomed = 0; spanUD = 0; spanLR = 0;
planeNum = get(handles.PlanePopup,'Value');
if get(handles.PCRadio,'Value')
    set(handles.DatasetPopup,'String',fieldnames(handles.pcDatasets(planeNum).Data));
end 
set(handles.DatasetPopup,'Value',1);
updateImages(handles);

% --- Executes during object creation, after setting all properties.
function PlanePopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PlanePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in DatasetPopup.
function DatasetPopup_Callback(hObject, eventdata, handles)
% hObject    handle to DatasetPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns DatasetPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from DatasetPopup
updateImages(handles);

% --- Executes during object creation, after setting all properties.
function DatasetPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DatasetPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function PlaneSlider_Callback(hObject, eventdata, handles)
% hObject    handle to PlaneSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
updateImages(handles);

% --- Executes during object creation, after setting all properties.
function PlaneSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PlaneSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in PCRadio.
function PCRadio_Callback(hObject, eventdata, handles)
% hObject    handle to PCRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PCRadio
set(handles.PlanePopup,'String',{handles.pcDatasets.Names});
set(handles.DatasetPopup,'Enable','on');
set(handles.DatasetPopup,'String',fieldnames(handles.pcDatasets(1).Data));
set(handles.MAGRadio,'Value',0);
set(handles.PlanePopup,'Value',1);
set(handles.DatasetPopup,'Value',1);
updateImages(handles);


% --- Executes on button press in MAGRadio.
function MAGRadio_Callback(hObject, eventdata, handles)
% hObject    handle to MAGRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of MAGRadio
set(handles.PlanePopup,'String',{handles.magDatasets.Names});
set(handles.DatasetPopup,'Enable','off');
set(handles.PCRadio,'Value',0);
set(handles.PlanePopup,'Value',1);
set(handles.DatasetPopup,'Value',1);
set(handles.DatasetPopup,'String',' ');
updateImages(handles);


% --- Executes during object creation, after setting all properties.
function PlanePlot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PlanePlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate PlanePlot


function varargout = updateImages(varargin)

global spanUD spanLR zoomed 
handles = varargin{1};
planeNum = get(handles.PlanePopup,'Value');
datasetNum = get(handles.DatasetPopup,'Value');
isPC = get(handles.PCRadio,'Value');

if isPC
    dataset = handles.pcDatasets;
else
    dataset = handles.magDatasets;
end 

if isstruct(dataset(planeNum).Data)
    imageSet = struct2cell(dataset(planeNum).Data);
else
    imageSet = dataset(planeNum).Data;
end 

if iscell(imageSet)
    images = imageSet(datasetNum);
    images = cell2mat(images);
else 
    images = imageSet;
end 

if ndims(images)<3
    steps = [1 1];
    set(handles.PlaneSlider, 'SliderStep', steps);
    slice = images;
else
    dim3size = size(images,3);
    steps = [1/(dim3size-1) 10/(dim3size-1)];
    set(handles.PlaneSlider, 'SliderStep', steps);
    sliceNum = 1+round( get(handles.PlaneSlider,'Value').*(dim3size-1) );
    slice = images(:,:,sliceNum);
end 

if zoomed
    slice = slice(spanUD,spanLR);
end 

varargout{1} = slice;   
axes(handles.PlanePlot);

if ~isempty(handles.pcDatasets(planeNum).ROI)
    hold on
    imshow(slice,[]);
else 
    imshow(slice,[])
end 

