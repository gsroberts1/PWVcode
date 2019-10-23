function varargout = EnterZslices(varargin)
% ENTERZSLICES MATLAB code for EnterZslices.fig
%      ENTERZSLICES, by itself, creates a new ENTERZSLICES or raises the existing
%      singleton*.
%
%      H = ENTERZSLICES returns the handle to a new ENTERZSLICES or the handle to
%      the existing singleton*.
%
%      ENTERZSLICES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ENTERZSLICES.M with the given input arguments.
%
%      ENTERZSLICES('Property','Value',...) creates a new ENTERZSLICES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EnterZslices_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EnterZslices_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EnterZslices

% Last Modified by GUIDE v2.5 22-Oct-2019 13:39:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EnterZslices_OpeningFcn, ...
                   'gui_OutputFcn',  @EnterZslices_OutputFcn, ...
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

% --- Executes just before EnterZslices is made visible.
function EnterZslices_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to EnterZslices (see VARARGIN)

    % Choose default command line output for EnterZslices
    handles.output = hObject;
    handles.flow = varargin{1};
    
     switch numel(handles.flow)
        case 1
            set(handles.slice2edit,'Enable','off');
            set(handles.slice3edit,'Enable','off');
            set(handles.slice4edit,'Enable','off');
            set(handles.slice5edit,'Enable','off');
            set(handles.slice6edit,'Enable','off');
            set(handles.slice7edit,'Enable','off');
            set(handles.slice8edit,'Enable','off');
            set(handles.slice2text,'Enable','off');
            set(handles.slice3text,'Enable','off');
            set(handles.slice4text,'Enable','off');
            set(handles.slice5text,'Enable','off');
            set(handles.slice6text,'Enable','off');
            set(handles.slice7text,'Enable','off');
            set(handles.slice8text,'Enable','off');
            set(handles.slice1text,'String',handles.flow(1).Name);
        case 2
            set(handles.slice3edit,'Enable','off');
            set(handles.slice4edit,'Enable','off');
            set(handles.slice5edit,'Enable','off');
            set(handles.slice6edit,'Enable','off');
            set(handles.slice7edit,'Enable','off');
            set(handles.slice8edit,'Enable','off');
            set(handles.slice3text,'Enable','off');
            set(handles.slice4text,'Enable','off');
            set(handles.slice5text,'Enable','off');
            set(handles.slice6text,'Enable','off');
            set(handles.slice7text,'Enable','off');
            set(handles.slice8text,'Enable','off');
            set(handles.slice1text,'String',handles.flow(1).Name);
            set(handles.slice2text,'String',handles.flow(2).Name);
        case 3
            set(handles.slice4edit,'Enable','off');
            set(handles.slice5edit,'Enable','off');
            set(handles.slice6edit,'Enable','off');
            set(handles.slice7edit,'Enable','off');
            set(handles.slice8edit,'Enable','off');
            set(handles.slice4text,'Enable','off');
            set(handles.slice5text,'Enable','off');
            set(handles.slice6text,'Enable','off');
            set(handles.slice7text,'Enable','off');
            set(handles.slice8text,'Enable','off');
            set(handles.slice1text,'String',handles.flow(1).Name);
            set(handles.slice2text,'String',handles.flow(2).Name);
            set(handles.slice3text,'String',handles.flow(3).Name);
        case 4
            set(handles.slice5edit,'Enable','off');
            set(handles.slice6edit,'Enable','off');
            set(handles.slice7edit,'Enable','off');
            set(handles.slice8edit,'Enable','off');
            set(handles.slice5text,'Enable','off');
            set(handles.slice6text,'Enable','off');
            set(handles.slice7text,'Enable','off');
            set(handles.slice8text,'Enable','off');
            set(handles.slice1text,'String',handles.flow(1).Name);
            set(handles.slice2text,'String',handles.flow(2).Name);
            set(handles.slice3text,'String',handles.flow(3).Name);
            set(handles.slice4text,'String',handles.flow(4).Name);
        case 5
            set(handles.slice6edit,'Enable','off');
            set(handles.slice7edit,'Enable','off');
            set(handles.slice8edit,'Enable','off');
            set(handles.slice6text,'Enable','off');
            set(handles.slice7text,'Enable','off');
            set(handles.slice8text,'Enable','off');
            set(handles.slice1text,'String',handles.flow(1).Name);
            set(handles.slice2text,'String',handles.flow(2).Name);
            set(handles.slice3text,'String',handles.flow(3).Name);
            set(handles.slice4text,'String',handles.flow(4).Name);
            set(handles.slice5text,'String',handles.flow(5).Name);       
        case 6
            set(handles.slice7edit,'Enable','off');
            set(handles.slice8edit,'Enable','off');
            set(handles.slice7text,'Enable','off');
            set(handles.slice8text,'Enable','off');
            set(handles.slice1text,'String',handles.flow(1).Name);
            set(handles.slice2text,'String',handles.flow(2).Name);
            set(handles.slice3text,'String',handles.flow(3).Name);
            set(handles.slice4text,'String',handles.flow(4).Name);
            set(handles.slice5text,'String',handles.flow(5).Name);  
            set(handles.slice6text,'String',handles.flow(6).Name);  
        case 7 
            set(handles.slice8edit,'Enable','off');
            set(handles.slice8text,'Enable','off');
            set(handles.slice1text,'String',handles.flow(1).Name);
            set(handles.slice2text,'String',handles.flow(2).Name);
            set(handles.slice3text,'String',handles.flow(3).Name);
            set(handles.slice4text,'String',handles.flow(4).Name);
            set(handles.slice5text,'String',handles.flow(5).Name);  
            set(handles.slice6text,'String',handles.flow(6).Name);  
            set(handles.slice7text,'String',handles.flow(7).Name);          
        case 8
            set(handles.slice1text,'String',handles.flow(1).Name);
            set(handles.slice2text,'String',handles.flow(2).Name);
            set(handles.slice3text,'String',handles.flow(3).Name);
            set(handles.slice4text,'String',handles.flow(4).Name);
            set(handles.slice5text,'String',handles.flow(5).Name);  
            set(handles.slice6text,'String',handles.flow(6).Name);  
            set(handles.slice7text,'String',handles.flow(7).Name);  
            set(handles.slice8text,'String',handles.flow(8).Name);  
        otherwise
    end 
    
    % Update handles structure
    guidata(hObject, handles);

    % UIWAIT makes EnterZslices wait for user response (see UIRESUME)
    uiwait(handles.zLocsGUI);


% --- Outputs from this function are returned to the command line.
function varargout = EnterZslices_OutputFcn(hObject, eventdata, handles) 
    % varargout  cell array for returning output args (see VARARGOUT);
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure
    varargout{1} = handles.zLocs;
    delete(handles.zLocsGUI);

%%% SLICE 1
function slice1edit_Callback(hObject, eventdata, handles)

function slice1edit_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function slice1text_CreateFcn(hObject, eventdata, handles)



%%% SLICE 2
function slice2edit_Callback(hObject, eventdata, handles)

function slice2edit_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function slice2text_CreateFcn(hObject, eventdata, handles)



%%% SLICE 3
function slice3edit_Callback(hObject, eventdata, handles)

function slice3edit_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function slice3text_CreateFcn(hObject, eventdata, handles)



%%% SLICE 4
function slice4edit_Callback(hObject, eventdata, handles)

function slice4edit_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function slice4text_CreateFcn(hObject, eventdata, handles)



%%% SLICE 5
function slice5edit_Callback(hObject, eventdata, handles)

function slice5edit_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function slice5text_CreateFcn(hObject, eventdata, handles)



%%% SLICE 6
function slice6edit_Callback(hObject, eventdata, handles)

function slice6edit_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function slice6text_CreateFcn(hObject, eventdata, handles)



%%% SLICE 7
function slice7edit_Callback(hObject, eventdata, handles)

function slice7edit_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function slice7text_CreateFcn(hObject, eventdata, handles)


%%% SLICE 8
function slice8edit_Callback(hObject, eventdata, handles)

function slice8edit_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function slice8text_CreateFcn(hObject, eventdata, handles)


%%% Complete Loading
function completeLoadingButton_Callback(hObject, eventdata, handles)
    switch numel(handles.flow)
        case 1
            zLocs(1) = str2num(get(handles.slice1edit,'String'));
        case 2
            zLocs(1) = str2num(get(handles.slice1edit,'String'));
            zLocs(2) = str2num(get(handles.slice2edit,'String'));
        case 3
            zLocs(1) = str2num(get(handles.slice1edit,'String'));
            zLocs(2) = str2num(get(handles.slice2edit,'String'));
            zLocs(3) = str2num(get(handles.slice3edit,'String'));
        case 4
            zLocs(1) = str2num(get(handles.slice1edit,'String'));
            zLocs(2) = str2num(get(handles.slice2edit,'String'));
            zLocs(3) = str2num(get(handles.slice3edit,'String'));
            zLocs(4) = str2num(get(handles.slice4edit,'String'));
        case 5
            zLocs(1) = str2num(get(handles.slice1edit,'String'));
            zLocs(2) = str2num(get(handles.slice2edit,'String'));
            zLocs(3) = str2num(get(handles.slice3edit,'String'));
            zLocs(4) = str2num(get(handles.slice4edit,'String'));
            zLocs(5) = str2num(get(handles.slice5edit,'String'));      
        case 6
            zLocs(1) = str2num(get(handles.slice1edit,'String'));
            zLocs(2) = str2num(get(handles.slice2edit,'String'));
            zLocs(3) = str2num(get(handles.slice3edit,'String'));
            zLocs(4) = str2num(get(handles.slice4edit,'String'));
            zLocs(5) = str2num(get(handles.slice5edit,'String')); 
            zLocs(6) = str2num(get(handles.slice6edit,'String'));
        case 7 
            zLocs(1) = str2num(get(handles.slice1edit,'String'));
            zLocs(2) = str2num(get(handles.slice2edit,'String'));
            zLocs(3) = str2num(get(handles.slice3edit,'String'));
            zLocs(4) = str2num(get(handles.slice4edit,'String'));
            zLocs(5) = str2num(get(handles.slice5edit,'String'));
            zLocs(6) = str2num(get(handles.slice6edit,'String'));
            zLocs(7) = str2num(get(handles.slice7edit,'String'));         
        case 8
            zLocs(1) = str2num(get(handles.slice1edit,'String'));
            zLocs(2) = str2num(get(handles.slice2edit,'String'));
            zLocs(3) = str2num(get(handles.slice3edit,'String'));
            zLocs(4) = str2num(get(handles.slice4edit,'String'));
            zLocs(5) = str2num(get(handles.slice5edit,'String'));
            zLocs(6) = str2num(get(handles.slice6edit,'String'));
            zLocs(7) = str2num(get(handles.slice7edit,'String'));
            zLocs(8) = str2num(get(handles.slice8edit,'String'));
        otherwise
    end 
handles.zLocs = zLocs;
guidata(hObject,handles);
close(handles.zLocsGUI);
 
%%% Close GUI
function zLocsGUI_CloseRequestFcn(hObject, eventdata, handles)
    if isequal(get(hObject, 'waitstatus'), 'waiting')
        % The GUI is still in UIWAIT, us UIRESUME
        uiresume(hObject);
    else
        % The GUI is no longer waiting, just close it
        delete(hObject);
    end   
