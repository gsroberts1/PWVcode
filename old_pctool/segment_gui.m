function varargout = segment_gui(varargin)
% SEGMENT_GUI M-file for segment_gui.fig
%      SEGMENT_GUI, by itself, creates a new SEGMENT_GUI or raises the existing
%      singleton*.
%
%      H = SEGMENT_GUI returns the handle to a new SEGMENT_GUI or the handle to
%      the existing singleton*.
%
%      SEGMENT_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEGMENT_GUI.M with the given input arguments.
%
%      SEGMENT_GUI('Property','Value',...) creates a new SEGMENT_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before segment_gui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to segment_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help segment_gui

% Last Modified by GUIDE v2.5 15-Jun-2010 11:08:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @segment_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @segment_gui_OutputFcn, ...
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


% --- Executes just before segment_gui is made visible.
function segment_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to segment_gui (see VARARGIN)

% Choose default command line output for segment_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes segment_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = segment_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function cdname_Callback(hObject, eventdata, handles)
% hObject    handle to cdname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cdname as text
%        str2double(get(hObject,'String')) returns contents of cdname as a double


% --- Executes during object creation, after setting all properties.
function cdname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cdname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function magname_Callback(hObject, eventdata, handles)
% hObject    handle to magname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of magname as text
%        str2double(get(hObject,'String')) returns contents of magname as a double


% --- Executes during object creation, after setting all properties.
function magname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to magname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function load_menu_Callback(hObject, eventdata, handles)
% hObject    handle to load_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function load_mag_Callback(hObject, eventdata, handles)
% hObject    handle to load_mag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)





function magnam_Callback(hObject, eventdata, handles)
% hObject    handle to magname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of magname as text
%        str2double(get(hObject,'String')) returns contents of magname as a double


% --- Executes during object creation, after setting all properties.
function magnam_CreateFcn(hObject, eventdata, handles)
% hObject    handle to magname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rcxres_Callback(hObject, eventdata, handles)
% hObject    handle to rcxres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rcxres as text
%        str2double(get(hObject,'String')) returns contents of rcxres as a double


% --- Executes during object creation, after setting all properties.
function rcxres_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rcxres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rcyres_Callback(hObject, eventdata, handles)
% hObject    handle to rcyres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rcyres as text
%        str2double(get(hObject,'String')) returns contents of rcyres as a double


% --- Executes during object creation, after setting all properties.
function rcyres_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rcyres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rczres_Callback(hObject, eventdata, handles)
% hObject    handle to rczres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rczres as text
%        str2double(get(hObject,'String')) returns contents of rczres as a double


% --- Executes during object creation, after setting all properties.
function rczres_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rczres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xfov_Callback(hObject, eventdata, handles)
% hObject    handle to xfov (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xfov as text
%        str2double(get(hObject,'String')) returns contents of xfov as a double


% --- Executes during object creation, after setting all properties.
function xfov_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xfov (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function yfov_Callback(hObject, eventdata, handles)
% hObject    handle to yfov (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yfov as text
%        str2double(get(hObject,'String')) returns contents of yfov as a double


% --- Executes during object creation, after setting all properties.
function yfov_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yfov (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function zfov_Callback(hObject, eventdata, handles)
% hObject    handle to zfov (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zfov as text
%        str2double(get(hObject,'String')) returns contents of zfov as a double


% --- Executes during object creation, after setting all properties.
function zfov_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zfov (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tframes_Callback(hObject, eventdata, handles)
% hObject    handle to tframes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tframes as text
%        str2double(get(hObject,'String')) returns contents of tframes as a double


% --- Executes during object creation, after setting all properties.
function tframes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tframes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tres_Callback(hObject, eventdata, handles)
% hObject    handle to tres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tres as text
%        str2double(get(hObject,'String')) returns contents of tres as a double


% --- Executes during object creation, after setting all properties.
function tres_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in load_button.
function load_button_Callback(hObject, eventdata, handles)
% hObject    handle to load_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global xfov;
global yfov;
global zfov;
global rcxres;
global rcyres;
global rczres;
global MAG;
global CD;
global delX;
global delY;
global delZ;
global tres;
global tframes;

xfov = str2double(get(handles.xfov,'String'));
yfov = str2double(get(handles.yfov,'String'));
zfov = str2double(get(handles.zfov,'String'));

rcxres = str2double(get(handles.rcxres,'String'));
rcyres = str2double(get(handles.rcyres,'String'));
rczres = str2double(get(handles.rczres,'String'));

tstr = get(handles.rcxres,'string');
set(handles.stopX,'string',tstr);
tstr = get(handles.rcyres,'string');
set(handles.stopY,'string',tstr);
tstr = get(handles.rczres,'string');
set(handles.stopZ,'string',tstr);

tres    = str2double(get(handles.tres,'String'));
tframes = str2double(get(handles.tframes,'String'));

delX = xfov / rcxres;
delY = yfov / rcyres;
delZ = zfov / rczres;

%%Matlab won't take runtime string or number of bits
cd_type = get(handles.cd_data_type,'Value');
if cd_type == 1 %short16 little endian
CD =single( reshape(fread(fopen(get(handles.cdname,'String'),'r'),rcxres*rcyres*rczres,'short'),[rcxres rcyres rczres]) );
elseif cd_type == 2 %float32 little endian
CD =single( reshape(fread(fopen(get(handles.cdname,'String'),'r'),rcxres*rcyres*rczres,'float'),[rcxres rcyres rczres]) );
elseif cd_type == 3 %short little endian 
CD =single( reshape(fread(fopen(get(handles.cdname,'String'),'rb'),rcxres*rcyres*rczres,'short','b'),[rcxres rcyres rczres]) );
elseif cd_type == 4 %short little endian 
CD =single( reshape(fread(fopen(get(handles.cdname,'String'),'rb'),rcxres*rcyres*rczres,'float','b'),[rcxres rcyres rczres]) );
end

mag_type = get(handles.mag_data_type,'Value');
if mag_type == 1 %short16 little endian
MAG =single( reshape(fread(fopen(get(handles.magname,'String'),'r'),rcxres*rcyres*rczres,'short'),[rcxres rcyres rczres]) );
elseif mag_type == 2 %float32 little endian
MAG =single( reshape(fread(fopen(get(handles.magname,'String'),'r'),rcxres*rcyres*rczres,'float'),[rcxres rcyres rczres]) );
elseif mag_type == 3 %short little endian 
MAG =single( reshape(fread(fopen(get(handles.magname,'String'),'rb'),rcxres*rcyres*rczres,'short','b'),[rcxres rcyres rczres]) );
elseif mag_type == 4 %short little endian 
MAG =single( reshape(fread(fopen(get(handles.magname,'String'),'rb'),rcxres*rcyres*rczres,'float','b'),[rcxres rcyres rczres]) );
end

CD = CD / max(CD(:));
MAG= MAG/ max(MAG(:));

%%Set mips
axes(handles.mipx)
imagesc(reshape(max(CD,[],1),[rcyres rczres]),[0 1]);
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
colormap('gray');
daspect([1 1 1]);

axes(handles.mipy)
imagesc(reshape(max(CD,[],2),[rcxres rczres]),[0 1]);
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
colormap('gray');
daspect([1 1 1]);

axes(handles.mipz)
imagesc(reshape(max(CD,[],3),[rcxres rcyres]),[0 1]);
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
colormap('gray');
daspect([1 1 1]);



function mask_alpha_Callback(hObject, eventdata, handles)
% hObject    handle to mask_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mask_alpha as text
%        str2double(get(hObject,'String')) returns contents of mask_alpha as a double


% --- Executes during object creation, after setting all properties.
function mask_alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mask_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mask_beta_Callback(hObject, eventdata, handles)
% hObject    handle to mask_beta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mask_beta as text
%        str2double(get(hObject,'String')) returns contents of mask_beta as a double


% --- Executes during object creation, after setting all properties.
function mask_beta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mask_beta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in mask_method.
function mask_method_Callback(hObject, eventdata, handles)
% hObject    handle to mask_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns mask_method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from mask_method

%%%%KMJ Segmentation Method %%%%
mask_method = get(handles.mask_method,'Value');
%%%Decoding
%% 1 = mask on CD
%% 2 = mask on Max
%% 3 = mixture model (not there yet)
%% 4 = load from mimics
%% 5 = load from mimics CE-MRA

if mask_method < 4 
    
   set(handles.mask_alpha_name,'Visible','on');
   set(handles.mask_alpha,'Visible','on');
   set(handles.mask_beta_name,'Visible','on');
   set(handles.mask_beta,'Visible','on');
   set(handles.mask_iterations_name,'Visible','on');
   set(handles.mask_iterations,'Visible','on');
        
   set(handles.mimics_name,'Visible','off');
   set(handles.mimics_push,'Visible','off');
   set(handles.mimics_list,'Visible','off');
else
    set(handles.mask_alpha_name,'Visible','off');
   set(handles.mask_alpha,'Visible','off');
   set(handles.mask_beta_name,'Visible','off');
   set(handles.mask_beta,'Visible','off');
   set(handles.mask_iterations_name,'Visible','off');
   set(handles.mask_iterations,'Visible','off');
   
    set(handles.mimics_name,'Visible','on');
   set(handles.mimics_push,'Visible','on');
   set(handles.mimics_list,'Visible','on');
end    


% --- Executes during object creation, after setting all properties.
function mask_method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mask_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mask_iterations_Callback(hObject, eventdata, handles)
% hObject    handle to mask_iterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mask_iterations as text
%        str2double(get(hObject,'String')) returns contents of mask_iterations as a double


% --- Executes during object creation, after setting all properties.
function mask_iterations_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mask_iterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in mask_button.
function mask_button_Callback(hObject, eventdata, handles)
% hObject    handle to mask_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global m_xstart;
global m_xstop;
global m_ystart;
global m_ystop;
global m_zstart;
global m_zstop;
global m_alpha;
global m_beta;

global m_xlength;
global m_ylength;
global m_zlength;

global xfov;
global yfov;
global zfov;
global rcxres;
global rcyres;
global rczres;
global MAG;
global CD;
global MASK;

global U_pcvipr;
global O_pcvipr;
    


%%%%KMJ Segmentation Method %%%%
mask_method = get(handles.mask_method,'Value');
%%%Decoding
%% 1 = mask on CD
%% 2 = mask on Max
%% 3 = mixture model (not there yet)
%% 4 = load from mimics
%% 5 = load from mimics MRA

if mask_method < 4

    m_alpha = str2double(get(handles.mask_alpha,'String'));
    m_beta  = str2double(get(handles.mask_beta,'String'));
    m_iter  = str2double(get(handles.mask_iterations,'String'));

    %%%%%SIMPLE THRESHOLD TO BEGING WITH %%%%%%%%%%
    if( mask_method == 1)
        MASK = single( CD > m_alpha);
    elseif( mask_method ==2)
        MASK = single( MAG > m_alpha);
    end

    %%%%%%ANTI ISLAND FILTER %%%%%%%%%%%%%%%%%%%%%%
    newpixels = 10;
    iter=0;
    level = 27*m_beta;

    while( (0 <newpixels) && (iter < m_iter) )
        newpixels=0;
        iter = iter +1;
        for x=m_xstart:m_xstop
            for y=m_ystart:m_ystop
                for z=m_zstart:m_zstop

                    if(MASK(x,y,z) > 0)
                        mag_area = sum(CD( x-1:x+1,y-1:y+1,z-1:z+1),'native');
                        mask_area= (MASK( (x-1):(x+1),(y-1):(y+1),(z-1):(z+1)));
                        mask_sum = sum(mask_area(:));

                        if(mask_sum < level)
                            MASK(x,y,z)=0;
                            newpixels = newpixels +1;
                        end
                    end
                end
            end
        end
    end

else
    %%%Have to convert to scanner coordinates
    MASK = zeros( size(CD));
        
    mimics_names = get(handles.mimics_list,'String');

    for num= 1:length(mimics_names)
    %%%Load up the mimics file
    [xM yM zM intensity] = textread(mimics_names{num},'%f,%f,%f,%f');
        
    U_pcvipr
    O_pcvipr
    
    %%Currently Mimics Apearrs to Convert X/Y from [0 FOV]. Only works for
    %%axial
    O_pcvipr(1)= 0.0;
    O_pcvipr(2)= 0.0;
    
    Inv_pcvipr = inv(U_pcvipr);
      
    for pos = 1:length(xM)
        temp = Inv_pcvipr*( [xM(pos); yM(pos); zM(pos)] - O_pcvipr);
        xM(pos)=temp(1);
        yM(pos)=temp(2);
        zM(pos)=temp(3);
    end
    
% ix 1.0937500000
% iy -0.0000000000
% iz -0.0000000000
% jx 0.0000000000
% jy 1.0937500000
% jz -0.0000000000
% kx -0.0000000000
% ky 0.0000000000
% kz -1.0937500000
% sx -140.0000000000
% sy -122.0470962524
% sz 166.7861022949
    
    idx = sub2ind(size(CD),round(xM),round(yM),round(zM));
    MASK(idx)=num;
    end
end
    


update_images(handles);

% --- Executes on selection change in visual_method.
function visual_method_Callback(hObject, eventdata, handles)
% hObject    handle to visual_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns visual_method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from visual_method

update_images(handles)

% --- Executes during object creation, after setting all properties.
function visual_method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to visual_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function mipz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mipz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate mipz
imagesc(zeros(10,10),[0 1]);
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
colormap('gray');



% --- Executes during object creation, after setting all properties.
function mipx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mipx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate mipx
imagesc(zeros(10,10),[0 1]);
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
colormap('gray');




% --- Executes during object creation, after setting all properties.
function mipy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mipy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate mipy

imagesc(zeros(10,10),[0 1]);
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
colormap('gray');

function startX_Callback(hObject, eventdata, handles)
% hObject    handle to startX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of startX as text
%        str2double(get(hObject,'String')) returns contents of startX as a double
update_images(handles);

% --- Executes during object creation, after setting all properties.
function startX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to startX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function stopX_Callback(hObject, eventdata, handles)
% hObject    handle to stopX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stopX as text
%        str2double(get(hObject,'String')) returns contents of stopX as a double
update_images(handles);

% --- Executes during object creation, after setting all properties.
function stopX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stopX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function startY_Callback(hObject, eventdata, handles)
% hObject    handle to startY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of startY as text
%        str2double(get(hObject,'String')) returns contents of startY as a double
update_images(handles);

% --- Executes during object creation, after setting all properties.
function startY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to startY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function stopY_Callback(hObject, eventdata, handles)
% hObject    handle to stopY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stopY as text
%        str2double(get(hObject,'String')) returns contents of stopY as a double

update_images(handles);

% --- Executes during object creation, after setting all properties.
function stopY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stopY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double


% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function stopZ_Callback(hObject, eventdata, handles)
% hObject    handle to stopZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stopZ as text
%        str2double(get(hObject,'String')) returns contents of stopZ as a double
update_images(handles);


% --- Executes during object creation, after setting all properties.
function stopZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stopZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function startZ_Callback(hObject, eventdata, handles)
% hObject    handle to startZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of startZ as text
%        str2double(get(hObject,'String')) returns contents of startZ as a double
update_images(handles);

% --- Executes during object creation, after setting all properties.
function startZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to startZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function update_images(handles)

global m_xstart;
global m_xstop;
global m_ystart;
global m_ystop;
global m_zstart;
global m_zstop;

global m_xlength;
global m_ylength;
global m_zlength;

global xfov;
global yfov;
global zfov;
global rcxres;
global rcyres;
global rczres;
global MAG;
global CD;
global MASK;

%%%COPY VALUES%%%%%%
gridlines =(get(handles.gridlines,'Value'));
m_xstart = str2double(get(handles.startX,'String'));
m_xstop  = str2double(get(handles.stopX,'String'));
m_ystart = str2double(get(handles.startY,'String'));
m_ystop  = str2double(get(handles.stopY,'String'));
m_zstart = str2double(get(handles.startZ,'String'));
m_zstop  = str2double(get(handles.stopZ,'String'));
rcxres = str2double(get(handles.rcxres,'String'));
rcyres = str2double(get(handles.rcyres,'String'));
rczres = str2double(get(handles.rczres,'String'));

gamma = get(handles.gamma_slider,'Value');

visual_method = get(handles.visual_method,'Value');
%%%VISUAL METHOD DECODING
%   1 = cd mip
%   2 = mag mip
%   3 = cd mip + overlay
%   4 = mag mip + overlay

m_xlength = m_xstop - m_xstart +1;
m_ylength = m_ystop - m_ystart +1;
m_zlength = m_zstop - m_zstart +1 ;

%%Set mips

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% XMIPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
blank=128*ones(rcyres+128,rczres+128);
axes(handles.mipx)
if( (visual_method ==1) || (visual_method==3) )
    im = (reshape(max(CD(m_xstart:m_xstop,m_ystart:m_ystop,m_zstart:m_zstop),[],1),[m_ylength m_zlength]));
elseif( (visual_method ==2) || (visual_method==4) )
    im = (reshape(max(MAG(m_xstart:m_xstop,m_ystart:m_ystop,m_zstart:m_zstop),[],1),[m_ylength m_zlength]));
end

im = im * 256 / max(im(:));

if( visual_method > 2)
    %mask_mip =(reshape(max(MASK(m_xstart:m_xstop,m_ystart:m_ystop,m_zstart:m_zstop),[],1),[m_ylength m_zlength]));
    for posy = m_ystart:m_ystop
        for posz = m_zstart:m_zstop 
           idx = find( squeeze(MASK(:,posy,posz)),1,'first');
           if numel(idx) ~= 0
               im(posy-m_ystart+1,posz-m_zstart+1) = 256 + MASK(idx,posy, posz); 
           end
        end
    end
    blank(m_ystart:m_ystop,m_zstart:m_zstop)=im;
    masks = max(im(:)) - 256
    map =[jet(masks); gray(256)];

    imagesc(blank,[0 (257+masks)]);
    colormap(map);
else
    blank(m_ystart:m_ystop,m_zstart:m_zstop)=im.^gamma;
    imagesc(blank,[0 256]);
    colormap('gray');
end
ylim([ m_ystart (m_ystart+max([m_ylength m_zlength])-1)]);
xlim([ m_zstart (m_zstart+max([m_ylength m_zlength])-1)]);
if(gridlines ==1 )
    grid on;
    set(gca,'XColor','r');
    set(gca,'YColor','r');
    xlabel('Z-Pos');
    ylabel('Y-Pos');
else
    grid off
    set(gca,'XTickLabel','')
    set(gca,'YTickLabel','')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% YMIPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(handles.mipy)
blank=128*ones(rcxres+128,rczres+128);
if( (visual_method ==1) || (visual_method==3) )
    im =(reshape(max(CD(m_xstart:m_xstop,m_ystart:m_ystop,m_zstart:m_zstop),[],2),[m_xlength m_zlength]));
elseif( (visual_method ==2) || (visual_method==4) )
    im =(reshape(max(MAG(m_xstart:m_xstop,m_ystart:m_ystop,m_zstart:m_zstop),[],2),[m_xlength m_zlength]));
end

im = im * 256 / max(im(:));

if( visual_method > 2)
   %mask_mip =(reshape(max(MASK(m_xstart:m_xstop,m_ystart:m_ystop,m_zstart:m_zstop),[],1),[m_ylength m_zlength]));
   for posx = m_xstart:m_xstop
        for posz = m_zstart:m_zstop 
           idx = find( squeeze(MASK(posx,:,posz)),1,'first');
           if numel(idx) ~= 0
               im(posx-m_xstart+1,posz-m_zstart+1) = 256 + MASK(posx,idx, posz); 
           end
        end
    end
    blank(m_xstart:m_xstop,m_zstart:m_zstop)=im;
    masks = max(im(:)) - 256
    map =[jet(masks); gray(256)];
    imagesc(blank,[0 (257+masks)]);
    colormap(map);
else
    blank(m_xstart:m_xstop,m_zstart:m_zstop)=im.^gamma;
    imagesc(blank);
    colormap('gray');
end

if(gridlines ==1 )
    grid on;
    set(gca,'XColor','r');  
    set(gca,'YColor','r');
    xlabel('Z-Pos');
    ylabel('X-Pos');
else
    grid off
    set(gca,'XTickLabel','')
    set(gca,'YTickLabel','')
end
xlim([m_zstart (m_zstart -1+ max([m_xlength m_zlength]))]);
ylim([m_xstart (m_xstart -1+ max([m_xlength m_zlength]))]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ZMIPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(handles.mipz)
blank=128*ones(rcxres+128,rcyres+128);

if( (visual_method ==1) || (visual_method==3) )
    im = (reshape(max(CD(m_xstart:m_xstop,m_ystart:m_ystop,m_zstart:m_zstop),[],3),[m_xlength m_ylength]));
elseif( (visual_method ==2) || (visual_method==4) )
    im = (reshape(max(MAG(m_xstart:m_xstop,m_ystart:m_ystop,m_zstart:m_zstop),[],3),[m_xlength m_ylength]));
end
im = im * 256 / max(im(:));

if( visual_method > 2)
   for posx = m_xstart:m_xstop
        for posy = m_ystart:m_ystop 
           idx = find( squeeze(MASK(posx,posy,:)),1,'first');
           if numel(idx) ~= 0
               im(posx-m_xstart+1,posy-m_ystart+1) = 256 + MASK(posx,posy,idx); 
           end
        end
    end
    %blank(m_ystart:m_ystop,m_zstart:m_zstop)=im.^gamma + mask_mip;
    blank(m_xstart:m_xstop,m_ystart:m_ystop)= im;
    masks = max(im(:)) - 256
    map =[jet(masks); gray(256)];
    imagesc(blank,[0 (257+masks)]);
    colormap(map);
else
    blank(m_xstart:m_xstop,m_ystart:m_ystop)=im.^gamma;
    imagesc(blank);
    colormap('gray');
end

if(gridlines ==1 )
    grid on;
    set(gca,'XColor','r');
    set(gca,'YColor','r');
    xlabel('Y-Pos');
    ylabel('X-Pos');
else
    grid off
    set(gca,'XTickLabel','')
    set(gca,'YTickLabel','')
end

xlim([(m_ystart) (m_ystart-1 +max([m_ylength m_xlength]))]);
ylim([(m_xstart) (m_xstart-1 +max([m_ylength m_xlength]))]);



% --- Executes on button press in gridlines.
function gridlines_Callback(hObject, eventdata, handles)
% hObject    handle to gridlines (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gridlines
update_images(handles)



% --- Executes on slider movement.
function gamma_slider_Callback(hObject, eventdata, handles)
% hObject    handle to gamma_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

update_images(handles)


% --- Executes during object creation, after setting all properties.
function gamma_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gamma_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end




% --- Executes on selection change in cd_data_type.
function cd_data_type_Callback(hObject, eventdata, handles)
% hObject    handle to cd_data_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns cd_data_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cd_data_type


% --- Executes during object creation, after setting all properties.
function cd_data_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cd_data_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in mag_data_type.
function mag_data_type_Callback(hObject, eventdata, handles)
% hObject    handle to mag_data_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns mag_data_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from mag_data_type


% --- Executes during object creation, after setting all properties.
function mag_data_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mag_data_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in parameter_read.
function parameter_read_Callback(hObject, eventdata, handles)
% hObject    handle to parameter_read (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global U_pcvipr;
global O_pcvipr;

%%%Reads external files for parameters
[parameter value]=textread('pcvipr_header.txt', '%s %s');  

%%%FOVS
idx = find(strcmp('fovx',parameter));
xfov = str2num( value{idx} );
set(handles.xfov,'string',num2str(xfov));

idx = find(strcmp('fovy',parameter));
yfov = str2num( value{idx} );
set(handles.yfov,'string',num2str(yfov));

idx = find(strcmp('fovz',parameter));
zfov = str2num( value{idx} );
set(handles.zfov,'string',num2str(zfov));

%%%Matrix
idx = find(strcmp('matrixx',parameter));
rcxres = str2num( value{idx} );
set(handles.rcxres,'string',num2str(rcxres));

idx = find(strcmp('matrixy',parameter));
rcyres = str2num( value{idx} );
set(handles.rcyres,'string',num2str(rcyres));

idx = find(strcmp('matrixz',parameter));
rczres = str2num( value{idx} );
set(handles.rczres,'string',num2str(rczres));

%%%Time Stuff
idx = find(strcmp('frames',parameter));
frames = str2num( value{idx} );
if(frames==-1)
    frames = 3;
end
set(handles.tframes,'string',num2str(frames));

idx = find(strcmp('timeres',parameter));
tres = str2num( value{idx} )/1000;
set(handles.tres,'string',num2str(tres));

idx = find(strcmp('version',parameter));
version = str2num( value{idx} );

if version > 1
    idx = find(strcmp('ix',parameter));
    ix = str2num( value{idx} );
    idx = find(strcmp('iy',parameter));
    iy = str2num( value{idx} );
    idx = find(strcmp('iz',parameter));
    iz = str2num( value{idx} );

    idx = find(strcmp('jx',parameter));
    jx = str2num( value{idx} );
    idx = find(strcmp('jy',parameter));
    jy = str2num( value{idx} );
    idx = find(strcmp('jz',parameter));
    jz = str2num( value{idx} );

    idx = find(strcmp('kx',parameter));
    kx = str2num( value{idx} );
    idx = find(strcmp('ky',parameter));
    ky = str2num( value{idx} );
    idx = find(strcmp('kz',parameter));
    kz = str2num( value{idx} );
    
    idx = find(strcmp('sx',parameter));
    sx = str2num( value{idx} );
    idx = find(strcmp('sy',parameter));
    sy = str2num( value{idx} );
    idx = find(strcmp('sz',parameter));
    sz = str2num( value{idx} );

    %%PC VIPR Orientation
    U_pcvipr = zeros(3,3);
    U_pcvipr(1,:) = [ ix  jx kx];
    U_pcvipr(2,:) = [ iy  jy ky];
    U_pcvipr(3,:) = [ iz  jz kz];
    O_pcvipr      = [ sx; sy; sz];
else
   disp('ERROR:Please Recon with Latest Recon to Get header!!');
end

% --- Executes on button press in cine_anatomy_checkbox.
function cine_anatomy_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to cine_anatomy_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cine_anatomy_checkbox
global cine_anatomy_flag;
cine_anatomy_flag = get(hObject,'Value');



function mimics_name_Callback(hObject, eventdata, handles)
% hObject    handle to mimics_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mimics_name as text
%        str2double(get(hObject,'String')) returns contents of mimics_name as a double


% --- Executes during object creation, after setting all properties.
function mimics_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mimics_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in mimics_push.
function mimics_push_Callback(hObject, eventdata, handles)
% hObject    handle to mimics_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[name path] = uigetfile( ...
{  '*.txt','Text (*.txt)'; ...
   '*.*',  'All Files (*.*)'}, ...
   'Select All Mimics Files', ...
   'MultiSelect', 'on');
if iscell(name)
    mimics_files = size(name,2)
    for ii = 1:mimics_files 
        mimics_names{ii} = [path name{ii}]
    end
else
    mimics_names{1} = [path name]
end
set(handles.mimics_list,'String',mimics_names)

% --- Executes on selection change in mimics_list.
function mimics_list_Callback(hObject, eventdata, handles)
% hObject    handle to mimics_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns mimics_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from mimics_list


% --- Executes during object creation, after setting all properties.
function mimics_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mimics_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in preesure_checkbox.
function pressure_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to preesure_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of preesure_checkbox
global pressure_flag;
pressure_flag = (get(handles.pressure_checkbox,'Value'));

% --- Executes on button press in wss_checkbox.
function wss_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to wss_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of wss_checkbox
global wss_flag;
wss_flag = (get(handles.wss_checkbox,'Value'));


% --- Executes on button press in flow_checkbox.
function flow_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to flow_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of flow_checkbox
global flow_flag;
flow_flag = (get(handles.flow_checkbox,'Value'));



% --- Executes on button press in ensight_checkbox.
function ensight_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to ensight_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ensight_checkbox
global ensight_flag;
ensight_flag = get(hObject,'Value');



% --- Executes on button press in frame_by_frame_checkbox.
function frame_by_frame_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to frame_by_frame_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of frame_by_frame_checkbox
global frame_by_frame;
frame_by_frame = get(hObject,'Value');

% --- Executes on button press in pwv.
function pwv_Callback(hObject, eventdata, handles)
% hObject    handle to pwv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pwv
global pwv_flag; pwv_flag = get(hObject,'Value'); %alw

% --- Executes on button press in load_mra_button.
function load_mra_button_Callback(hObject, eventdata, handles)
% hObject    handle to load_mra_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%This loads a CE-MRA Dataset
path = uigetdir('','Select MRA Dicom Folder');

global mra_image;
global mra_info;
global mra_meta_data;

h = waitbar(0,'Reading DICOM Files...','WindowStyle','modal');

dicom_files = dir(path);

slice_pos = 1;
N = numel(dicom_files);
info_exists =0;
for file_pos = 1 : N
    if dicom_files(file_pos,1).isdir == 0
        if info_exists == 0
            mra_info = dicominfo([path,'\',dicom_files(file_pos,1).name]);            
            mra_image =zeros(mra_info.Width,mra_info.Height,N-2);
            info_exists =1;
        end
        mra_image(:,:,slice_pos) = dicomread([path,'\',dicom_files(file_pos,1).name]);
        slice_pos = slice_pos + 1;
        waitbar(file_pos/N*0.3,h);
    end
end


%%%Convert 
global MRA;
global CD;
MRA = zeros(size(CD));

global rcxres;
global rcyres;
global rczres;

global O_pcvipr
global U_pcvipr

%%%MRA Unit Vectors
OV_mra = mra_info.ImagePositionPatient;
XV_mra = mra_info.ImageOrientationPatient(1:3)*mra_info.ReconstructionDiameter/mra_info.Width;
YV_mra = mra_info.ImageOrientationPatient(4:6)*mra_info.ReconstructionDiameter/mra_info.Height;
ZV_mra = -cross(XV_mra,YV_mra).*mra_info.SliceThickness;
Orientation = mra_info.PatientPosition;

U_mra(1,:) = [ XV_mra(1) YV_mra(1) ZV_mra(1)];
U_mra(2,:) = [ XV_mra(2) YV_mra(2) ZV_mra(2)];
U_mra(3,:) = [ XV_mra(3) YV_mra(3) ZV_mra(3)];
O_mra      = [ OV_mra(1); OV_mra(2); OV_mra(3)];
UI_mra = inv(U_mra);

U_mra
UI_mra
O_mra

U_pcvipr
O_pcvipr

[xt yt] = meshgrid(0:rcxres-1,0:rcyres-1);
for slice = 1:6:rczres
    waitbar(slice/rczres*0.7+0.3,h);
    
    %%Raw Scanner Coordinates for PC VIPR
    XX_pcvipr_raw = O_pcvipr(1) + xt*U_pcvipr(1,1) + yt*U_pcvipr(1,2) + (slice-1)*U_pcvipr(1,3);
    YY_pcvipr_raw = O_pcvipr(2) + xt*U_pcvipr(2,1) + yt*U_pcvipr(2,2) + (slice-1)*U_pcvipr(2,3);
    ZZ_pcvipr_raw = O_pcvipr(3) + xt*U_pcvipr(3,1) + yt*U_pcvipr(3,2) + (slice-1)*U_pcvipr(3,3);

    disp(['Slice:',num2str(slice),'  Range X [',num2str(min(XX_pcvipr_raw(:))),' ',num2str(max(XX_pcvipr_raw(:))),' ] Y [',num2str(min(YY_pcvipr_raw(:))),' ',num2str(max(YY_pcvipr_raw(:))),' ] Z [ ',num2str(min(ZZ_pcvipr_raw(:))),' ',num2str(max(ZZ_pcvipr_raw(:))),' ] ']);
    
    %%%Convert to MRA coordinates 
    %%% A* x-mra + b = x-scanner
    %%% A^-1*( x-scanner - b) = x-mra
    XX_pcvipr_raw = XX_pcvipr_raw - O_mra(1);
    YY_pcvipr_raw = YY_pcvipr_raw - O_mra(2);
    ZZ_pcvipr_raw = ZZ_pcvipr_raw - O_mra(3);

    XX_pcvipr_mra = UI_mra(1,1)*XX_pcvipr_raw + UI_mra(1,2)*YY_pcvipr_raw + UI_mra(1,3)*ZZ_pcvipr_raw;
    YY_pcvipr_mra = UI_mra(2,1)*XX_pcvipr_raw + UI_mra(2,2)*YY_pcvipr_raw + UI_mra(2,3)*ZZ_pcvipr_raw;
    ZZ_pcvipr_mra = UI_mra(3,1)*XX_pcvipr_raw + UI_mra(3,2)*YY_pcvipr_raw + UI_mra(3,3)*ZZ_pcvipr_raw;
    
    disp(['Slice:',num2str(slice),'  Range X [',num2str(min(XX_pcvipr_mra(:))),' ',num2str(max(XX_pcvipr_mra(:))),' ] Y [',num2str(min(YY_pcvipr_mra(:))),' ',num2str(max(YY_pcvipr_mra(:))),' ] Z [ ',num2str(min(ZZ_pcvipr_mra(:))),' ',num2str(max(ZZ_pcvipr_mra(:))),' ] ']);

    MRA(:,:,slice) = volume_interp( mra_image,XX_pcvipr_mra,YY_pcvipr_mra,ZZ_pcvipr_mra,[rcxres rcyres]);
end
close(h);

global MAG
MAG = MRA;