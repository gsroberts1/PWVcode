function varargout = flow_gui(varargin)
% FLOW_GUI M-file for flow_gui.fig
%      FLOW_GUI, by itself, creates a new FLOW_GUI or raises the existing
%      singleton*.
%
%      H = FLOW_GUI returns the handle to a new FLOW_GUI or the handle to
%      the existing singleton*.
%
%      FLOW_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FLOW_GUI.M with the given input arguments.
%
%      FLOW_GUI('Property','Value',...) creates a new FLOW_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before flow_gui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to flow_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help flow_gui

% Last Modified by GUIDE v2.5 14-Mar-2008 16:52:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @flow_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @flow_gui_OutputFcn, ...
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


% --- Executes just before flow_gui is made visible.
function flow_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to flow_gui (see VARARGIN)

% Choose default command line output for flow_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes flow_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = flow_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in mask_update.
function mask_update_Callback(hObject, eventdata, handles)
% hObject    handle to mask_update (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

update_mask(handles);


function update_mask(handles)

global sMAG;
global sCD;
global sMASK;
global vis_axis;
global vis_alpha;
global wss_axis;
global vis_thresh;
global VELX;
global VELY;
global VELZ;
global m_xlength;
global m_ylength;
global m_zlength;
global XPTS;
global YPTS;
global ZPTS;

global verts;
global norms;
global norm_handle;
global hpatch;

mask_type = get(handles.mask_type,'Value');
%%%VISUAL METHOD DECODING
%   1 = CD
%   2 = MAG

vis_thresh = str2double(get(handles.lumen_thresh,'String'));

new_fig = ishandle(wss_axis);
figure(wss_axis);

if new_fig == 1
    cmpos = campos;
    cmva  = camva;
    zoom reset 
end  

set(gca,'CameraPositionMode','manual');
set(gca,'CameraTargetMode','manual');
set(gca,'CameraUpVectorMode','manual');
set(gca,'CameraViewAngleMode','manual');
clf;

if mask_type == 1
    hpatch = patch(isosurface(sCD,vis_thresh));
else
    hpatch = patch(isosurface(sMAG,vis_thresh));
end
colormap('jet');
reducepatch(hpatch,0.1);
set(hpatch,'FaceColor','red','EdgeColor', 'none');

if mask_type == 1
    isonormals(sCD,hpatch)
else
    isonormals(sMAG,hpatch)
end

camlight right; 
lighting flat
alpha(0.9)
set(vis_axis, 'Renderer','OpenGL')
set(vis_axis, 'RendererMode','Manual');
set(gca,'color','black');
set(gcf,'color','black');
daspect([1 1 1])

if new_fig ==1
    campos([cmpos]);
    camva([cmva]);
end

if(new_fig==0)

view([-1 -1 0]);
zoom(0.8);
xlim([1 (m_ylength)]);
ylim([1 (m_xlength)]);
zlim([1 (m_zlength)]);
scrsz = get(0,'ScreenSize');
SIZE=[(scrsz(3)*1/2-64) scrsz(4)*1/4 scrsz(3)*1/2 scrsz(4)*1/2];
set(gcf,'Position',SIZE);
set(gcf,'Name','Visualization Window');
end

hold on
axis equal tight off vis3d;
% norms = (get(hpatch,'VertexNormals'));
% verts = (get(hpatch,'Vertices'));
% 
% scl = 15
% norm_handle = plot3( [verts(:,1)  verts(:,1) + scl*norms(:,1)]',[verts(:,2) verts(:,2)+scl*norms(:,2)]',[verts(:,3)  verts(:,3) + scl*norms(:,3)]','b','LineWidth',2)

norms = (get(hpatch,'VertexNormals'));
verts = (get(hpatch,'Vertices'));


function lumen_thresh_Callback(hObject, eventdata, handles)
% hObject    handle to lumen_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lumen_thresh as text
%        str2double(get(hObject,'String')) returns contents of lumen_thresh as a double


% --- Executes during object creation, after setting all properties.
function lumen_thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lumen_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in mask_type.
function mask_type_Callback(hObject, eventdata, handles)
% hObject    handle to mask_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns mask_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from mask_type


% --- Executes during object creation, after setting all properties.
function mask_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mask_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in normal_type.
function normal_type_Callback(hObject, eventdata, handles)
% hObject    handle to normal_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns normal_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from normal_type


% --- Executes during object creation, after setting all properties.
function normal_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to normal_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in normal_plot_type.
function normal_plot_type_Callback(hObject, eventdata, handles)
% hObject    handle to normal_plot_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns normal_plot_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from normal_plot_type


% --- Executes during object creation, after setting all properties.
function normal_plot_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to normal_plot_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function poly_num_Callback(hObject, eventdata, handles)
% hObject    handle to poly_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of poly_num as text
%        str2double(get(hObject,'String')) returns contents of poly_num as a double


% --- Executes during object creation, after setting all properties.
function poly_num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to poly_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function poly_order_Callback(hObject, eventdata, handles)
% hObject    handle to poly_order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of poly_order as text
%        str2double(get(hObject,'String')) returns contents of poly_order as a double


% --- Executes during object creation, after setting all properties.
function poly_order_CreateFcn(hObject, eventdata, handles)
% hObject    handle to poly_order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function thresh_Callback(hObject, eventdata, handles)
% hObject    handle to thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thresh as text
%        str2double(get(hObject,'String')) returns contents of thresh as a double


% --- Executes during object creation, after setting all properties.
function thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function box_xslide_Callback(hObject, eventdata, handles)
% hObject    handle to box_xslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

update_box(handles);

% --- Executes during object creation, after setting all properties.
function box_xslide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to box_xslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function box_yslide_Callback(hObject, eventdata, handles)
% hObject    handle to box_yslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
update_box(handles);

% --- Executes during object creation, after setting all properties.
function box_yslide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to box_yslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function box_zslide_Callback(hObject, eventdata, handles)
% hObject    handle to box_zslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

update_box(handles);

% --- Executes during object creation, after setting all properties.
function box_zslide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to box_zslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function box_xsize_Callback(hObject, eventdata, handles)
% hObject    handle to box_xsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of box_xsize as text
%        str2double(get(hObject,'String')) returns contents of box_xsize as a double
update_box(handles);


% --- Executes during object creation, after setting all properties.
function box_xsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to box_xsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function box_ysize_Callback(hObject, eventdata, handles)
% hObject    handle to box_ysize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of box_ysize as text
%        str2double(get(hObject,'String')) returns contents of box_ysize as a double
update_box(handles);

% --- Executes during object creation, after setting all properties.
function box_ysize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to box_ysize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function box_zsize_Callback(hObject, eventdata, handles)
% hObject    handle to box_zsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of box_zsize as text
%        str2double(get(hObject,'String')) returns contents of box_zsize as a double
update_box(handles);

% --- Executes during object creation, after setting all properties.
function box_zsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to box_zsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function update_box(handles)

%%Anatomical Data
global sMAG;
global sCD;
global sMASK;

%%Visualization Stuff
global vis_axis;
global vis_alpha;
global wss_axis;
global vis_thresh;

%Velocity data
global VELX;
global VELY;
global VELZ;
global VELXt;
global VELYt;
global VELZt;

%%Matrix Size
global m_xlength;
global m_ylength;
global m_zlength;
global tframes;

%%Pressure Stuff??
global XPTS;
global YPTS;
global ZPTS;

%%Spacings
global delX;
global delY;
global delZ;
global tres;

%%Mask Stuff
global verts;
global norms;
global norm_handle;
global norm_mag_handle;
global norm_cd_handle;
global color_range;
global hpatch;
global Cdata;

%Box stuff
global box_idx;
global box_handle;
global hist_xout;
global hist_nout;

%%%Flow Parameters (control)
global CINE_flag;
global slice_CD;
global slice_FLOW;
global slice_MASK;
global thresh;

%%%%%%Delete Object or Not %%%%%%%%%%%%%%
new_fig = ishandle(wss_axis);
figure(wss_axis);

%Delete Box
if ishandle(box_handle)
    delete(box_handle);
end

%%%CINE Flag
CINE_flag = get(handles.cine_flag,'Value');

%%Threshold
thresh = str2num(get(handles.thresh,'String'));

%%Interpolation Number
pol = str2num(get(handles.interp_num,'String'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  get Box Position and Coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phix = ( (get(handles.box_xrot,'Value')*pi - pi/2 ));
phiy = ( (get(handles.box_yrot,'Value')*pi - pi/2 ));
phiz = ( (get(handles.box_zrot,'Value')*pi - pi/2 ));

xpos = floor(1 + (get(handles.box_xslide,'Value')*( m_ylength-1) ));
ypos = floor(1 + (get(handles.box_yslide,'Value')*( m_xlength-1) ));
zpos = floor(1 + (get(handles.box_zslide,'Value')*( m_zlength-1) ));

xsize = str2num(get(handles.box_xsize,'String'));
ysize = str2num(get(handles.box_ysize,'String'));
zsize = str2num(get(handles.box_zsize,'String'));

[xB,yB,zB] = meshgrid(linspace(-xsize,xsize,2), ...
                      linspace(-ysize,ysize,2), ...
                      linspace(-zsize,zsize,2));
Tes = [ 1 2 3 4;
        5 6 7 8;
        1 3 5 7;
        2 4 6 8;
        3 4 7 8; 
        1 2 5 6];

Rx = [ 1 0        0; 
       0 cos(phix) sin(phix);   
       0 -sin(phix) cos(phix)];
Ry = [ cos(phiy) 0   -sin(phiy); 
       0         1        0;   
       sin(phiy) 0    cos(phiy)];
Rz = [ cos(phiz) sin(phiz) 0;   
      -sin(phiz) cos(phiz) 0;
      0 0 1];
  
RT = Rx*Ry*Rz;      

XB = [xB(:) yB(:) zB(:)]';
size(XB);
XB = ( RT*XB )' + repmat([xpos ypos zpos],[8 1]);;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Visualize The Box
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Visualize Control
box_on = get(handles.box_on,'Value');

%%%VISUAL Box Color DECODING
%   1 = White
%   2 = Blue
%   3 = Yellow
%   4 = Green
%   5 = Red
box_color= get(handles.box_color,'Value');
if box_color == 1 
    col = 'w';
elseif box_color == 2
    col = 'b';
elseif box_color == 3
    col = 'y';
elseif box_color == 4
    col = 'g';
elseif box_color == 5
    col = 'r';
end 

box_alpha=1.0 - 0.25*(get(handles.box_alpha,'Value')-1);

if(box_on)
    box_handle = tetramesh(Tes,XB,'FaceAlpha',box_alpha,'FaceColor',col,'EdgeColor','none');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Resample the Image in Box Coordinates (No Cine Yet)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%Find Indices Within the box
[xB,yB,zB] = meshgrid((-xsize:xsize),(-ysize:ysize),(-zsize:zsize));

Rx = [ 1 0        0; 
       0 cos(phix) sin(phix);   
       0 -sin(phix) cos(phix)];
Ry = [ cos(phiy) 0   -sin(phiy); 
       0         1        0;   
       sin(phiy) 0    cos(phiy)];
Rz = [ cos(phiz) sin(phiz) 0;   
      -sin(phiz) cos(phiz) 0;
      0 0 1];
  
RT = Rx*Ry*Rz;      

XB = [xB(:) yB(:) zB(:)]';
XB = ( RT*XB )' + repmat([xpos ypos zpos],[numel(XB)/3 1]);;

slice_CD = zeros(size(xB));
vx_slice = zeros(size(xB));
vy_slice = zeros(size(xB));
vz_slice = zeros(size(xB));

for pos=1:numel(xB)
    X = XB(pos,:);
    slice_CD(pos) = lin3d(sCD,X(2),X(1),X(3));
    vx_slice(pos) = lin3d(VELX,X(2),X(1),X(3));
    vy_slice(pos) = lin3d(VELY,X(2),X(1),X(3));
    vz_slice(pos) = lin3d(VELZ,X(2),X(1),X(3));
end

if size(xB,3)==1
slice_CD = interp2(slice_CD,pol);
vx_slice = interp2(vx_slice,pol);
vy_slice = interp2(vy_slice,pol);
vz_slice = interp2(vz_slice,pol);
else
slice_CD = interp3(slice_CD,pol);
vx_slice = interp3(vx_slice,pol);
vy_slice = interp3(vy_slice,pol);
vz_slice = interp3(vz_slice,pol);
end

%%%Test Plot Normal
Norm = RT*[0; 0; 1];
slice_FLOW = Norm(1)*vy_slice + Norm(2)*vx_slice + Norm(3)*vz_slice;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Any Processing 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Get the Mask
threshM = thresh/100*max(slice_CD(:));
slice_MASK = slice_CD > threshM;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Update the Figures in the GUI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Update Figure
mapbw = gray(128);
mapcol = jet(128);
map = [mapbw; mapcol];
axes(handles.cd_axes);
colormap(map);
im = mean(slice_CD,3);
im = 511*im/max(im(:));
im = (mean(slice_MASK,3) == 0).*im +  (mean(slice_MASK,3) > 0)*768;
imagesc(im,[0 1024]);
set(gca,'YTick',[])
set(gca,'XTick',[])
xlabel('');
ylabel('');
daspect([1 1 1]);

%%%Update Figure
axes(handles.flow_axes);
im = mean(slice_FLOW,3);
im = im/max(abs(im(:)));
im = 768 + im*255;
imagesc(im,[0 1024]);
set(gca,'YTick',[])
set(gca,'XTick',[])
xlabel('');
ylabel('');
colorbar;
daspect([1 1 1]);


% 
% box_idx = find( (abs(rverts(:,1)) < xsize) & (abs(rverts(:,2)) < ysize) & (abs(rverts(:,3)) < zsize));
% 
% if length(box_idx) > 0
% [hist_nout,hist_xout] = hist(wss(1,box_idx),str2num(get(handles.hist_bins,'String')));
% axes(handles.cd_axes)
% bar(hist_xout,hist_nout);
% xlim([0 max(wss(1,box_idx))]);
% xlabel('Avg WSS');
% ylabel('Number of Points');
% end


% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5


% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function box_filename_Callback(hObject, eventdata, handles)
% hObject    handle to box_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of box_filename as text
%        str2double(get(hObject,'String')) returns contents of box_filename as a double


% --- Executes during object creation, after setting all properties.
function box_filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to box_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in box_save.
function box_save_Callback(hObject, eventdata, handles)
% hObject    handle to box_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global wss;
global wsst;
global osi;
global box_idx;
global hist_nout;
global hist_xout;


base_name = get(handles.box_filename,'string');

hist_name = [base_name,'.hist'];
raw_name  = [base_name,'.box_vals'];

wss_out = wss(1,box_idx);


dlmwrite(hist_name,[hist_xout; hist_nout]','\t');
dlmwrite(raw_name,[1:length(box_idx); wss_out]','\t');


% fid = fopen(hist_name,'w');
% fwrite(fid,wss(box_idx),'float');
% fclose(fid);
% 
% fid = fopen(raw_name,'w');
% fwrite(fid,wss(box_idx),'float');
% fclose(fid);



% --- Executes on selection change in box_color.
function box_color_Callback(hObject, eventdata, handles)
% hObject    handle to box_color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns box_color contents as cell array
%        contents{get(hObject,'Value')} returns selected item from box_color


update_box(handles);


% --- Executes during object creation, after setting all properties.
function box_color_CreateFcn(hObject, eventdata, handles)
% hObject    handle to box_color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in box_on.
function box_on_Callback(hObject, eventdata, handles)
% hObject    handle to box_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of box_on
update_box(handles)




function hist_bin_Callback(hObject, eventdata, handles)
% hObject    handle to hist_bin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hist_bin as text
%        str2double(get(hObject,'String')) returns contents of hist_bin as a double
update_box(handles)

% --- Executes during object creation, after setting all properties.
function hist_bin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hist_bin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in box_alpha.
function box_alpha_Callback(hObject, eventdata, handles)
% hObject    handle to box_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns box_alpha contents as cell array
%        contents{get(hObject,'Value')} returns selected item from box_alpha
update_box(handles)

% --- Executes during object creation, after setting all properties.
function box_alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to box_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function box_xrot_Callback(hObject, eventdata, handles)
% hObject    handle to box_xrot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
update_box(handles)

% --- Executes during object creation, after setting all properties.
function box_xrot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to box_xrot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

set(hObject,'Value',0.5);


% --- Executes on slider movement.
function box_yrot_Callback(hObject, eventdata, handles)
% hObject    handle to box_yrot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
update_box(handles)

% --- Executes during object creation, after setting all properties.
function box_yrot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to box_yrot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


set(hObject,'Value',0.5);

% --- Executes on slider movement.
function box_zrot_Callback(hObject, eventdata, handles)
% hObject    handle to box_zrot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
update_box(handles)

% --- Executes during object creation, after setting all properties.
function box_zrot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to box_zrot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

set(hObject,'Value',0.5);



% --- Executes on button press in auto_center.
function auto_center_Callback(hObject, eventdata, handles)
% hObject    handle to auto_center (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global sMAG;
global sCD;
global sMASK;
global vis_axis;
global vis_alpha;
global wss_axis;
global vis_thresh;
global VELX;
global VELY;
global VELZ;
global m_xlength;
global m_ylength;
global m_zlength;
global XPTS;
global YPTS;
global ZPTS;
global tframes;
global delX;
global delY;
global delZ;
global verts;
global norms;
global norm_handle;
global norm_mag_handle;
global norm_cd_handle;
global color_range;
global hpatch;
global Cdata;
global visc;
global wss;
global wsst;
global osi;
global box_idx;
global box_handle;
global hist_xout;
global hist_nout;

xpos = floor(1 + (get(handles.box_xslide,'Value')*( m_ylength-1) ));
ypos = floor(1 + (get(handles.box_yslide,'Value')*( m_xlength-1) ));
zpos = floor(1 + (get(handles.box_zslide,'Value')*( m_zlength-1) ));

xsize = str2num(get(handles.box_xsize,'String'));
ysize = str2num(get(handles.box_ysize,'String'));
zsize = str2num(get(handles.box_zsize,'String'));

phix = ( (get(handles.box_xrot,'Value')*pi/2 ));
phiy = ( (get(handles.box_yrot,'Value')*pi/2 ));
phiz = ( (get(handles.box_zrot,'Value')*pi/2 ));

[xB,yB,zB] = meshgrid((-xsize:xsize),(-ysize:ysize),(-zsize:zsize));

Rx = [ 1 0        0; 
       0 cos(phix) sin(phix);   
       0 -sin(phix) cos(phix)];
Ry = [ cos(phiy) 0   -sin(phiy); 
       0         1        0;   
       sin(phiy) 0    cos(phiy)];
Rz = [ cos(phiz) sin(phiz) 0;   
      -sin(phiz) cos(phiz) 0;
      0 0 1];
  
RT = Rx*Ry*Rz;      

XB = [xB(:) yB(:) zB(:)]';
size(XB)
XB = ( RT*XB )' + repmat([xpos ypos zpos],[numel(XB)/3 1]);;

%%%%Find Values in Box
COMx=0;
COMy=0;
COMz=0;
M0=0;

for pos=1:numel(xB)
    X = XB(pos,:);
    V = lin3d(sCD,X(2),X(1),X(3));
    M0 = M0 + V;
    COMx = X(1)*V + COMx;
    COMy = X(2)*V + COMy;
    COMz = X(3)*V + COMz;
end
COMx= COMx/M0
COMy= COMy/M0
COMz= COMz/M0
set(handles.box_xslide,'Value',round(COMx-1)/( m_xlength-1));
set(handles.box_yslide,'Value',round(COMy-1)/( m_ylength-1));
set(handles.box_zslide,'Value',round(COMz-1)/( m_zlength-1));
update_box(handles)

% --- Executes on button press in auto_align.
function auto_align_Callback(hObject, eventdata, handles)
% hObject    handle to auto_align (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global sMAG;
global sCD;
global sMASK;
global vis_axis;
global vis_alpha;
global wss_axis;
global vis_thresh;
global VELX;
global VELY;
global VELZ;
global m_xlength;
global m_ylength;
global m_zlength;
global XPTS;
global YPTS;
global ZPTS;
global tframes;
global delX;
global delY;
global delZ;
global verts;
global norms;
global norm_handle;
global norm_mag_handle;
global norm_cd_handle;
global color_range;
global hpatch;
global Cdata;
global visc;
global wss;
global wsst;
global osi;
global box_idx;
global box_handle;
global hist_xout;
global hist_nout;

xpos = floor(1 + (get(handles.box_xslide,'Value')*( m_ylength-1) ));
ypos = floor(1 + (get(handles.box_yslide,'Value')*( m_xlength-1) ));
zpos = floor(1 + (get(handles.box_zslide,'Value')*( m_zlength-1) ));

xsize = str2num(get(handles.box_xsize,'String'));
ysize = str2num(get(handles.box_ysize,'String'));
zsize = str2num(get(handles.box_zsize,'String'));

phix = ( (get(handles.box_xrot,'Value')*pi - pi/2 ));
phiy = ( (get(handles.box_yrot,'Value')*pi - pi/2 ));
phiz = ( (get(handles.box_zrot,'Value')*pi - pi/2 ));

MIN_ERR = 1e99;

for pass =1:2
    if pass ==1
        res = 5;
        range = linspace( -pi/4,pi/4,res);
        phix0 = phix;
        phiy0 = phiy;
    else
        res = 5;
        maxp = pi / res / 4;
        range = linspace( -maxp,maxp,res);
        phix0 = phix;
        phiy0 = phiy;
    end

for cphix = range + phix0
 for cphiy = range + phiy0
  
    [xB,yB,zB] = meshgrid((-xsize:xsize),(-ysize:ysize),(-zsize:zsize));
      
    Rx = [ 1 0        0; 
       0 cos(cphix) sin(cphix);   
       0 -sin(cphix) cos(cphix)];
    Ry = [ cos(cphiy) 0   -sin(cphiy); 
       0         1        0;   
       sin(cphiy) 0    cos(cphiy)];
    Rz = [ cos(phiz) sin(phiz) 0;   
      -sin(phiz) cos(phiz) 0;
      0 0 1];
    RT = Rx*Ry*Rz;      

    XB = [xB(:) yB(:) zB(:)]';
    XB = ( RT*XB )' + repmat([xpos ypos zpos],[numel(XB)/3 1]);;
    
    for slice = 1:size(xB,3)
        COMx(slice)=0;
        COMy(slice)=0;
        M0=0;
        
        for x=1:xsize*2+1
         for y=1:ysize*2+1
            pos = sub2ind(size(xB),y,x,slice);
            X = XB(pos,:);
            V = lin3d(sCD,X(2),X(1),X(3));
            M0 = M0 + V;
            COMx(slice) =x*V + COMx(slice);
            COMy(slice) =y*V + COMy(slice);
         end
        end
        COMx(slice) = COMx(slice)/M0;
        COMy(slice) = COMy(slice)/M0;
    end
    ERR = abs(std(COMx)^2 + std(COMy)^2);
    
    if ERR < MIN_ERR
        phix = cphix;
        phiy = cphiy;
        MIN_ERR =ERR;
    end
 end    
end
disp(['Pass ',int2str(pass),': ',num2str(MIN_ERR)]);

end

set(handles.box_xrot,'Value',(phix+pi/2)/(pi));
set(handles.box_yrot,'Value',(phiy+pi/2)/(pi));

update_box(handles)


% --- Executes on selection change in cine_flag.
function cine_flag_Callback(hObject, eventdata, handles)
% hObject    handle to cine_flag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns cine_flag contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cine_flag


% --- Executes during object creation, after setting all properties.
function cine_flag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cine_flag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function interp_num_Callback(hObject, eventdata, handles)
% hObject    handle to interp_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of interp_num as text
%        str2double(get(hObject,'String')) returns contents of interp_num as a double


% --- Executes during object creation, after setting all properties.
function interp_num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to interp_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


