function varargout = pwv_tool(varargin)
% PWV_TOOL M-file for pwv_tool.fig
%      PWV_TOOL, by itself, creates a new PWV_TOOL or raises the existing
%      singleton*.
%
%      H = PWV_TOOL returns the handle to a new PWV_TOOL or the handle to
%      the existing singleton*.
%
%      PWV_TOOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PWV_TOOL.M with the given input
%      arguments.
%
%      PWV_TOOL('Property','Value',...) creates a new PWV_TOOL or raises
%      the
%      existing singleton*.  Starting from the left, property value pairs
%      are
%      applied to the GUI before pwv_tool_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pwv_tool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pwv_tool

% Last Modified by GUIDE v2.5 03-Sep-2010 11:52:38
% Written by Andrew Wentland, June 7, 2010

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pwv_tool_OpeningFcn, ...
                   'gui_OutputFcn',  @pwv_tool_OutputFcn, ...
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


% --- Executes just before pwv_tool is made visible.
function pwv_tool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pwv_tool (see VARARGIN)

% Choose default command line output for pwv_tool
handles.output = hObject;

global plane_counter; plane_counter = 0;
global check_counter; check_counter = 1;
global pwv_3d_computed; pwv_3d_computed = 0;

% Update handles structure
guidata(hObject, handles);
% UIWAIT makes pwv_tool wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = pwv_tool_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in exit.
function exit_Callback(hObject, eventdata, handles)
% hObject    handle to exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close

function lumen_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to lumen_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lumen_threshold as text
%        str2double(get(hObject,'String')) returns contents of lumen_threshold as a double

% --- Executes during object creation, after setting all properties.
function lumen_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lumen_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in update_image.
function update_image_Callback(hObject, eventdata, handles)
% hObject    handle to update_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update_mask(handles);

function update_mask(handles)
global sCD vis_thresh verts norms hpatch;

vis_thresh = str2double(get(handles.lumen_threshold,'String'));
set(handles.mip,'Visible','on');
axes(handles.mip)
set(gcf, 'Renderer','OpenGL') %2010-10-19 renderer change
opengl software  %2010-10-19 renderer change
% 2010-09-30 correct 3D display to show normal anatomic display
set(gca,'YDir','reverse','ZDir','reverse');

hpatch = patch(isosurface(sCD,vis_thresh));

colormap('jet'); % might not be alright
reducepatch(hpatch,0.1); % might not be alright
set(hpatch,'FaceColor','red','EdgeColor', 'none');

isonormals(sCD,hpatch)

camlight right; 
lighting gouraud
alpha(0.9)
set(gcf, 'Renderer','OpenGL')
opengl software % 2010-10-06 Text is upside down likely due to a rendering problem with OpenGL. Setting the renderer to software avoids this problem. -alw
set(gcf, 'RendererMode','Manual');
set(gca,'Color',[0.6,0.8,1]);
daspect([1 1 1])

norms = (get(hpatch,'VertexNormals'));
verts = (get(hpatch,'Vertices'));

xlabel('X');ylabel('Y');zlabel('Z');

function update_box(handles)

% Anatomical Data
global sMAG sCD area;

% Velocity data
global VELXt VELYt VELZt;

% Matrix Size
global m_xlength m_ylength m_zlength tframes;

% Spacings (voxel sizes)
global delX delY;
global tres; %% length of tframes in ms

% Plane Variable
global plane_handle;
global check_counter;
global plane_counter plane_color;
global overall_position plane_coordinates;
global centerline_points_3d;

% Flow Parameters (control)
global slice_CD slice_VEL slice_MASK;
global thresh;

% Delete Box
if check_counter == plane_counter
    if isempty(plane_handle)
    else
        delete(plane_handle(plane_counter).handle);
    end
else
    check_counter = plane_counter;
end

plane_color = ['k' 'b' 'g' 'r' 'c' 'm']; % black, blue, green, red, light blue, and purple, respectively

% Threshold
thresh = str2double(get(handles.thresh,'String'));

% Interpolation Number
pol = str2double(get(handles.interp_num,'String'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Get Box Position and Coordinates      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phix = ((get(handles.box_xrot,'Value')*pi - pi/2));
phiy = ((get(handles.box_yrot,'Value')*pi - pi/2));
phiz = ((get(handles.box_zrot,'Value')*pi - pi/2));

xpos = floor(1 + (get(handles.box_xslide,'Value')*( m_ylength-1)));
ypos = floor(1 + (get(handles.box_yslide,'Value')*( m_xlength-1)));
zpos = floor(1 + (get(handles.box_zslide,'Value')*( m_zlength-1)));

xsize = str2double(get(handles.box_xsize,'String'));
ysize = str2double(get(handles.box_ysize,'String'));
zsize = str2double(get(handles.box_zsize,'String'));

hold on
axes(handles.mip);
set(gcf, 'Renderer','OpenGL') %2010-10-19 renderer change
opengl software  %2010-10-19 renderer change
[xB,yB,zB] = meshgrid(linspace(-xsize,xsize,2), ...
             linspace(-ysize,ysize,2), ...
             linspace(-zsize,zsize,2));
Tes = [ 1 2 3 4;
        5 6 7 8;
        1 3 5 7;
        2 4 6 8;
        3 4 7 8;
        1 2 5 6];

Rx = [ 1 0 0;
    0 cos(phix) sin(phix);
    0 -sin(phix) cos(phix)];
Ry = [cos(phiy) 0 -sin(phiy);
    0         1        0;
    sin(phiy) 0    cos(phiy)];
Rz = [ cos(phiz) sin(phiz) 0;
    -sin(phiz) cos(phiz) 0;
    0 0 1];

RT = Rx*Ry*Rz;

XB = [xB(:) yB(:) zB(:)]';
size(XB);
XB = (RT*XB )' + repmat([xpos ypos zpos],[8 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       Create the Plane for Analysis       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plane_alpha = 1;

% Create the plane
axes(handles.mip)
set(gcf, 'Renderer','OpenGL') %2010-10-19 renderer change
opengl software  %2010-10-19 renderer change
if plane_counter < 7
    plane_handle(plane_counter).handle = tetramesh(Tes,XB,'FaceAlpha',plane_alpha,'FaceColor',plane_color(plane_counter),'EdgeColor','none');
elseif plane_counter > 6
    plane_handle(plane_counter).handle = tetramesh(Tes,XB,'FaceAlpha',plane_alpha,'FaceColor',plane_color(plane_counter-6),'EdgeColor','none');
elseif plane_counter > 12
    plane_handle(plane_counter).handle = tetramesh(Tes,XB,'FaceAlpha',plane_alpha,'FaceColor',plane_color(plane_counter-12),'EdgeColor','none');
end

plane_handle(plane_counter).xloc = get(handles.box_xslide,'Value');
plane_handle(plane_counter).yloc = get(handles.box_yslide,'Value');
plane_handle(plane_counter).zloc = get(handles.box_zslide,'Value');
plane_handle(plane_counter).xrot = get(handles.box_xrot,'Value');
plane_handle(plane_counter).yrot = get(handles.box_yrot,'Value');
plane_handle(plane_counter).zrot = get(handles.box_zrot,'Value');
plane_handle(plane_counter).xsize = str2double(get(handles.box_xsize,'String'));
plane_handle(plane_counter).ysize = str2double(get(handles.box_ysize,'String'));
plane_handle(plane_counter).zsize = str2double(get(handles.box_zsize,'String'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Resample the Image in Box Coordinates (No Cine Yet)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%Find Indices Within the box
[xB,yB,zB] = meshgrid(linspace(-xsize,xsize,pol*xsize*2+1),linspace(-ysize,ysize,pol*ysize*2+1),linspace(-zsize,zsize,pol*zsize*2+1));

Rx = [ 1 0 0;
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
XB = ( RT*XB )' + repmat([xpos ypos zpos],[numel(XB)/3 1]);

% CINE data set processing
slice_CD = zeros(size(xB));
slice_VEL = zeros(size(xB,1),size(xB,2),size(xB,3),tframes);
vx_slice = zeros(size(xB,1),size(xB,2),size(xB,3),tframes);
vy_slice = zeros(size(xB,1),size(xB,2),size(xB,3),tframes);
vz_slice = zeros(size(xB,1),size(xB,2),size(xB,3),tframes);
vtemp = zeros(size(xB,1),size(xB,2),size(xB,3));
cdtemp = zeros(size(xB,1),size(xB,2),size(xB,3));
magtemp = zeros(size(xB,1),size(xB,2),size(xB,3));

% Get a starting point
XF = floor(XB(:,1)); YF = floor(XB(:,2)); ZF = floor(XB(:,3));

% Find points outside the matrix
MASK = (XF < 1 ) | (YF < 1) | (ZF < 1) | (XF > m_ylength-1) | (YF > m_xlength-1) | (ZF > m_zlength-1);
idx2 = find(MASK == 0);
XF = XF(idx2); YF = YF(idx2); ZF = ZF(idx2);

dx = XB(idx2,1) - XF; dy = XB(idx2,2) - YF; dz = XB(idx2,3) - ZF;
mdx = 1 - dx; mdy = 1 - dy; mdz = 1 - dz;

CDtemp = sCD(:,:,:);
idx = sub2ind(size(CDtemp),YF,XF,ZF);
cdtemp(idx2)= cdtemp(idx2)+  mdy.*mdx.*mdz.*CDtemp(idx);
idx = sub2ind(size(CDtemp),YF,XF+1,ZF);
cdtemp(idx2)= cdtemp(idx2)+  mdy.*dx.*mdz.*CDtemp(idx);
idx = sub2ind(size(CDtemp),YF+1,XF,ZF);
cdtemp(idx2)= cdtemp(idx2)+  dy.*mdx.*mdz.*CDtemp(idx);
idx = sub2ind(size(CDtemp),YF+1,XF+1,ZF);
cdtemp(idx2)= cdtemp(idx2)+  dy.*dx.*mdz.*CDtemp(idx);
idx = sub2ind(size(CDtemp),YF,XF,ZF+1);
cdtemp(idx2)= cdtemp(idx2)+  mdy.*mdx.*dz.*CDtemp(idx);
idx = sub2ind(size(CDtemp),YF,XF+1,ZF+1);
cdtemp(idx2)= cdtemp(idx2)+  mdy.*dx.*dz.*CDtemp(idx);
idx = sub2ind(size(CDtemp),YF+1,XF,ZF+1);
cdtemp(idx2)= cdtemp(idx2)+  dy.*mdx.*dz.*CDtemp(idx);
idx = sub2ind(size(CDtemp),YF+1,XF+1,ZF+1);
cdtemp(idx2)= cdtemp(idx2)+  dy.*dx.*dz.*CDtemp(idx);
slice_CD(:,:,:) = cdtemp;

MAGtemp = sMAG(:,:,:);
idx = sub2ind(size(MAGtemp),YF,XF,ZF);
magtemp(idx2)= magtemp(idx2)+  mdy.*mdx.*mdz.*MAGtemp(idx);
idx = sub2ind(size(MAGtemp),YF,XF+1,ZF);
magtemp(idx2)= magtemp(idx2)+  mdy.*dx.*mdz.*MAGtemp(idx);
idx = sub2ind(size(MAGtemp),YF+1,XF,ZF);
magtemp(idx2)= magtemp(idx2)+  dy.*mdx.*mdz.*MAGtemp(idx);
idx = sub2ind(size(MAGtemp),YF+1,XF+1,ZF);
magtemp(idx2)= magtemp(idx2)+  dy.*dx.*mdz.*MAGtemp(idx);
idx = sub2ind(size(MAGtemp),YF,XF,ZF+1);
magtemp(idx2)= magtemp(idx2)+  mdy.*mdx.*dz.*MAGtemp(idx);
idx = sub2ind(size(MAGtemp),YF,XF+1,ZF+1);
magtemp(idx2)= magtemp(idx2)+  mdy.*dx.*dz.*MAGtemp(idx);
idx = sub2ind(size(MAGtemp),YF+1,XF,ZF+1);
magtemp(idx2)= magtemp(idx2)+  dy.*mdx.*dz.*MAGtemp(idx);
idx = sub2ind(size(MAGtemp),YF+1,XF+1,ZF+1);
magtemp(idx2)= magtemp(idx2)+  dy.*dx.*dz.*MAGtemp(idx);
slice_MAG(:,:,:)=magtemp;

for t=1:tframes

    %interpolation in x
	vtemp(:) = 0;
	Vtemp = VELXt(:,:,:,t);
	idx = sub2ind(size(Vtemp),YF,XF,ZF);
	vtemp(idx2)= vtemp(idx2)+  mdy.*mdx.*mdz.*Vtemp(idx);
	idx = sub2ind(size(Vtemp),YF,XF+1,ZF);
	vtemp(idx2)= vtemp(idx2)+  mdy.*dx.*mdz.*Vtemp(idx);
	idx = sub2ind(size(Vtemp),YF+1,XF,ZF);
	vtemp(idx2)= vtemp(idx2)+  dy.*mdx.*mdz.*Vtemp(idx);
	idx = sub2ind(size(Vtemp),YF+1,XF+1,ZF);
	vtemp(idx2)= vtemp(idx2)+  dy.*dx.*mdz.*Vtemp(idx);
	idx = sub2ind(size(Vtemp),YF,XF,ZF+1);
	vtemp(idx2)= vtemp(idx2)+  mdy.*mdz.*dz.*Vtemp(idx);
	idx = sub2ind(size(Vtemp),YF,XF+1,ZF+1);
	vtemp(idx2)= vtemp(idx2)+  mdy.*dx.*dz.*Vtemp(idx);
	idx = sub2ind(size(Vtemp),YF+1,XF,ZF+1);
	vtemp(idx2)= vtemp(idx2)+  dy.*mdx.*dz.*Vtemp(idx);
	idx = sub2ind(size(Vtemp),YF+1,XF+1,ZF+1);
    vtemp(idx2)= vtemp(idx2)+  dy.*dx.*dz.*Vtemp(idx);
    vx_slice(:,:,:,t)=vtemp;

    %interpolation in y
    vtemp(:) = 0;
    Vtemp = VELYt(:,:,:,t);
    idx = sub2ind(size(Vtemp),YF,XF,ZF);
    vtemp(idx2)= vtemp(idx2)+  mdy.*mdx.*mdz.*Vtemp(idx);
    idx = sub2ind(size(Vtemp),YF,XF+1,ZF);
    vtemp(idx2)= vtemp(idx2)+  mdy.*dx.*mdz.*Vtemp(idx);
    idx = sub2ind(size(Vtemp),YF+1,XF,ZF);
    vtemp(idx2)= vtemp(idx2)+  dy.*mdx.*mdz.*Vtemp(idx);
    idx = sub2ind(size(Vtemp),YF+1,XF+1,ZF);
    vtemp(idx2)= vtemp(idx2)+  dy.*dx.*mdz.*Vtemp(idx);
    idx = sub2ind(size(Vtemp),YF,XF,ZF+1);
    vtemp(idx2)= vtemp(idx2)+  mdy.*mdx.*dz.*Vtemp(idx);
    idx = sub2ind(size(Vtemp),YF,XF+1,ZF+1);
    vtemp(idx2)= vtemp(idx2)+  mdy.*dx.*dz.*Vtemp(idx);
    idx = sub2ind(size(Vtemp),YF+1,XF,ZF+1);
    vtemp(idx2)= vtemp(idx2)+  dy.*mdx.*dz.*Vtemp(idx);
    idx = sub2ind(size(Vtemp),YF+1,XF+1,ZF+1);
    vtemp(idx2)= vtemp(idx2)+  dy.*dx.*dz.*Vtemp(idx);
    vy_slice(:,:,:,t)=vtemp;

    %interpolation in z
    vtemp(:) = 0;
    Vtemp = VELZt(:,:,:,t);
    idx = sub2ind(size(Vtemp),YF,XF,ZF);
    vtemp(idx2)= vtemp(idx2)+  mdx.*mdy.*mdz.*Vtemp(idx);
    idx = sub2ind(size(Vtemp),YF,XF+1,ZF);
    vtemp(idx2)= vtemp(idx2)+  mdy.*dx.*mdz.*Vtemp(idx);
    idx = sub2ind(size(Vtemp),YF+1,XF,ZF);
    vtemp(idx2)= vtemp(idx2)+  dy.*mdx.*mdz.*Vtemp(idx);
    idx = sub2ind(size(Vtemp),YF+1,XF+1,ZF);
    vtemp(idx2)= vtemp(idx2)+  dy.*dx.*mdz.*Vtemp(idx);
    idx = sub2ind(size(Vtemp),YF,XF,ZF+1);
    vtemp(idx2)= vtemp(idx2)+  mdy.*mdx.*dz.*Vtemp(idx);
    idx = sub2ind(size(Vtemp),YF,XF+1,ZF+1);
    vtemp(idx2)= vtemp(idx2)+  mdy.*dx.*dz.*Vtemp(idx);
    idx = sub2ind(size(Vtemp),YF+1,XF,ZF+1);
    vtemp(idx2)= vtemp(idx2)+  dy.*mdx.*dz.*Vtemp(idx);
    idx = sub2ind(size(Vtemp),YF+1,XF+1,ZF+1);
    vtemp(idx2)= vtemp(idx2)+  dy.*dx.*dz.*Vtemp(idx);
    vz_slice(:,:,:,t)=vtemp;
end

% Test Plot Normal
Norm = RT*[0; 0; 1];
slice_VEL = Norm(1)*vy_slice + Norm(2)*vx_slice + Norm(3)*vz_slice;

threshM = thresh/100*max(slice_CD(:));
slice_MASK = slice_CD > threshM;
%roi_idx = find(slice_MASK > 0);
%voxels = length(roi_idx(:,:,1));

size(slice_CD);
%size_mask = size(slice_MASK);
%slice_mask = slice_MASK;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       Update the Figures in the GUI        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update anatomical figure
mapbw = gray(128);
mapcol = jet(128);
map = [mapbw; mapcol];
axes(handles.cd_axes);
set(gcf, 'Renderer','OpenGL') %2010-10-19 renderer change
opengl software  %2010-10-19 renderer change
axis square; % alw 2010-10-06
colormap(map);
im = mean(slice_CD,3);
im = 511*im/max(im(:));
im = (mean(slice_MASK,3) == 0).*im +  (mean(slice_MASK,3) > 0)*768;
cla;
imagesc(im,[0 1024]);
set(gca,'YTick',[])
set(gca,'XTick',[])
xlabel('');
ylabel('');
daspect([1 1 1]);

% Centroid calculations.
% This is extremely convoluted; there's got to be a better way of doing
% this.
% Since MATLAB won't automatically decrease if the value on the left side
% of the colon is greater than the value on the right side. Thus it is
% necessary to ask first if the left side is less than the right side. If
% not, then a decrease by -1 in the for loop needs to be specified, as seen
% below.
% The xpos, ypos, and zpos values are the coordinates of the center of the
% plane. Therefore, to scan across the region of the plane and where that
% plane intersects the surface rendering, it is necessary to back up from
% that each coordinate by half of the length of the plane in each
% dimension.
% alw
x = 0;y = 0;z = 0;
threshM = thresh/100*max(max(max(sCD)));
box_size = RT*[2*xsize; 2*ysize; 2*zsize];
if xpos < xpos+box_size(1)
    for j = xpos-floor(box_size(1)/2):xpos+floor(box_size(1)/2)
        if ypos < ypos+box_size(2)
            for i = ypos-floor(box_size(2)/2):ypos+floor(box_size(2)/2)
                if zpos < zpos+box_size(3)
                    for k = zpos-floor(box_size(3)/2):zpos+floor(box_size(3)/2)
                        if sCD(i,j,k) > threshM
                            x = [x i];y = [y j]; z = [z k];
                        end
                    end
                else
                    for k = zpos-floor(box_size(3)/2):-1:zpos+floor(box_size(3)/2)
                        if sCD(i,j,k) > threshM
                            x = [x i];y = [y j]; z = [z k];
                        end
                    end                
                end
            end
        else
            for i = ypos-floor(box_size(2)/2):-1:ypos+floor(box_size(2)/2)
                if zpos < zpos+box_size(3)
                    for k = zpos-floor(box_size(3)/2):zpos+floor(box_size(3)/2)
                        if sCD(i,j,k) > threshM
                            x = [x i];y = [y j]; z = [z k];
                        end
                    end
                else
                    for k = zpos-floor(box_size(3)/2):-1:zpos+floor(box_size(3)/2)
                        if sCD(i,j,k) > threshM
                            x = [x i];y = [y j]; z = [z k];
                        end
                    end            
                end
            end
        end
    end
else
    for j = xpos-floor(box_size(1)/2):-1:xpos+floor(box_size(1)/2)
        if ypos < ypos+box_size(2)
            for i = ypos-floor(box_size(2)/2):ypos+floor(box_size(2)/2)
                if zpos < zpos+box_size(3)
                    for k = zpos-floor(box_size(3)/2):zpos+floor(box_size(3)/2)
                        if sCD(i,j,k) > threshM
                            x = [x i];y = [y j]; z = [z k];
                        end
                    end
                else
                    for k = zpos-floor(box_size(3)/2):-1:zpos+floor(box_size(3)/2)
                        if sCD(i,j,k) > threshM
                            x = [x i];y = [y j]; z = [z k];
                        end
                    end                
                end
            end
        else
            for i = ypos-floor(box_size(2)/2):-1:ypos+floor(box_size(2)/2)
                if zpos < zpos+box_size(3)
                    for k = zpos-floor(box_size(3)/2):zpos+floor(box_size(3)/2)
                        if sCD(i,j,k) > threshM
                            x = [x i];y = [y j]; z = [z k];
                        end
                    end
                else
                    for k = zpos-floor(box_size(3)/2):-1:zpos+floor(box_size(3)/2)
                        if sCD(i,j,k) > threshM
                            x = [x i];y = [y j]; z = [z k];
                        end
                    end            
                end
            end
        end
    end
end

x = x(2:length(x)); y = y(2:length(y)); z = z(2:length(z));
overall_position = [floor(mean(y)) floor(mean(x)) floor(mean(z))];
plane_coordinates(plane_counter,:,:) = overall_position;

% Update velocity figure
axes(handles.flow_axes);
set(gcf, 'Renderer','OpenGL') %2010-10-19 renderer change
opengl software  %2010-10-19 renderer change
axis square; % alw 2010-10-06
vel = mean(slice_VEL,4);
im = mean(vel,3);
im = im/max(abs(im(:)));
im = 768 + im*255;
imagesc(im,[0 1024]);
set(gca,'YTick',[])
set(gca,'XTick',[])
xlabel('');
ylabel('');
daspect([1 1 1]);

% Update mean flow figure
axes(handles.mean_flow);
set(gcf, 'Renderer','OpenGL') %2010-10-19 renderer change
opengl software  %2010-10-19 renderer change
opengl software % 2010-10-11
delX_pol = delX / pol;
delY_pol = delY / pol;

global flow_t flow_array;

[o,p] = size(plane_coordinates); clear p; % 2010-10-22

area = size(slice_MASK,1)*size(slice_MASK,2)*delX_pol*delY_pol/100; %2011-06-07 alw

if plane_counter == o
    for t = 1:tframes
        masked_vel_t = slice_VEL(:,:,:,t).*double(slice_MASK);
        
        [C,I] = max(abs(masked_vel_t)); clear C;
        if masked_vel_t(I(1)) < 0
            masked_vel_t = masked_vel_t.*(-1);
        end
        
        flow_t(t) = sum(masked_vel_t(:))/size(slice_MASK,3)*delX_pol*delY_pol/1000; % 2010-11-23 removed abs
        time(t) = t*tres;
    end
else
    for t = 1:tframes
        masked_vel_t = slice_VEL(:,:,:,t).*double(slice_MASK);
        
        [C,I] = max(abs(masked_vel_t)); clear C;
        if masked_vel_t(I(1)) < 0
            masked_vel_t = masked_vel_t.*(-1);
        end
        
        flow_array(plane_counter+1,t) = sum(masked_vel_t(:))/size(slice_MASK,3)*delX_pol*delY_pol/1000; % 2010-11-23 removed abs
    end
end

% alw 2011-03-17 inverts the function if the flow curve is inverted
if plane_counter == o
else
    [C,I] = max(abs(flow_array(plane_counter+1,:))); clear C;
    
    if flow_array(plane_counter+1,I) < 0
        flow_array(plane_counter+1,:) = -1.*flow_array(plane_counter+1,:);
    end
end
clear I;

[C,I] = max(abs(flow_t)); clear C;
if flow_t(I) < 0
    flow_t = -1.*flow_t;
end
clear I;

% Needed so that only the current plane information is plotted; thus,
% the mean flow plot is created for each call to this update_box function.
cla;

% Handles six different colors for up to 18 planes.
% Plots flow waveforms from all planes except for the current plane.
if size(flow_array,1) > 1
    for i = 2:size(flow_array,1)
        if i < 8
            plot(flow_array(1,:),flow_array(i,:),'Color',plane_color(i-1));% removed abs from all of the plot calls 2011-03-17
        elseif i > 7
            plot(flow_array(1,:),flow_array(i,:),'Color',plane_color(i-5));
        elseif i > 13
            plot(flow_array(1,:),flow_array(i,:),'Color',plane_color(i-11));
        end
        hold on
    end
end

% Handles six different colors for up to 18 planes.
% Plots the current flow waveform from the active plane.
o = size(plane_coordinates,1);
if o < 7
    plot(flow_array(1,:),flow_t,'Color',plane_color(o));% removed abs from all of the plot calls 2011-03-17
elseif o > 6
    plot(flow_array(1,:),flow_t,'Color',plane_color(o-6));
elseif o > 12
    plot(flow_array(1,:),flow_t,'Color',plane_color(o-12));
end
xlabel('time [ms]');
ylabel('velocity [ml/s]');

% Updates the static text boxes in the GUI with the current location and
% rotation information of the plane.
get(handles.box_xslide,'Value')
get(handles.box_yslide,'Value')
get(handles.box_zslide,'Value')
 set(handles.x_loc,'String',['X Loc: ' num2str(get(handles.box_xslide,'Value'),2)]);
 set(handles.y_loc,'String',['Y Loc: ' num2str(get(handles.box_yslide,'Value'),2)]);
 set(handles.z_loc,'String',['Z Loc: ' num2str(get(handles.box_zslide,'Value'),2)]);
 set(handles.rotx,'String',['X Rot: ' num2str(get(handles.box_xrot,'Value'),2)]);
 set(handles.roty,'String',['Y Rot: ' num2str(get(handles.box_yrot,'Value'),2)]);
 set(handles.rotz,'String',['Z Rot: ' num2str(get(handles.box_zrot,'Value'),2)]);
% End of update_box function

% --- Executes on button press in new_plane.
function new_plane_Callback(hObject, eventdata, handles)
% hObject    handle to new_plane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global plane_counter plane_info plane_handle;
global flow_t time;
global flow_array; % array of mean flow over time for each plane
global tres tframes;
global overall_position plane_coordinates;

%plane_counter = plane_counter + 1;
[x,y] = size(plane_coordinates); clear y;
plane_counter = x + 1; clear x; % 2010-10-21

if plane_counter == 2
    plane_info = [...
        get(handles.box_xslide,'Value') get(handles.box_yslide,'Value')...
        get(handles.box_zslide,'Value') get(handles.box_xrot,'Value')...
        get(handles.box_yrot,'Value') get(handles.box_zrot,'Value')...
        str2num(get(handles.box_xsize,'String')) str2num(get(handles.box_ysize,'String'))...
        str2num(get(handles.box_zsize,'String'))];
elseif plane_counter > 2
    plane_info = [plane_info;...
        get(handles.box_xslide,'Value') get(handles.box_yslide,'Value')...
        get(handles.box_zslide,'Value') get(handles.box_xrot,'Value')...
        get(handles.box_yrot,'Value') get(handles.box_zrot,'Value')...
        str2num(get(handles.box_xsize,'String')) str2num(get(handles.box_ysize,'String'))...
        str2num(get(handles.box_zsize,'String'))];
end

for t = 1:tframes
    time(t) = t*tres;
end

if plane_counter == 1
    flow_array = time;
elseif plane_counter > 1
    flow_array = [flow_array;flow_t];
end

if plane_counter == 2
    plane_coordinates = overall_position;
elseif plane_counter > 2
    plane_coordinates = [plane_coordinates;overall_position];
end

if plane_counter < 10
    text = ['0' num2str(plane_counter)];
else
    text = num2str(plane_counter);
end

a = get(handles.plane_select,'String');
% text = [num2str(plane_counter) ', ' text];
text = [a; text];

if plane_counter > 1
    set(handles.plane_select,'String',text)
end

set(handles.plane_select,'Value',plane_counter)

update_box(handles)

% --- Executes on slider movement.
function box_xslide_Callback(hObject, eventdata, handles)
% hObject    handle to box_xslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global plane_counter;

if plane_counter == 0
    disp('You have not yet created a plane.');
else
    update_box(handles);
end

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
global plane_counter;

if plane_counter == 0
    disp('You have not yet created a plane.');
else
    update_box(handles);
end

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
global plane_counter;

if plane_counter == 0
    disp('You have not yet created a plane.');
else
    update_box(handles);
end

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
global plane_counter;

if plane_counter == 0
    disp('You have not yet created a plane.');
else
    update_box(handles);
end

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
global plane_counter;

if plane_counter == 0
    disp('You have not yet created a plane.');
else
    update_box(handles);
end

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
global plane_counter;

if plane_counter == 0
    disp('You have not yet created a plane.');
else
    update_box(handles);
end

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

% --- Executes on slider movement.
function box_xrot_Callback(hObject, eventdata, handles)
% hObject    handle to box_xrot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global plane_counter;

if plane_counter == 0
    disp('You have not yet created a plane.');
else
    update_box(handles);
end

% --- Executes during object creation, after setting all properties.
function box_xrot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to box_xrot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function box_yrot_Callback(hObject, eventdata, handles)
% hObject    handle to box_yrot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global plane_counter;

if plane_counter == 0
    disp('You have not yet created a plane.');
else
    update_box(handles);
end

% --- Executes during object creation, after setting all properties.
function box_yrot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to box_yrot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function box_zrot_Callback(hObject, eventdata, handles)
% hObject    handle to box_zrot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global plane_counter;

if plane_counter == 0
    disp('You have not yet created a plane.');
else
    update_box(handles);
end

% --- Executes during object creation, after setting all properties.
function box_zrot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to box_zrot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on button press in save_plane_info.
function save_plane_info_Callback(hObject, eventdata, handles)
% hObject    handle to save_plane_info (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global plane_counter plane_info;

if plane_counter ~= 0
    plane_info = [plane_info;...
        get(handles.box_xslide,'Value') get(handles.box_yslide,'Value')...
        get(handles.box_zslide,'Value') get(handles.box_xrot,'Value')...
        get(handles.box_yrot,'Value') get(handles.box_zrot,'Value')...
        str2num(get(handles.box_xsize,'String')) str2num(get(handles.box_ysize,'String'))...
        str2num(get(handles.box_zsize,'String'))];
end

plane_info_file = get(handles.save_name,'string');

if plane_counter ~= 0
    fid = fopen([plane_info_file,'.txt'],'w');
    fprintf(fid,'Plane;X-Loc;Y-Loc;Z-Loc;Rot-X;Rot-Y;Rot-Z;Sx;Sy;Sz\n');
    
    for i = 1:size(plane_info,1) % rows
        fprintf(fid,'Plane %d;',i);
        for j = 1:size(plane_info,2) % columns
            fprintf(fid,'%f;',plane_info(i,j));
        end
        
        if i < size(plane_info,1)
            fprintf(fid,'\n');
        end
    end    
    
    fclose(fid);
else
    warndlg('You have not yet created any planes.','!! Warning !!');
end

% --- Executes on button press in load_saved_planes.
function load_saved_planes_Callback(hObject, eventdata, handles)
% hObject    handle to load_saved_planes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global plane_counter; plane_counter = 0;

[filename, pathname] = uigetfile('*.*', 'Select the text file with your plane information.');
fid = fopen([pathname filename]);
C = textscan(fid, '%s %s %s %s %s %s %s %s %s %s', 'delimiter', ';');
fclose(fid);
x = size(C{1},1);

for i = 2:x 
    set(handles.box_xslide,'Value',str2num(C{2}{i}));
    set(handles.box_yslide,'Value',str2num(C{3}{i}));
    set(handles.box_zslide,'Value',str2num(C{4}{i}));
    set(handles.box_xrot,'Value',str2num(C{5}{i}));
    set(handles.box_yrot,'Value',str2num(C{6}{i}));
    set(handles.box_zrot,'Value',str2num(C{7}{i}));
    set(handles.box_xsize,'String',num2str(str2num(C{8}{i})))
    set(handles.box_ysize,'String',num2str(str2num(C{9}{i})))
    set(handles.box_zsize,'String',num2str(str2num(C{10}{i})))
    set(handles.x_loc,'String',['X Loc: ' num2str(get(handles.box_xslide,'Value'),2)]);
    set(handles.y_loc,'String',['Y Loc: ' num2str(get(handles.box_yslide,'Value'),2)]);
    set(handles.z_loc,'String',['Z Loc: ' num2str(get(handles.box_zslide,'Value'),2)]);
    set(handles.rotx,'String',['X Rot: ' num2str(get(handles.box_xrot,'Value'),2)]);
    set(handles.roty,'String',['Y Rot: ' num2str(get(handles.box_yrot,'Value'),2)]);
    set(handles.rotz,'String',['Z Rot: ' num2str(get(handles.box_zrot,'Value'),2)]);
    
    axes(handles.mip)
    set(gcf, 'Renderer','OpenGL')
    opengl software
    hold on
    new_plane_Callback(hObject, eventdata, handles);
end

function thresh_Callback(hObject, eventdata, handles)
% hObject    handle to thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thresh as text
%        str2double(get(hObject,'String')) returns contents of thresh as a
%        double

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

function interp_num_Callback(hObject, eventdata, handles)
% hObject    handle to interp_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of interp_num as text
%        str2double(get(hObject,'String')) returns contents of interp_num
%        as a double

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

function save_name_Callback(hObject, eventdata, handles)
% hObject    handle to save_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of save_name as text
%        str2double(get(hObject,'String')) returns contents of save_name as a double


% --- Executes during object creation, after setting all properties.

function save_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to save_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in load_2D_slice.
function load_2D_slice_Callback(hObject, eventdata, handles)
% hObject    handle to load_2D_slice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global data pathname PD MAG time_current header_2d OV_2d XV_2d YV_2d;
load_type = get(handles.load_type, 'Value');
if load_type == 1
    [filename_PD, pathname] = uigetfile('*.dat*', 'Select the phase difference .dat file');
elseif load_type == 2
    [filename_PD, pathname] = uigetfile('*.*', 'Select any image from the directory containing your 2D slice of interest');
elseif load_type == 3
    [filename_PD, pathname] = uigetfile('*.mat*', 'Select the phase difference .mat file');
end

filename_check_mat = regexp(filename_PD,'mat','once');
filename_check_dat = regexp(filename_PD, 'dat', 'once');
if load_type == 1
    filename_check_dcm = dicomread([pathname filename_PD]);
end

if filename_PD ~= 0
    if pathname ~= 0
        if ~isempty(filename_check_dat)
            
            % Load header information
            headerDir = fullfile(pathname,'pcvipr_header.txt');
            data = parameter_read(headerDir);
            numFrames = data.nT;
            timeFrames = 1:numFrames;
            time_current = (timeFrames*data.dT)/1000;
            header_2d = struct('ImagePositionPatient', data.s, 'ImageOrientationPatient', [data.i; data.j], 'PixelSpacing', [data.resx; data.resy]); 
            OV_2d = [OV_2d; data.s];
            XV_2d = [XV_2d; data.i];
            YV_2d = [YV_2d; data.j];
            
            % Load magnitude information
            [filename_mag, pathname] = uigetfile('*.dat*', 'Select the magnitude file');     
            fid = fopen([pathname filename_mag],'r');
            MAG = reshape(fread(fid,data.xSize*data.ySize*(data.nT+1),'short'),[data.xSize data.ySize data.nT+1]);
            MAG(:,:,1) = [];
            fclose(fid);
            
            % Load phase difference information
            fid = fopen([pathname filename_PD],'r');
            PD = reshape(fread(fid,data.xSize*data.ySize*(data.nT+1),'short'),[data.xSize data.ySize data.nT+1]);
            PD(:,:,1) = [];
            fclose(fid);  
            
            % Display image
            frame_to_display = round(data.nT/3);
            image = MAG(:,:,frame_to_display);   
            
        elseif ~isempty(filename_check_dcm)
            
            % Load dicom header information
            [Selection, ~] = listdlg('PromptString','Select the magnitude images.','ListString',fileNames,'ListSize',[300 300]);    
            numFrames = numel(Selection);   
            header_2d = dicominfo([pathname filename]);            
            OV_2d = [OV_2d; header_2d.ImagePositionPatient'];
            XV_2d = [XV_2d; header_2d.ImageOrientationPatient(1:3)'*xres];
            YV_2d = [YV_2d; header_2d.ImageOrientationPatient(4:6)'*yres];

            % Load magnitude information
            frame_to_display = round(numFrames/3); % 2010-11-23
            image = dicomread([pathname fileNames{Selection(frame_to_display)}]);
            MAG = zeros([size(image) numFrames]);
            MAG(:,:,1) = image;
            for i = 2:numFrames
                MAG(:,:,i) = dicomread([pathname fileNames{Selection(i)}]);
            end

            % Load phase difference information
            [Selection, ~] = listdlg('PromptString','Select the phase contrast images.','ListString',fileNames,'ListSize',[300 300]);
            if size(Selection) ~= numFrames
                warndlg('You did not select the same number of PC images as magnitude images. Please try again.','!! Warning !!');
            end
            PD = zeros([size(image) numFrames]);
            for i = 1:numFrames
                PD(:,:,i) = dicomread([pathname fileNames{Selection(i)}]);
                info = dicominfo([pathname fileNames{Selection(i)}]); % added 2010-08-24 to get time info
                time_current(i) = info.TriggerTime;
            end  
                % SCALING
                % 2010-12-17 On 11/29/10 Andrew Wentland, Liz Nett,and Mike
                % Loecher scanned a phantom to determine the header fields
                % indicating if the mag mask is off. Here is what we learned:
                % VasCollapse (header.Private_0043_1030), 2 – magnitude, 3 – R/L flow
                % 4 – A/P flow, 5 – S/I flow, 8 – Oblique magnitude, 9 – Oblique R/L flow
                % 10 – Oblique A/P flow, 11 – Oblique S/I flow
                %
                % Venc (header.Private_0019_10cc)
                % VasCollapse (header.Private_0043_1032)
                %       8 – Flow Analysis Off, Mag Mask Off
                %       10 – Flow Analysis Off, Mag Mask On
                %       24 – Flow Analysis On, Mag Mask Off
                %       26 – Flow Analysis On, Mag Mask On
                % VelEncodeScale (header.Private_0019_10e2)
                % You can check VasCollapse to determine if the second bit is
                % turned on; if it is, mag mask is on.
                % 
                % If mag mask is off, we don't need to scale it. If mag mask is
                % on, the phase difference images should be multiplied by VENC,
                % divided by VelEncodeScale, divided by pi, and divided by the
                % magnitude images.
            velEncodeScale = header_2d.Private_0019_10e2;
            venc = double(header_2d.Private_0019_10cc);
            vasCollapse = double(header_2d.Private_0043_1032);
            if bitget(vasCollapse, 2) == 1
                nonZero = find(MAG ~= 0);
                PD(nonZero) = double(PD(nonZero)).*venc./velEncodeScale./pi./double(MAG(nonZero));
            end 
            
        elseif ~isempty(filename_check_mat)
            
            %Load header information
            headerDir = fullfile(pathname,'pcvipr_header.txt');
            data = parameter_read(headerDir);
            numFrames = data.nT;
            timeFrames = 1:numFrames;
            time_current = (timeFrames*data.dT)/1000;
            header_2d = struct('ImagePositionPatient', data.s, 'ImageOrientationPatient', [data.i; data.j], 'PixelSpacing', [data.resx; data.resy]); 
            OV_2d = [OV_2d; data.s];
            XV_2d = [XV_2d; data.i];
            YV_2d = [YV_2d; data.j];
            
            % Load magnitude information
            [filename_mag, pathname] = uigetfile('*.mat*', 'Select the magnitude file');                  
            load(filename_mag,'-mat');
            MAG = eval(filename_mag(1:end-4));
            
            % Load phase difference information
            load(filename_PD,'-mat');
            PD = eval(filename_PD(1:end-4));
            
            % Display image
            frame_to_display = round(data.nT/3);
            image = MAG(:,:,frame_to_display);   
            
        else
            warndlg('An error occurred loading the images.','!! Warning !!');
        end       
        
    % Display the initial image
    axes(handles.slice_axes);
    set(gcf, 'Renderer','OpenGL') %2010-10-19 renderer change
    opengl software  %2010-10-19 renderer change
    imagesc(double(image));
    colormap('gray');
    axis off
    set(handles.slice_axes,'Visible','on');
    set(handles.slice_slider, 'String', num2str(frame_to_display));
    set(handles.slice_slider,'Visible','on');
    set(handles.slice_slider_right,'Visible','on');
    set(handles.slice_slider_left,'Visible','on');
    set(handles.slice_slider_text,'Visible','on');
    set(handles.uipanel6,'Visible','on');
    set(handles.EditSpline,'Visible','on');
    set(handles.slice_slider_text,'String',[num2str(frame_to_display) '/' num2str(numFrames)]);
    set(handles.load_2D_slice,'BackgroundColor',[0.87 0.92 0.98])
    set(handles.load_2D_slice,'ForegroundColor','black')
    set(handles.EditSpline,'BackgroundColor',[0.847 0.1608 0])
    set(handles.EditSpline,'ForegroundColor','white')
    
    end
end
clear filename_check_dat filename_check_dcm filename_check_mat;


function update_2D_slice(handles)

global MAG numFrames;

set(handles.slice_slider_text,'String',[get(handles.slice_slider,'String') '/' num2str(numFrames)]);

% Display the image
axes(handles.slice_axes);
set(gcf, 'Renderer','OpenGL') %2010-10-19 renderer change
opengl software  %2010-10-19 renderer change
imagesc(double(MAG(:,:,str2double(get(handles.slice_slider,'String')))));
colormap('gray');
axis off;

% --- Executes on slider movement.
function slice_slider_Callback(hObject, eventdata, handles)
% hObject    handle to slice_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global numFrames;

current_num = str2double(get(handles.slice_slider,'String'));
if current_num < numFrames+1
    if current_num > 0
        update_2D_slice(handles);
    else
        error('You are outside the range of images.');
    end
else
    error('You are outside the range of images.');
end

% --- Executes during object creation, after setting all properties.
function slice_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slice_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on button press in slice_slider_right.
function slice_slider_right_Callback(hObject, eventdata, handles)
% hObject    handle to slice_slider_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global numFrames;

current_num = str2double(get(handles.slice_slider,'String'));

if current_num+1 > numFrames
    set(handles.slice_slider,'String',num2str(current_num));
    error('You are already at the final image.');
else
    set(handles.slice_slider,'String',num2str(current_num+1));
end
update_2D_slice(handles);

% --- Executes on button press in slice_slider_left.
function slice_slider_left_Callback(hObject, eventdata, handles)
% hObject    handle to slice_slider_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
current_num = str2double(get(handles.slice_slider,'String'));

if current_num-1 < 1
    set(handles.slice_slider,'String',num2str(current_num));
    error('You are already at the first image.');
else
    set(handles.slice_slider,'String',num2str(current_num-1));
end    
update_2D_slice(handles);

% --- Executes on button press in EditSpline.
function EditSpline_Callback(hObject, eventdata, handles)
% hObject    handle to EditSpline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CreateSpline(handles);

set(handles.EditSpline,'BackgroundColor',[0.87 0.92 0.98])
set(handles.EditSpline,'ForegroundColor','black')
set(handles.plotROI,'BackgroundColor',[0.847 0.1608 0])
set(handles.plotROI,'ForegroundColor','white')

function CreateSpline(handles)

global ellipse_h;

axes(handles.slice_axes);
set(gcf, 'Renderer','OpenGL') %2010-10-19 renderer change
opengl software  %2010-10-19 renderer change
ellipse_h = imellipse;

function plotROI(handles)

global PD MAG numFrames ellipse_h;
global time_current flow_current;
global time_2d flow_2d area_2d;
global xres yres;
global header_2d;

flow_current = 0;
BW = createMask(ellipse_h);
[m,n] = size(BW);

% Count the number of voxels within the ROI.
count = 0;
for i = 1:m
    for j = 1:n
        if BW(i,j) ~= 0
            count = count+1;
        end
    end
end

% Scaled area according to the resolution
area_2d = xres*yres*count;

% 2010-09-03
% mag mask is probably on; ask Oliver about getting rid of it; simple
% division?
% mag_mask_off = PD./MAG;
% for k = 1:numFrames
%     for i = 1:m
%         for j = 1:n
%             if isnan(mag_mask_off(i,j,k))
%                 mag_mask_off(i,j,k) = 0;
%             elseif isinf(mag_mask_off(i,j,k))
%                 mag_mask_off(i,j,k) = 0;
%             end
%         end
%     end
% end

for k = 1:numFrames
    
    %multiplied = BW.*mag_mask_off(:,:,k);
    multiplied = BW.*PD(:,:,k);    
    sum = 0;

    for i = 1:m
        for j = 1:n
            sum = sum + multiplied(i,j);
        end
    end
    
    flow_current(k) = sum/count*area_2d/1000; % Divide by 1000 to convert from mm^3 to cm^3

end

% 2010-11-23 Flips the curve if the peak is negative. Doing it this way
% preserves retrograde flow when plotted.
[C,I] = max(abs(flow_current)); clear C;
if flow_current(I(1)) < 0
    flow_current = flow_current.*(-1);
end

set(handles.flow_axes_2d,'Visible','on');
set(handles.add_flow_info,'Visible','on');
set(handles.show_ttf_plot,'Visible','on');
set(handles.spline_interpolate_2d,'Visible','on');
axes(handles.flow_axes_2d);
set(gcf, 'Renderer','OpenGL') %2010-10-19 renderer change
opengl software  %2010-10-19 renderer change

% [x,y] = size(time_2d);
% if time_2d ~= 0
%     for i = 1:y
%         plot(time_2d(:,i),flow_2d(:,i))
%         hold on;
%     end
% end
plot(time_current,flow_current) % 2010-11-23 Removed abs

xlabel('Time (ms)');
ylabel('Flow (ml/s)');

function [ttp, ttf, ttu] = compute_tt_2d(handles, time, flow)

% Time to peak
[c,index] = max(flow);
clear c;
if length(index) > 1
    index = index(1);
end

ttp = time(index);
% End time to peak calculation

% Added 2012-08-01 to correct for when an entire waveform has a DC offset
flow = flow - min(flow(1:index));

% Time to foot
end_point = 0;
counter = 1;
data_points = 0;

for i = index-1:-1:1
    if flow(i) < 0.8*flow(index) % 2010-11-23 removed abs
        end_point(counter) = i;
        counter = counter+1;
    end
end

end_point = end_point(1);
begin_point = 0;

i = end_point-1;
while i > 0 && flow(i) > 0.2*flow(index)
    begin_point = i;
    i = i-1;
end 

% Probably not the best way to code this; i f begin_point is zero, it's
% because the point after end point is below the 20% mark.
if begin_point == 0
    begin_point = end_point-1;
end

data_points = begin_point:1:end_point; clear begin_point; clear end_point;

time_data = zeros(length(data_points),1); flow_data = zeros(length(data_points),1);
for i = 1:length(data_points)
    time_data(i) = time(data_points(i));
    flow_data(i) = flow(data_points(i));
end

p = polyfit(time_data,flow_data,1); % Fits the upstroke linearly, y = mx + b
ttf = -1*p(2)/p(1); % Solves for x in the equation y = mx + b, where y is 0, as we want to find the intersection of the slope of the upstroke with the x-axis
% End time to foot calculation

show_ttf = get(handles.show_ttf_plot,'Value');
if show_ttf
    figure; plot(time,flow); hold on;
    string = [num2str(p(1)) '*x + ' num2str(p(2))];
    fplot(string,[0 time(index)])
    xlim([0 time(end)])
end

% 2nd derivative calculation
%flow_padded = [0 flow 0];
points = fnplt(cscvn([time; flow]));
flow_padded = [0 points(2,:) 0];
second_der = -flow_padded(1:end-2) + 2*flow_padded(2:end-1) -flow_padded(3:end); %/ delt^2;
%figure; plot(time,second_der)
%figure; plot(points(1,:),second_der)

% [c,upstroke] = min(abs(second_der(1:index))); clear c;

[c,index] = max(abs(points(2,:)));
clear c;
if length(index) > 1
    index = index(1);
end

counter = 1;
for i = (index-1):-1:1
    if abs(points(2,i)) > 0.2*abs(points(2,index))
        begin_point(counter) = i;
        counter = counter+1;
    end
end

begin_point = begin_point(end); clear counter;

counter = 1;
upstroke = [];
%for i = 1:index
%for i = data_points(1):index
for i = begin_point:index
    %p = polyfit([time(i) time(i+1)],[second_der(i) second_der(i+1)],1);
    p = polyfit([points(1,i) points(1,i+1)],[second_der(i) second_der(i+1)],1);
    if isinf(p(1)) || isinf(p(2))
    else
        f = @(x)p(1)*x+p(2);
    %x_intercept = fzero(f,time(i));
        x_intercept = fzero(f,points(1,i+1));
%     if x_intercept > time(i)
%         if x_intercept < time(i+1)
%             upstroke(counter) = x_intercept;
%             counter = counter + 1;
%         end
%     end
        if x_intercept > points(1,i)
            if x_intercept < points(1,i+1)
                upstroke(counter) = x_intercept;
                counter = counter + 1;
            end
        end
    end
end
clear counter x_intercept p f;

%ttu = time(upstroke(1)); % ttu = time to upstroke
if ~isempty(upstroke)
    ttu = upstroke(1); % ttu = time to upstroke
else
    [c,upstroke] = min(abs(second_der(1:index))); clear c;
    ttu = time(upstroke(1));
end
clear upstroke index;
% End 2nd derivative calculation

% --- Executes on button press in plotROI.
function plotROI_Callback(hObject, eventdata, handles)
% hObject    handle to plotROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

plotROI(handles);

set(handles.plotROI,'BackgroundColor',[0.87 0.92 0.98])
set(handles.plotROI,'ForegroundColor','black')
set(handles.add_flow_info,'BackgroundColor',[0.847 0.1608 0])
set(handles.add_flow_info,'ForegroundColor','white')

% --- Executes on button press in add_flow_info.
function add_flow_info_Callback(hObject, eventdata, handles)
% hObject    handle to add_flow_info (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global time_current flow_current;
global ttp_2d ttf_2d ttu_2d xcorr_2d;
global xcorr_flow xcorr_time;
global header_2d points_2d normal_2d;

% Need to get a consistent number of points for the time data; alw
% 2010-09-17
time_interpolant = time_current(1):(time_current(end) - time_current(1))/400:time_current(end);
pp = spline(time_current, flow_current, time_interpolant);
xcorr_flow = [xcorr_flow; pp];
xcorr_time = [xcorr_time; time_interpolant];

xcorr_2d = [xcorr_2d; compute_xcorr_2d(handles)];

if get(handles.spline_interpolate_2d,'Value')
    [ttp, ttf, ttu] = compute_tt_2d(handles, xcorr_time(end,:), xcorr_flow(end,:));
else
    [ttp, ttf, ttu] = compute_tt_2d(handles, time_current, flow_current);
end
set(handles.spline_interpolate_2d,'Enable','off');
ttp_2d = [ttp_2d; ttp]
ttf_2d = [ttf_2d; ttf]
ttu_2d = [ttu_2d; ttu]

% Getting the normals and points to define the various 2D slices
xres_2d = header_2d.PixelSpacing(1);
yres_2d = header_2d.PixelSpacing(2);

OV_2d = header_2d.ImagePositionPatient;
XV_2d = header_2d.ImageOrientationPatient(1:3)*xres_2d;
YV_2d = header_2d.ImageOrientationPatient(4:6)*yres_2d;

P1 = [OV_2d(1), OV_2d(2), OV_2d(3)];
P2 = [OV_2d(1)+XV_2d(1), OV_2d(2)+XV_2d(2), OV_2d(3)+XV_2d(3)];
P3 = [OV_2d(1)+YV_2d(1), OV_2d(2)+YV_2d(2), OV_2d(3)+YV_2d(3)];

points_2d = [points_2d; P1];

normal = cross(P1-P2, P1-P3);
normal_2d = [normal_2d; normal];

set(handles.add_flow_info,'BackgroundColor',[0.87 0.92 0.98])
set(handles.add_flow_info,'ForegroundColor','black')
set(handles.load_2D_slice,'BackgroundColor',[0.847 0.1608 0])
set(handles.load_2D_slice,'ForegroundColor','white')

function corr_time = compute_xcorr_2d(handles)

global time_current xcorr_flow;

[x,y] = size(xcorr_flow); clear y;
[c, lags] = xcorr(xcorr_flow(1,:),xcorr_flow(x,:));
[C,I] = max(abs(c)); %added abs
clear C; clear c;

corr_time = max(time_current)/401*abs(lags(I));

% --- Executes on button press in load_candy_cane.
function load_candy_cane_Callback(hObject, eventdata, handles)
% hObject    handle to load_candy_cane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global header_candy;
global xsize_candy;

[filename, pathname] = uigetfile('*.*', 'Select a candy cane view.');

if filename ~= 0
    if pathname ~= 0
        image = dicomread([pathname filename]);
        [xsize_candy,ysize_candy] = size(image);
        clear ysize_candy;
        header_candy = dicominfo([pathname filename]);
        axes(handles.candy_cane);
        set(gcf, 'Renderer','OpenGL') %2010-10-19 renderer change
        opengl software  %2010-10-19 renderer change
        imagesc(image);
        colormap('gray');
        axis off
        
        set(handles.candy_cane,'Visible','on');
        set(handles.place_point,'Visible','on');
        set(handles.load_2d_planes,'Visible','on');
    end
end

set(handles.load_2D_slice,'BackgroundColor',[0.87 0.92 0.98])
set(handles.load_2D_slice,'ForegroundColor','black')
set(handles.place_point,'BackgroundColor',[0.847 0.1608 0])
set(handles.place_point,'ForegroundColor','white')
set(handles.load_2d_planes,'BackgroundColor',[0.847 0.1608 0])
set(handles.load_2d_planes,'ForegroundColor','white')

% --- Executes on button press in place_point.
function place_point_Callback(hObject, eventdata, handles)
% hObject    handle to place_point (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global point_locs centerline_points;

h = impoint(handles.candy_cane);

point_locs = [point_locs; getPosition(h)];

[x,y] = size(point_locs);
clear y;
x_locs = 0; y_locs = 0;
for i = 1:x
    x_locs(i) = point_locs(i,1);
    y_locs(i) = point_locs(i,2);
end
clear x;

hold on;fnplt(cscvn([x_locs; y_locs]))
centerline_points = fnplt(cscvn([x_locs; y_locs])); % This doesn't plot anything, but it gives all
% of the points that would otherwise be plotted with fnplt. I should be
% able to use these points to calculate distance.

% --- Executes on button press in load_2d_planes.
function load_2d_planes_Callback(hObject, eventdata, handles)
% hObject    handle to load_2d_planes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global header_candy points_2d normal_2d;
global OV_candy XV_candy YV_candy;

set(handles.identify_2d_planes,'Visible','on');

% First get an equation for the plane of the candy cane
xres_candy = header_candy.PixelSpacing(1);
yres_candy = header_candy.PixelSpacing(2);

OV_candy = header_candy.ImagePositionPatient;
XV_candy = header_candy.ImageOrientationPatient(1:3)*xres_candy;
YV_candy = header_candy.ImageOrientationPatient(4:6)*yres_candy;

P1 = [OV_candy(1), OV_candy(2), OV_candy(3)];
P2 = [OV_candy(1)+XV_candy(1), OV_candy(2)+XV_candy(2), OV_candy(3)+XV_candy(3)];
P3 = [OV_candy(1)+YV_candy(1), OV_candy(2)+YV_candy(2), OV_candy(3)+YV_candy(3)];
P_candy = P1;

normal_candy = cross(P1-P2, P1-P3);

[a,b] = size(points_2d); clear b;

% See:
% http://www.mathworks.com/matlabcentral/fileexchange/17618-plane-intersection
% Thank you, Kevin Johnson, for helping me to figure this out! 

axes(handles.candy_cane)
set(gcf, 'Renderer','OpenGL') %2010-10-19 renderer change
opengl software  %2010-10-19 renderer change

for i = 1:a
    [P,N,check] = plane_intersect(normal_candy,P_candy,normal_2d(i,:),points_2d(i,:)); clear check;

    A = [XV_candy YV_candy];
    point1 = linsolve(A,(P-N*500)' - OV_candy);
    point2 = linsolve(A,(P+N*500)' - OV_candy);

    line([point1(1) point2(1)], [point1(2) point2(2)])
end

set(handles.load_2d_planes,'BackgroundColor',[0.87 0.92 0.98])
set(handles.load_2d_planes,'ForegroundColor','black')
set(handles.identify_2d_planes,'BackgroundColor',[0.847 0.1608 0])
set(handles.identify_2d_planes,'ForegroundColor','white')

% --- Executes on button press in identify_2d_planes.
function identify_2d_planes_Callback(hObject, eventdata, handles)
% hObject    handle to identify_2d_planes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global plane_locs centerline_points distances_2d;
global ttp_2d;
global header_candy;

xres_candy = header_candy.PixelSpacing(1);

if centerline_points
    h = impoint(handles.candy_cane);

    plane_locs = [plane_locs; getPosition(h)];

    [x,y] = size(centerline_points); clear x;

    [o,p] = size(plane_locs); clear p;

    closest_points = zeros(1,o);
    for i = 1:o
        for j = 1:y
            x_diff(j) = abs(plane_locs(i,1) - centerline_points(1,j));
            y_diff(j) = abs(plane_locs(i,2) - centerline_points(2,j));
            dist(j) = hypot(x_diff(j),y_diff(j));
        end
    
        clear x_diff y_diff;
        [C,I] = min(dist); clear C;
        closest_points(i) = I(1); clear I;
    end

    distances_2d = 0;
    
    % Compute the distance between planes
    if length(closest_points) > 1
   
        for i = 2:o
       
            dist = 0;
       
            for j = closest_points(1):closest_points(i)-1
                x_diff = abs(centerline_points(1,j) - centerline_points(1,j+1));
                y_diff = abs(centerline_points(2,j) - centerline_points(2,j+1));
                dist = dist + hypot(x_diff,y_diff);
            end
       
            distances_2d(i) = dist;% removed the -1 in distances_2d(i-1) so that the first point would be at 0 mm; 2010-09-03
            distances_2d = distances_2d.*xres_candy; % 2010-05-10; This converts the pixel distance of distances_2d into distance in cm.
            % This assumes, however, that the voxels are isotropic, as it only multiplies by the x-resolution.
        end
    end
    
else
    warndlg('You have not yet defined a centerline.','!! Warning !!');
end

a = length(ttp_2d);
b = length(distances_2d);
if b == a
    set(handles.identify_2d_planes,'BackgroundColor',[0.87 0.92 0.98])
    set(handles.identify_2d_planes,'ForegroundColor','black')
    set(handles.identify_2d_planes,'Enable','off');
    
    set(handles.show_PWV,'Visible','on');
    set(handles.show_PWV,'BackgroundColor',[0.847 0.1608 0])
    set(handles.show_PWV,'ForegroundColor','white')
end

clear a b;

set(handles.place_point,'BackgroundColor',[0.87 0.92 0.98])
set(handles.place_point,'ForegroundColor','black')

% --- Calculate the distance along the centerline on the 3D plot.
function compute_3d_distance()

global centerline_points_3d distance_3d;
global plane_counter plane_coordinates;
global flow_array;

% Compute the distance along the centerline in the 3D plot
if centerline_points_3d

    distance_3d = 0;
    [x,y] = size(centerline_points_3d); clear x; % Number of points along the centerline
    
    numPlanes = size(flow_array,1)-1; % check this
    for i = 2:numPlanes % Number of planes, starting from plane #2
        
        diff = zeros(y);
        
        for j = 1:y
            diff(j) = sqrt((plane_coordinates(i,1)-centerline_points_3d(1,j))^2 + (plane_coordinates(i,2)-centerline_points_3d(2,j))^2 + (plane_coordinates(i,3)-centerline_points_3d(3,j))^2);
        end
        
        [C,I] = min(diff); I = I(1); clear C;
        
        total = 0;
        
        for k = 1:I-1
            x_diff = abs(centerline_points_3d(1,k) - centerline_points_3d(1,k+1));
            y_diff = abs(centerline_points_3d(2,k) - centerline_points_3d(2,k+1));
            z_diff = abs(centerline_points_3d(3,k) - centerline_points_3d(3,k+1));
            
            total = total + sqrt(x_diff^2 + y_diff^2 + z_diff^2);
        end
        
        distance_3d(i) = total;
    end
end

% --- Executes on button press in show_PWV.
function show_PWV_Callback(hObject, eventdata, handles)
% hObject    handle to show_PWV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global distances_2d;
global ttp_2d ttf_2d ttu_2d xcorr_2d;

set(handles.show_ttf_plot,'Visible','off')
set(handles.spline_interpolate_2d,'Visible','off')

% Plot the Time to Upstroke data
set(handles.ttu_2d_plot,'Visible','on');
axes(handles.ttu_2d_plot);
set(gcf, 'Renderer','OpenGL') %2010-10-19 renderer change
opengl software  %2010-10-19 renderer change
plot(distances_2d,ttu_2d); ylabel('ttu [ms]');
hold on;
p = polyfit(distances_2d,ttu_2d',1); % Find the slope and intercept of a best fit line to the data
f = @(x) p(1)*x + p(2); % Define an anonymous function, y = mx + b
fplot(f,[0 max(distances_2d)],':r')
set(handles.ttu_2d_pwv,'Visible','on');
set(handles.ttu_2d_pwv,'String',[num2str(round2(1/p(1),0.1)) ' m/s']); % Display PWV

% Plot the Time to Peak data
set(handles.ttp_2d_plot,'Visible','on');
axes(handles.ttp_2d_plot);
set(gcf, 'Renderer','OpenGL') %2010-10-19 renderer change
opengl software  %2010-10-19 renderer change
plot(distances_2d,ttp_2d); ylabel('ttp [ms]');
hold on;
p = polyfit(distances_2d,ttp_2d',1); % Find the slope and intercept of a best fit line to the data
f = @(x) p(1)*x + p(2); % Define an anonymous function, y = mx + b
fplot(f,[0 max(distances_2d)],':r')
set(handles.ttp_2d_pwv,'Visible','on');
set(handles.ttp_2d_pwv,'String',[num2str(round2(1/p(1),0.1)) ' m/s']); % Display PWV

% Plot the Time to Foot data
set(handles.ttf_2d_plot,'Visible','on');
axes(handles.ttf_2d_plot);
set(gcf, 'Renderer','OpenGL') %2010-10-19 renderer change
opengl software  %2010-10-19 renderer change
plot(distances_2d,ttf_2d); ylabel('ttf [ms]');
hold on;
p = polyfit(distances_2d,ttf_2d',1); % Find the slope and intercept of a best fit line to the data
f = @(x) p(1)*x + p(2); % Define an anonymous function, y = mx + b
fplot(f,[0 max(distances_2d)],':r')
set(handles.ttf_2d_pwv,'Visible','on');
set(handles.ttf_2d_pwv,'String',[num2str(round2(1/p(1),0.1)) ' m/s']); % Display PWV

% Plot the Cross Correlation data
set(handles.xcorr_2d_plot,'Visible','on');
axes(handles.xcorr_2d_plot);
set(gcf, 'Renderer','OpenGL') %2010-10-19 renderer change
opengl software  %2010-10-19 renderer change
plot(distances_2d,xcorr_2d); ylabel('xcorr [ms]'); xlabel('distance [mm]');
hold on;
p = polyfit(distances_2d,xcorr_2d',1); % Find the slope and intercept of a best fit line to the data
f = @(x) p(1)*x + p(2); % Define an anonymous function, y = mx + b
fplot(f,[0 max(distances_2d)],':r')
set(handles.xcorr_2d_pwv,'Visible','on');
set(handles.xcorr_2d_pwv,'String',[num2str(round2(1/p(1),0.1)) ' m/s']); % Display PWV

set(handles.show_PWV,'BackgroundColor',[0.87 0.92 0.98])
set(handles.show_PWV,'ForegroundColor','black')

% --- Executes on button press in show_PWV_3d.
function show_PWV_3d_Callback(hObject, eventdata, handles)
% hObject    handle to show_PWV_3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global flow_t flow_array distance_3d distances_2d;
global plane_counter;
global pwv_3d_computed;

num2cut = str2double(get(handles.num2cut,'String'));

flow_t = flow_t(1:end-num2cut);
flow_array = flow_array(:,1:end-num2cut);
% for i = 1:size(flow_array,1)
%     flow_array = flow_array(i,1:end-num2cut);

pwv_3d_computed = 1;

set(handles.new_plane,'Enable','off');
set(handles.plane_select,'Enable','off');

set(handles.box_xslide,'Enable','off');
set(handles.box_yslide,'Enable','off');
set(handles.box_zslide,'Enable','off');

set(handles.box_xrot,'Enable','off');
set(handles.box_yrot,'Enable','off');
set(handles.box_zrot,'Enable','off');

set(handles.box_xsize,'Enable','off');
set(handles.box_ysize,'Enable','off');
set(handles.box_zsize,'Enable','off');

flow_array = [flow_array; flow_t];

if get(handles.use_2d_distance, 'Value') == 1
    if not(isempty(distances_2d))
        dist = distances_2d;
    else
        warndlg('ERROR: No distances and/or centerline were defined from the 2D images!!','!! Warning !!');
    end
else
    compute_3d_distance();
    dist = distance_3d;
end
    
num_planes = size(flow_array,1) - 1; % Subtract 1 since the first row is time information

% 2011-02-03
% This is code to use the check box that allows only the first and last
% slices to be used for the PWV calculations. Saving the plane information
% should include all planes, so it should be simple to reload the data and
% compute PWV with or without all slices.
if get(handles.first_last,'Value')
    if num_planes > 2
        num_planes = 2;
        flow_array = [flow_array(1:2,:); flow_array(end,:)];
        dist
        dist = [dist(1) dist(end)];
    end
end

% 2011-02-04
% This is code to use the check box that allows for a spline interpolation
% of the flow waveforms.
if get(handles.spline_interpolate,'Value')
    time_interpolant = flow_array(1,1):(flow_array(1,end) - flow_array(1,1))/400:flow_array(1,end);
    
    for i = 1:(size(flow_array,1)-1)
    	flow_array_temp(i,:) = spline(flow_array(1,:), flow_array(i+1,:), time_interpolant);
        figure; plot(flow_array(1,:),flow_array(i+1,:)); hold on; plot(time_interpolant,flow_array_temp(i,:))
    end
    
    flow_array = [time_interpolant; flow_array_temp]; clear time_interpolant flow_array_temp;
end

for i = 1:num_planes
    
    % Begin time-to-peak (ttp) calculation
    [c,index] = max(abs(flow_array(i+1,:))); clear c;
    index = index(1);
    
    ttp_3d(i) = flow_array(1,index);
    % End time-to-peak calculation

    % Begin time-to-foot (ttf) calculation
    end_point = 0;
    counter = 1;
    data_points = 0;
    
    minimum = min(flow_array(i+1,1:index));
    minimum = minimum(1);

    for m = index-1:-1:1
        if abs(flow_array(i+1,m)-minimum) < 0.8*abs(flow_array(i+1,index)-minimum)
            end_point(counter) = m;
            counter = counter+1;
        end
    end

    end_point = end_point(1);
    
    begin_point = 0;
    %counter = 1;
    
    n = end_point-1;
    while n > 0 && abs(flow_array(i+1,n)-minimum) > 0.2*abs(flow_array(i+1,index)-minimum)
        begin_point = n;
        n = n-1;
    end
    
    if begin_point == 0
        %begin_point = 1;
        if end_point > 1
            begin_point = end_point-1; % If the data point after end_point already falls below the 20% threshold, just set the beginning point of
        else
            begin_point = 1;
        end
        % the upstroke to be the next point
    else
        begin_point = begin_point(1);
    end
    clear counter;

    data_points = begin_point:1:end_point; clear begin_point; clear end_point;

    time_data = zeros(1,length(data_points)); flow_data = zeros(1,length(data_points));
    for o = 1:length(data_points)
        time_data(o) = flow_array(1,data_points(o));
        flow_data(o) = flow_array(i+1,data_points(o))-minimum;
%         flow_data(o) = flow_array(i+1,data_points(o))-flow_array(i+1,data_points(1));
    end

    p = polyfit(time_data,abs(flow_data),1);
    ttf_3d(i) = -1*p(2)/p(1);
    
    figure; plot(flow_array(1,:),abs(flow_array(i+1,:)-minimum)); hold on;
%     figure; plot(flow_array(1,:),abs(flow_array(i+1,:)-flow_array(i+1,data_points(1)))); hold on;
    string = [num2str(p(1)) '*x + ' num2str(p(2))];
    fplot(string,[0 flow_array(1,index)])
    set(gca,'XLim',[0 flow_array(1,end)]);
    
    %clear index; 2010-11-07
    % End time-to-foot calculation
    
    % 2nd derivative calculation
    points = fnplt(cscvn([flow_array(1,:); flow_array(i+1,:)-minimum]));
    flow_padded = [0 points(2,:) 0];
    second_der = -flow_padded(1:end-2) + 2*flow_padded(2:end-1) -flow_padded(3:end); %/ delt^2;
    %figure; plot(points(1,:),second_der)

    [c,index] = max(abs(points(2,:)));
    clear c;
    if length(index) > 1
        index = index(1);
    end

    counter = 1;
    for q = (index-1):-1:1
        if abs(points(2,q)) > 0.2*abs(points(2,index))
            begin_point(counter) = q;
            counter = counter+1;
        end
    end

    begin_point = begin_point(end); clear counter;

    counter = 1;
    upstroke = [];
    
    for r = begin_point:index
        p = polyfit([points(1,r) points(1,r+1)],[second_der(r) second_der(r+1)],1);
        if isinf(p(1)) || isinf(p(2))
        else
            f = @(x)p(1)*x+p(2);
            x_intercept = fzero(f,points(1,r+1));

            if x_intercept > points(1,r)
                if x_intercept < points(1,r+1)
                    upstroke(counter) = x_intercept;
                    counter = counter + 1;
                end
            end
        end
    end
    clear counter x_intercept p f;

    if ~isempty(upstroke)
        ttu_3d(i) = upstroke(1); % ttu = time to upstroke
    else
        [c,upstroke] = min(abs(second_der(1:index))); clear c;
        ttu_3d(i) = time(upstroke(1));
    end
    clear upstroke index;
    % End 2nd derivative calculation
end

% The following commented-out code is for computing xcorr with interpolated
% flow waveforms.
% time = flow_array(1,:);
% time_interpolant = time(1):(time(end) - time(1))/400:time(end);
% pp = spline(time, flow_array(2,:), time_interpolant);
% 
% % Begin cross-correlation (xcorr) calculation
% for k = 2:plane_counter
%     pp2 = spline(time,flow_array(k+1,:),time_interpolant);
%     
%     [c, lags] = xcorr(pp,pp2);
%     [C,I] = max(abs(c)); clear C; clear c;
%     xcorr_3d(k) = time(end)/401*abs(lags(I));
% end
% End cross-correlation calculation

% Begin cross-correlation (xcorr) calculation
% This is computed without interpolation of the flow waveforms
for k = 2:(size(flow_array,1)-1)
    [c, lags] = xcorr(flow_array(2,:)-minimum,flow_array(k+1,:)-minimum);
    [C,I] = max(abs(c)); clear C; clear c;
    xcorr_3d(k) = (flow_array(1,2)-flow_array(1,1))*abs(lags(I));
end
% End cross-correlation calculation

% Plot the Time to Upstroke data
set(handles.ttu_3d_plot,'Visible','on');
axes(handles.ttu_3d_plot);
set(gcf, 'Renderer','OpenGL')
opengl software
plot(dist,ttu_3d); ylabel('ttu [ms]');
hold on;
p = polyfit(dist,ttu_3d,1); % Find the slope and intercept of a best fit line to the data
f = @(x) p(1)*x + p(2); % Define an anonymous function, y = mx + b
fplot(f,[0 max(dist)],':r')
set(handles.ttu_3d_pwv,'Visible','on');
set(handles.ttu_3d_pwv,'String',[num2str(round2(1/p(1),0.1)) ' m/s']); % Display PWV

% Plot the Time to Peak data
set(handles.ttp_3d_plot,'Visible','on');
axes(handles.ttp_3d_plot);
set(gcf, 'Renderer','OpenGL') %2010-10-19 renderer change
opengl software  %2010-10-19 renderer change
plot(dist,ttp_3d); ylabel('ttp [ms]');
hold on;
p = polyfit(dist,ttp_3d,1); % Find the slope and intercept of a best fit line to the data
f = @(x) p(1)*x + p(2); % Define an anonymous function, y = mx + b
fplot(f,[0 max(dist)],':r')
set(handles.ttp_3d_pwv,'Visible','on');
set(handles.ttp_3d_pwv,'String',[num2str(round2(1/p(1),0.1)) ' m/s']); % Display PWV

% Plot the Time to Foot data
set(handles.ttf_3d_plot,'Visible','on');
axes(handles.ttf_3d_plot);
set(gcf, 'Renderer','OpenGL') %2010-10-19 renderer change
opengl software  %2010-10-19 renderer change
plot(dist,ttf_3d); ylabel('ttf [ms]');
hold on;
p = polyfit(dist,ttf_3d,1); % Find the slope and intercept of a best fit line to the data
f = @(x) p(1)*x + p(2); % Define an anonymous function, y = mx + b
fplot(f,[0 max(dist)],':r')
set(handles.ttf_3d_pwv,'Visible','on');
set(handles.ttf_3d_pwv,'String',[num2str(round2(1/p(1),0.1)) ' m/s']); % Display PWV

% Plot the Cross Correlation data
set(handles.xcorr_3d_plot,'Visible','on');
axes(handles.xcorr_3d_plot);
set(gcf, 'Renderer','OpenGL') %2010-10-19 renderer change
opengl software  %2010-10-19 renderer change
plot(dist,xcorr_3d); ylabel('xcorr [ms]'); xlabel('distance [mm]');
hold on;
p = polyfit(dist,xcorr_3d,1); % Find the slope and intercept of a best fit line to the data
f = @(x) p(1)*x + p(2); % Define an anonymous function, y = mx + b
fplot(f,[0 max(dist)],':r')
set(handles.xcorr_3d_pwv,'Visible','on');
set(handles.xcorr_3d_pwv,'String',[num2str(round2(1/p(1),0.1)) ' m/s']); % Display PWV

ttu_3d
ttp_3d
ttf_3d
xcorr_3d

% --- Executes when the Update Image button is pressed.
function data = parameter_read(headerDir)

fid = fopen(headerDir);
C = textscan(fid,'%s %s'); % Reads PCVIPR header
parameter = C{1};
value = C{2};
fclose(fid);

global U_pcvipr O_pcvipr;

% grab all header information using findVal function.
data.fovx = findVal(parameter,value,'fovx',1);
data.fovy = findVal(parameter,value,'fovy',1);
data.fovz = findVal(parameter,value,'fovz',1);
data.xSize = findVal(parameter,value,'matrixx',1);
data.ySize = findVal(parameter,value,'matrixy',1);
data.zSize = findVal(parameter,value,'matrixz',1);
data.nT = findVal(parameter,value,'frames',1);
data.dT = findVal(parameter,value,'timeres',1);
data.VENC = findVal(parameter,value,'VENC',1);
data.version = findVal(parameter,value,'version',1);
data.resx = data.fovx./data.xSize;
data.resy = data.fovy./data.ySize;
data.resz = data.fovz./data.zSize;
data.i = findVal(parameter,value,'ix',3);
data.j = findVal(parameter,value,'jx',3);
data.k = findVal(parameter,value,'kx',3);
data.s = findVal(parameter,value,'sx',3);

if data.version > 1
    %%PC VIPR Orientation
    U_pcvipr = zeros(3,3);
    U_pcvipr(1,:) = data.i;
    U_pcvipr(2,:) = data.j;
    U_pcvipr(3,:) = data.k;
    O_pcvipr      = data.s;
else
   warndlg('ERROR: VERSION OUT OF DATE. Please recon with latest recon to get header.');
end


% --- Executes on button press in start_over_3d.
function start_over_3d_Callback(hObject, eventdata, handles)
% hObject    handle to show_PWV_3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global U_pcvipr O_pcvipr;
global res_pcvipr;

global flow_t flow_array distance_3d time;
global plane_counter plane_info;

global vis_thresh;
global verts norms norm_handle hpatch;

global plane_handle check_counter plane_color;
global overall_position plane_coordinates;
global centerline_points_3d;

global slice_CD slice_VEL slice_MASK;
global thresh;

clear U_pcvipr O_pcvipr res_pcvipr;

clear flow_t flow_array distance_3d time;
clear plane_info; plane_counter = 0;
check_counter = 1;

clear vis_thresh;
clear verts norms norm_handle hpatch;

plane_handle = [];
clear plane_color;
clear overall_position plane_coordinates;
clear centerline_points_3d;

clear slice_CD slice_VEL slice_MASK;
clear thresh;

set(handles.new_plane,'Enable','on');
set(handles.plane_select,'Enable','on');

axes(handles.mip); cla;
axes(handles.cd_axes); cla reset;
axes(handles.flow_axes); cla reset;
axes(handles.mean_flow); cla;
axes(handles.ttu_3d_plot); cla;
axes(handles.ttp_3d_plot); cla;
axes(handles.ttf_3d_plot); cla;
axes(handles.xcorr_3d_plot); cla;

set(handles.ttu_3d_plot,'Visible','off');
set(handles.ttp_3d_plot,'Visible','off');
set(handles.ttf_3d_plot,'Visible','off');
set(handles.xcorr_3d_plot,'Visible','off');
set(handles.ttu_3d_pwv,'Visible','off');
set(handles.ttp_3d_pwv,'Visible','off');
set(handles.ttf_3d_pwv,'Visible','off');
set(handles.xcorr_3d_pwv,'Visible','off');

set(handles.use_2d_distance,'Value',0);

set(handles.box_xslide,'Value',0.5);
set(handles.box_yslide,'Value',0.5);
set(handles.box_zslide,'Value',0.5);

set(handles.box_xrot,'Value',0.5);
set(handles.box_yrot,'Value',0.5);
set(handles.box_zrot,'Value',0.5);

set(handles.x_loc,'String','X Loc: 0.5');
set(handles.y_loc,'String','Y Loc: 0.5');
set(handles.z_loc,'String','Z Loc: 0.5');

set(handles.rotx,'String','Rot X: 0.5');
set(handles.roty,'String','Rot Y: 0.5');
set(handles.rotz,'String','Rot Z: 0.5');

set(handles.box_xslide,'Enable','on');
set(handles.box_yslide,'Enable','on');
set(handles.box_zslide,'Enable','on');

set(handles.box_xrot,'Enable','on');
set(handles.box_yrot,'Enable','on');
set(handles.box_zrot,'Enable','on');

set(handles.box_xsize,'Enable','on');
set(handles.box_ysize,'Enable','on');
set(handles.box_zsize,'Enable','on');

set(handles.box_xsize,'String','3');
set(handles.box_ysize,'String','3');
set(handles.box_zsize,'String','3');

set(handles.lumen_threshold,'String','0.2');
set(handles.thresh,'String','40');
set(handles.interp_num,'String','2');

set(handles.plane_select,'String','01')

% --- Executes when the Update Image button is pressed.
function correlate(hObject, eventdata, handles)

global OV_candy XV_candy YV_candy;
global OV_2d XV_2d YV_2d;
global plane_locs;
global true_coordinates; % true_coordinates are the scanner x,y,z points for
% where the 2D slices intersect the aorta based off of the candy cane view

global U_pcvipr O_pcvipr;
global res_pcvipr;

global m_xstart m_xstop;
global m_ystart m_ystop;
global m_zstart m_zstop;

global plane_counter; plane_counter = 0;

[x,y] = size(plane_locs); clear y;

vipr_start(1) = O_pcvipr(1) + m_xstart*U_pcvipr(1,1) + m_ystart*U_pcvipr(1,2) + m_zstart*U_pcvipr(1,3);
vipr_start(2) = O_pcvipr(2) + m_xstart*U_pcvipr(2,1) + m_ystart*U_pcvipr(2,2) + m_zstart*U_pcvipr(2,3);
vipr_start(3) = O_pcvipr(3) + m_xstart*U_pcvipr(3,1) + m_ystart*U_pcvipr(3,2) + m_zstart*U_pcvipr(3,3);

% Not really needed; in here for bebugging purposes, as
% true_coordinates should fall within the range between vipr_start and
% vipr_end.
% vipr_end(1) = vipr_start(1) + (m_xstop - m_xstart)*U_pcvipr(1,1) + (m_ystop - m_ystart)*U_pcvipr(1,2) + (m_zstop - m_zstart)*U_pcvipr(1,3);
% vipr_end(2) = vipr_start(2) + (m_xstop - m_xstart)*U_pcvipr(2,1) +
% (m_ystop - m_ystart)*U_pcvipr(2,2) + (m_zstop - m_zstart)*U_pcvipr(2,3);
% vipr_end(3) = vipr_start(3) + (m_xstop - m_xstart)*U_pcvipr(3,1) + (m_ystop - m_ystart)*U_pcvipr(3,2) + (m_zstop - m_zstart)*U_pcvipr(3,3);

for i = 1:x
    % Scanner coordinates for the points of the 2D planes intersecting the
    % centerline of the candy cane view of the aorta.
    % I'm fairly confident that true_coordinates is correct.
    true_coordinates(i,1) = OV_candy(1) + XV_candy(1)*plane_locs(i,1) + YV_candy(1)*plane_locs(i,2);
    true_coordinates(i,2) = OV_candy(2) + XV_candy(2)*plane_locs(i,1) + YV_candy(2)*plane_locs(i,2);
    true_coordinates(i,3) = OV_candy(3) + XV_candy(3)*plane_locs(i,1) + YV_candy(3)*plane_locs(i,2);
    
    % Gives a different coordinate on each 2D slice by adding 5 units in each direction from
    % XV_2d and YV_2d to the points identified from the intersection of the 2D slices with the
    % candy cane view of the aorta.
    diff_2d(i,1) = true_coordinates(i,1) + 5*XV_2d(i,1) + 5*YV_2d(i,1);
    diff_2d(i,2) = true_coordinates(i,2) + 5*XV_2d(i,2) + 5*YV_2d(i,2);
    diff_2d(i,3) = true_coordinates(i,3) + 5*XV_2d(i,3) + 5*YV_2d(i,3);
    
    % VIPR coordinates for the points of the 2D planes intersecting the
    % centerline of the candy cane view of the aorta.
    vipr_coordinates(i,1) = abs((true_coordinates(i,1) - vipr_start(1))/(U_pcvipr(1,1)+U_pcvipr(1,2)+U_pcvipr(1,3)));
    vipr_coordinates(i,2) = abs((true_coordinates(i,2) - vipr_start(2))/(U_pcvipr(2,1)+U_pcvipr(2,2)+U_pcvipr(2,3)));
    vipr_coordinates(i,3) = abs((true_coordinates(i,3) - vipr_start(3))/(U_pcvipr(3,1)+U_pcvipr(3,2)+U_pcvipr(3,3)));

    % The fraction of the coordinates along the length of the segmented
    % VIPR data set. This is used for setting box_(x,y,z)slide.
    vipr_fraction(i,1) = vipr_coordinates(i,1)/(m_xstop - m_xstart);
    vipr_fraction(i,2) = vipr_coordinates(i,2)/(m_ystop - m_ystart);
    vipr_fraction(i,3) = vipr_coordinates(i,3)/(m_zstop - m_zstart);
    
    % This gives the VIPR coordinates for the second set of coordinates
    % above.
    vipr_coordinates_diff(i,1) = abs((diff_2d(i,1) - vipr_start(1))/(U_pcvipr(1,1)+U_pcvipr(1,2)+U_pcvipr(1,3)));
    vipr_coordinates_diff(i,2) = abs((diff_2d(i,2) - vipr_start(2))/(U_pcvipr(2,1)+U_pcvipr(2,2)+U_pcvipr(2,3)));
    vipr_coordinates_diff(i,3) = abs((diff_2d(i,3) - vipr_start(3))/(U_pcvipr(3,1)+U_pcvipr(3,2)+U_pcvipr(3,3)));
    
    % Gives the angle that the 2d slices make in each direction. These
    % angles are then used below to compute the slope for the planes
    % created in the VIPR surface rendering.
    angle(i,1) = atan((vipr_coordinates_diff(i,3)-vipr_coordinates(i,3))/(vipr_coordinates_diff(i,2)-vipr_coordinates(i,2)))/pi*180;
    angle(i,2) = atan((vipr_coordinates_diff(i,3)-vipr_coordinates(i,3))/(vipr_coordinates_diff(i,1)-vipr_coordinates(i,1)))/pi*180;
    angle(i,3) = atan((vipr_coordinates_diff(i,2)-vipr_coordinates(i,2))/(vipr_coordinates_diff(i,1)-vipr_coordinates(i,1)))/pi*180;
    
    for a = 1:3
        if angle(i,a) < 0
            angle(i,a) = angle(i,a) + 360;
        end
    end
    clear a;
    
    % The rotation slider bars on the GUI begin from 90 degrees. Therefore,
    % below, if the angle is in the first quadrant, subtract from 90, if
    % the angle is in the second quadrant, subtract 180, take the absolute
    % value, and add 90 to get it to be the angle from 90 degrees, and so
    % on.
    for a = 1:3
        if angle(i,a) < 90
            slope(i,a) = (90 - angle(i,a))/180;
        elseif angle(i,a) < 180
            slope(i,a) = (abs(angle(i,a) - 180) + 90)/180;
        elseif angle(i,a) < 270
            slope(i,a) = (270 - angle(i,a))/180;
        else
            slope(i,a) = (450 - angle(i,a))/180;
        end
    end
    
    set(handles.box_xslide,'Value',vipr_fraction(i,2));% I think x and y are backwards above in update_box, so the sliders need to be set with x and y swapped
    set(handles.box_yslide,'Value',vipr_fraction(i,1));
    set(handles.box_zslide,'Value',vipr_fraction(i,3));
    set(handles.box_xrot,'Value',slope(i,1));
    set(handles.box_yrot,'Value',slope(i,2));
    set(handles.box_zrot,'Value',slope(i,3));
    set(handles.box_xsize,'String','8')
    set(handles.box_ysize,'String','8')
    set(handles.box_zsize,'String','1')
    set(handles.x_loc,'String',['X Loc: ' num2str(get(handles.box_xslide,'Value'),2)]);
    set(handles.y_loc,'String',['Y Loc: ' num2str(get(handles.box_yslide,'Value'),2)]);
    set(handles.z_loc,'String',['Z Loc: ' num2str(get(handles.box_zslide,'Value'),2)]);
    set(handles.rotx,'String',['X Rot: ' num2str(get(handles.box_xrot,'Value'),2)]);
    set(handles.roty,'String',['Y Rot: ' num2str(get(handles.box_yrot,'Value'),2)]);
    set(handles.rotz,'String',['Z Rot: ' num2str(get(handles.box_zrot,'Value'),2)]);
    
    axes(handles.mip)
    set(gcf, 'Renderer','OpenGL') %2010-10-19 renderer change
    opengl software  %2010-10-19 renderer change
    hold on
    new_plane_Callback(hObject, eventdata, handles);
end


% --- Executes on button press in 2D --> 3D.
function twoD_3D_compare_Callback(hObject, eventdata, handles)
% hObject    handle to show_PWV_3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

correlate(hObject, eventdata, handles);

% --- Executes on the drop-down menu, plane_select.
function plane_select_Callback(hObject, eventdata, handles)
% hObject    handle to show_PWV_3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global plane_handle plane_counter;

plane_counter = get(handles.plane_select,'Value');

set(handles.box_xslide,'Value',plane_handle(plane_counter).xloc)
set(handles.box_yslide,'Value',plane_handle(plane_counter).yloc)
set(handles.box_zslide,'Value',plane_handle(plane_counter).zloc)
set(handles.box_xrot,'Value',plane_handle(plane_counter).xrot)
set(handles.box_yrot,'Value',plane_handle(plane_counter).yrot)
set(handles.box_zrot,'Value',plane_handle(plane_counter).zrot)
set(handles.box_xsize,'String',num2str(plane_handle(plane_counter).xsize))
set(handles.box_ysize,'String',num2str(plane_handle(plane_counter).ysize))
set(handles.box_zsize,'String',num2str(plane_handle(plane_counter).zsize))
set(handles.x_loc,'String',['X Loc: ' num2str(plane_handle(plane_counter).xloc,2)]);
set(handles.y_loc,'String',['Y Loc: ' num2str(plane_handle(plane_counter).yloc,2)]);
set(handles.z_loc,'String',['Z Loc: ' num2str(plane_handle(plane_counter).zloc,2)]);
set(handles.rotx,'String',['X Rot: ' num2str(plane_handle(plane_counter).xrot,2)]);
set(handles.roty,'String',['Y Rot: ' num2str(plane_handle(plane_counter).yrot,2)]);
set(handles.rotz,'String',['Z Rot: ' num2str(plane_handle(plane_counter).zrot,2)]);

axes(handles.mip)
delete(plane_handle(plane_counter).handle)

update_box(handles)

% --- Executes on the button press add_centerline.
function add_centerline_Callback(hObject, eventdata, handles)
% hObject    handle to show_PWV_3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global plane_coordinates centerline_points_3d;

plane_coordinates

o = size(plane_coordinates,1);
% Draw a centerline that connects the centroids for each plane. This is a
% cubic spline interpolation. -alw
if o > 1
    curve = cscvn(plane_coordinates');
    %curve.coefs
    axes(handles.mip);
    set(gcf, 'Renderer','OpenGL') %2010-10-19 renderer change
    opengl software  %2010-10-19 renderer change
    fnplt(curve,1,'y'); % 1 is the thickness of the plotted line; y is the color yellow for the line
    centerline_points_3d = fnplt(cscvn(plane_coordinates'));
    % This doesn't plot anything, but it gives all
    % of the points that would otherwise be plotted with fnplt. I should be
    % able to use these points to calculate distance.
end
clear o;


% --- Executes on button press in start_over_2d.
function start_over_2d_Callback(hObject, eventdata, handles)
% hObject    handle to show_PWV_3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global PD MAG numFrames;
global xres yres;
global time_current flow_current;
% global header_2d;
global OV_2d XV_2d YV_2d;

global ellipse_h;

global time_2d flow_2d;
global ttp_2d ttf_2d ttu_2d xcorr_2d;

global xcorr_flow;
global header_2d points_2d normal_2d;

global header_candy xsize_candy;
global point_locs centerline_points;

global OV_candy XV_candy YV_candy;
global distances_2d true_coordinates;

clear PD MAG numFrames;
clear xres yres;
clear time_current flow_current;
clear header_2d;
clear OV_2d XV_2d YV_2d;

clear ellipse_h;

clear time_2d flow_2d;
clear ttp_2d ttf_2d ttu_2d xcorr_2d;

clear xcorr_flow;
clear header_2d points_2d normal_2d;

clear header_candy xsize_candy;
clear point_locs centerline_points;

clear OV_candy XV_candy YV_candy;
clear distances_2d true_coordinates;

axes(handles.slice_axes); cla;
axes(handles.candy_cane); cla;
axes(handles.flow_axes_2d); cla;
axes(handles.ttu_2d_plot); cla;
axes(handles.ttp_2d_plot); cla;
axes(handles.ttf_2d_plot); cla;
axes(handles.xcorr_2d_plot); cla;

set(handles.slice_axes, 'Visible','off')
set(handles.flow_axes_2d, 'Visible','off')
set(handles.uipanel6, 'Visible','off')
set(handles.slice_slider_text, 'Visible','off')
set(handles.slice_slider_left, 'Visible','off')
set(handles.slice_slider_right, 'Visible','off')
set(handles.slice_slider, 'Visible','off')
set(handles.candy_cane, 'Visible','off')
set(handles.add_flow_info, 'Visible','off')
set(handles.show_ttf_plot, 'Visible','off')
set(handles.ttu_2d_plot, 'Visible','off')
set(handles.ttp_2d_plot, 'Visible','off')
set(handles.ttf_2d_plot, 'Visible','off')
set(handles.xcorr_2d_plot, 'Visible','off')
set(handles.place_point, 'Visible','off')
set(handles.load_2d_planes, 'Visible','off')
set(handles.identify_2d_planes, 'Visible','off')
set(handles.show_PWV, 'Visible','off')
set(handles.ttu_2d_pwv, 'Visible','off')
set(handles.ttp_2d_pwv, 'Visible','off')
set(handles.ttf_2d_pwv, 'Visible','off')
set(handles.xcorr_2d_pwv, 'Visible','off')

set(handles.load_2D_slice,'BackgroundColor',[0.847 0.1608 0])
set(handles.load_2D_slice,'ForegroundColor','white')
set(handles.load_candy_cane,'BackgroundColor',[0.87 0.92 0.98])
set(handles.load_candy_cane,'ForegroundColor','black')

% --- Executes on button press save_figure.
function save_figure_Callback(hObject, eventdata, handles)
% hObject    handle to save_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

selected_axes = get(handles.figure_list,'Value');
if selected_axes == 1
    newFig = figure('Parent',0,'Visible','on');
    newAxes = axes('Parent',newFig);
    clone = copyobj(get(handles.slice_axes, 'Children'), newAxes);
    axis ij
    axis off
    colormap(gray);
    axis image
elseif selected_axes == 2
    newFig = figure('Parent',0,'Visible','on');
    newAxes = axes('Parent',newFig);
    clone = copyobj(get(handles.candy_cane, 'Children'), newAxes);
    axis ij
    axis off
    colormap(gray);
    axis image
elseif selected_axes == 3
    newFig = figure('Parent',0,'Visible','on');
    newAxes = axes('Parent',newFig);
    clone = copyobj(get(handles.flow_axes_2d, 'Children'), newAxes);
    axis auto
elseif selected_axes == 4
    newFig = figure('Parent',0,'Visible','on');
    newAxes = axes('Parent',newFig);
    clone = copyobj(get(handles.ttu_2d_plot, 'Children'), newAxes);
    axis auto
elseif selected_axes == 5
    newFig = figure('Parent',0,'Visible','on');
    newAxes = axes('Parent',newFig);
    clone = copyobj(get(handles.ttp_2d_plot, 'Children'), newAxes);
    axis auto
elseif selected_axes == 6
    newFig = figure('Parent',0,'Visible','on');
    newAxes = axes('Parent',newFig);
    clone = copyobj(get(handles.ttf_2d_plot, 'Children'), newAxes);
    axis auto
elseif selected_axes == 7
    newFig = figure('Parent',0,'Visible','on');
    newAxes = axes('Parent',newFig);
    clone = copyobj(get(handles.xcorr_2d_plot, 'Children'), newAxes);
    axis auto
elseif selected_axes == 8
    axis off
    axis auto
    axis ij
    
    set(gcf, 'Renderer','OpenGL')
    opengl software

    colormap('jet');

    camlight right; 
    lighting gouraud
    alpha(0.9)
    set(gcf, 'RendererMode','Manual');
    set(gca,'Color',[0.6,0.8,1]);
    daspect([1 1 1])
    newFig = figure('Parent',0,'Visible','on');
    newAxes = axes('Parent',newFig);
    set(gca,'YDir','reverse','ZDir','reverse');
    clone = copyobj(get(handles.mip, 'Children'), newAxes);
elseif selected_axes == 9
    newFig = figure('Parent',0,'Visible','on');
    newAxes = axes('Parent',newFig);
    clone = copyobj(get(handles.mean_flow, 'Children'), newAxes);
    axis auto
elseif selected_axes == 10
    newFig = figure('Parent',0,'Visible','on');
    newAxes = axes('Parent',newFig);
    clone = copyobj(get(handles.ttu_3d_plot, 'Children'), newAxes);
    axis auto
elseif selected_axes == 11
    newFig = figure('Parent',0,'Visible','on');
    newAxes = axes('Parent',newFig);
    clone = copyobj(get(handles.ttp_3d_plot, 'Children'), newAxes);
    axis auto
elseif selected_axes == 12
    newFig = figure('Parent',0,'Visible','on');
    newAxes = axes('Parent',newFig);
    clone = copyobj(get(handles.ttf_3d_plot, 'Children'), newAxes);
    axis auto
elseif selected_axes == 13
    newFig = figure('Parent',0,'Visible','on');
    newAxes = axes('Parent',newFig);
    clone = copyobj(get(handles.xcorr_3d_plot, 'Children'), newAxes);
    axis auto
end



name = inputdlg('Name your file:','Name Your File',1);
saveas(newFig, name{1}, 'tif');

% --- Executes on drop down menu figure_list.
function figure_list_Callback(hObject, eventdata, handles)
% hObject    handle to save_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on wrap_slider.
function wrap_slider_Callback(hObject, eventdata, handles)
% hObject    handle to wrap_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global flow_t flow_array;
global plane_color;
global pwv_3d_computed;
global distance_3d distances_2d;

wrap_num = str2double(get(handles.wrap_slider,'String'));
if (wrap_num > -1) && (wrap_num < length(flow_t)-1)
    
%     set(handles.wrap_time_frames,'String',['Wrap ' num2str(wrap_num) ' Time Frames'])

    e = size(flow_array,2);

    for i = 2:size(flow_array,1)
        a = flow_array(i, e-(wrap_num-1):e);
        b = flow_array(i,1:e-wrap_num);
        flow_array(i,:) = [a b];
    end
    
    a = flow_t(e-(wrap_num-1):e);
    b = flow_t(1:e-wrap_num);
    flow_t = [a b];
    clear a b e;

    % Update mean flow figure
    axes(handles.mean_flow);
    set(gcf, 'Renderer','OpenGL')
    opengl software

    cla;

    % Handles six different colors for up to 18 planes.
    % Plots flow waveforms from all planes except for the current plane.
    % alw
    if size(flow_array,1) > 1
        for i = 2:size(flow_array,1)
            if i < 8
                plot(flow_array(1,:),flow_array(i,:),'Color',plane_color(i-1)); % removed abs
            elseif i > 7
                plot(flow_array(1,:),flow_array(i,:),'Color',plane_color(i-5));
            elseif i > 13
                plot(flow_array(1,:),flow_array(i,:),'Color',plane_color(i-11));
            end

            hold on
        end
    end

    % Handles six different colors for up to 18 planes.
    % Plots the current flow waveform from the active plane.
    % alw
    o = size(flow_array,1);
    if o < 7
        plot(flow_array(1,:),flow_t,'Color',plane_color(o)); % removed abs
    elseif o > 6
        plot(flow_array(1,:),flow_t,'Color',plane_color(o-6));
    elseif o > 12
        plot(flow_array(1,:),flow_t,'Color',plane_color(o-12));
    end
    xlabel('time [ms]');
    ylabel('velocity [ml/s]');
else
    warndlg(['There are not that many images to wrap. Pick a new number less than ' num2str(length(flow_t)-1)],'!! Warning !!');
end

if pwv_3d_computed
    num_planes = size(flow_array,1) - 1; % Subtract 1 since the first row is time information

    for i = 1:num_planes
        [c,index] = max(abs(flow_array(i+1,:))); clear c;
        index = index(1);
        
        % Begin time-to-foot (ttf) calculation
        end_point = 0;
        counter = 1;
        data_points = 0;

        for m = index-1:-1:1
            if abs(flow_array(i+1,m)) < 0.8*abs(flow_array(i+1,index))
                end_point(counter) = m;
                counter = counter+1;
            end
        end

        end_point = end_point(1);
        
        if end_point < index-1
            end_point = index-1;
        end %alw 2011-01-18
        
        begin_point = 0;

        n = end_point-1;
        while n > 0 && abs(flow_array(i+1,n)) > 0.2*abs(flow_array(i+1,index))
            begin_point = n;
            n = n-1;
        end

        if begin_point == 0
            %begin_point = 1;
            begin_point = end_point-1; % If the data point after end_point already falls below the 20% threshold, just set the beginning point of
            % the upstroke to be the next point
        else
            begin_point = begin_point(1);
        end
        clear counter;

        data_points = begin_point:1:end_point; clear begin_point; clear end_point;

        time_data = zeros(1,length(data_points)); flow_data = zeros(1,length(data_points));
        for o = 1:length(data_points)
            time_data(o) = flow_array(1,data_points(o));
            flow_data(o) = flow_array(i+1,data_points(o));
        end

        p = polyfit(time_data,abs(flow_data),1);
        ttf(i) = -1*p(2)/p(1);
        % End time-to-foot calculation
        
        figure;plot(flow_array(1,:),abs(flow_array(i+1,:)));hold on;
        string = [num2str(p(1)) '*x + ' num2str(p(2))];
        fplot(string,[0 flow_array(1,end)])
    end
    
    if get(handles.use_2d_distance, 'Value') == 1
        if not(isempty(distances_2d))
            dist = distances_2d;
        else
            warndlg('ERROR: No distances and/or centerline were defined from the 2D images!!','!! Warning !!');
        end
    else
        dist = distance_3d;
    end
    
    % Plot the Time to Foot data
    set(handles.ttf_3d_plot,'Visible','on');
    axes(handles.ttf_3d_plot);cla;
    set(gcf, 'Renderer','OpenGL') %2010-10-19 renderer change
    opengl software  %2010-10-19 renderer change
    plot(dist,ttf); ylabel('ttf [ms]');
    hold on;
    p = polyfit(dist,ttf,1); % Find the slope and intercept of a best fit line to the data
    f = @(x) p(1)*x + p(2); % Define an anonymous function, y = mx + b
    fplot(f,[0 max(dist)],':r')
    set(handles.ttf_3d_pwv,'Visible','on');
    num2str(round2(1/p(1),0.1))
    set(handles.ttf_3d_pwv,'String',[num2str(round2(1/p(1),0.1)) ' m/s']); % Display PWV
end

% --- Executes on first_last.
function first_last_Callback(hObject, eventdata, handles)
% hObject    handle to first_last (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

warndlg('If you wish to use this feature, make sure your plane information is saved prior to computing PWV.','!! Warning !!');

% --- Executes on get_flow.
function get_flow_Callback(hObject, eventdata, handles)
% hObject    handle to get_flow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% From CV Flow and using 2D slices in the ascending aorta, normal values
% are ~5500 ml/min, ~11 cm/s, and ~700 mm^2 for total flow, average
% velocity, and area in the ascending aorta, respectively. For the
% descending aorta, normal values are ~3000 ml/min, 10 cm/s, and 450 mm^2
% for total flow, average velocity, and area, respectively.

global area flow_array flow_t plane_counter;

plane_num = get(handles.plane_select,'Value');

time_interpolant = flow_array(1,1):(flow_array(1,end) - flow_array(1,1))/200:flow_array(1,end);
time_interpolant = time_interpolant/1000;
if plane_num == size(flow_array,1)
    flow_array_temp = spline(flow_array(1,:)/1000, flow_t, time_interpolant);
else
    flow_array_temp = spline(flow_array(1,:)/1000, flow_array(plane_num+1,:), time_interpolant);
end
vel_array_temp = flow_array_temp/area;

total_flow = trapz(time_interpolant,flow_array_temp)*60/(time_interpolant(end) - time_interpolant(1));
disp(['Total flow = ' num2str(total_flow) ' ml/min']);
peak_flow = max(flow_array_temp);
disp(['Peak flow = ' num2str(peak_flow) ' ml/s']);
ave_vel = mean(vel_array_temp);
disp(['Mean velocity = ' num2str(ave_vel) ' cm/s']);
vel_flow = max(vel_array_temp);
disp(['Peak velocity = ' num2str(vel_flow) ' cm/s']);

% --- Executes on unwrap.
function unwrap_Callback(hObject, eventdata, handles)
% hObject    handle to unwrap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global PD;

p = impoint;

points = p; clear p;

negative = 0;
slice_num = str2num(getString(slice_slider));

if PD(points(1),points(2),slice_num) < 1
    negative = 1;
end

done = 0;
point_count = 1;

while not(done)
    points_to_analyze(1) = PD(points(point_count,1)-1,points(point_count,2)+1,slice_num);
    points_to_analyze(2) = PD(points(point_count,1),points(point_count,2)+1,slice_num);
    points_to_analyze(3) = PD(points(point_count,1)+1,points(point_count,2)+1,slice_num);
    points_to_analyze(4) = PD(points(point_count,1)-1,points(point_count,2),slice_num);
    points_to_analyze(5) = PD(points(point_count,1)+1,points(point_count,2),slice_num);
    points_to_analyze(6) = PD(points(point_count,1)-1,points(point_count,2)-1,slice_num);
    points_to_analyze(7) = PD(points(point_count,1),points(point_count,2)-1,slice_num);
    points_to_analyze(8) = PD(points(point_count,1)+1,points(point_count,2)-1,slice_num);
    
    for i = 1:8
        if negative
            if points_to_analyze(i) < 0
                point_count = point_count + 1;
                if i == 1
                    points(point_count,:) = [(points(point_count-1,1)-1) (points(point_count-1,2)+1)];
                    if point_count > 1
                        for j = point_count-1:-1:1
                            if points(point_count,:) == points(k,:)
                                points(point_count,:) = [];
                                point_count = point_count - 1;
                            end
                        end
                    end                    
                elseif i == 2
                    points(point_count,:) = [points(point_count-1,1) (points(point_count-1,2)+1)];
                    if point_count > 1
                        for j = point_count-1:-1:1
                            if points(point_count,:) == points(k,:)
                                points(point_count,:) = [];
                                point_count = point_count - 1;
                            end
                        end
                    end   
                elseif i == 3
                    points(point_count,:) = [(points(point_count-1,1)+1) (points(point_count-1,2)+1)];
                    if point_count > 1
                        for j = point_count-1:-1:1
                            if points(point_count,:) == points(k,:)
                                points(point_count,:) = [];
                                point_count = point_count - 1;
                            end
                        end
                    end   
                elseif i == 4
                    points(point_count,:) = [(points(point_count-1,1)-1) points(point_count-1,2)];
                    if point_count > 1
                        for j = point_count-1:-1:1
                            if points(point_count,:) == points(k,:)
                                points(point_count,:) = [];
                                point_count = point_count - 1;
                            end
                        end
                    end   
                elseif i == 5
                    points(point_count,:) = [(points(point_count-1,1)+1) points(point_count-1,2)];
                    if point_count > 1
                        for j = point_count-1:-1:1
                            if points(point_count,:) == points(k,:)
                                points(point_count,:) = [];
                                point_count = point_count - 1;
                            end
                        end
                    end   
                elseif i == 6
                    points(point_count,:) = [(points(point_count-1,1)-1) (points(point_count-1,2)-1)];
                    if point_count > 1
                        for j = point_count-1:-1:1
                            if points(point_count,:) == points(k,:)
                                points(point_count,:) = [];
                                point_count = point_count - 1;
                            end
                        end
                    end   
                elseif i == 7
                    points(point_count,:) = [points(point_count-1,1) (points(point_count-1,2)-1)];
                    if point_count > 1
                        for j = point_count-1:-1:1
                            if points(point_count,:) == points(k,:)
                                points(point_count,:) = [];
                                point_count = point_count - 1;
                            end
                        end
                    end   
                elseif i == 8
                    points(point_count,:) = [(points(point_count-1,1)+1) (points(point_count-1,2)-1)];
                    if point_count > 1
                        for j = point_count-1:-1:1
                            if points(point_count,:) == points(k,:)
                                points(point_count,:) = [];
                                point_count = point_count - 1;
                            end
                        end
                    end
                end
            end
        else
            if points_to_analyze(i) > 0
            end
        end
    end 
end
    
% --- Executes on getFlow2d.
function getFlow2d_Callback(hObject, eventdata, handles)
% hObject    handle to getFlow2d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global time_current flow_current area_2d;

time = time_current/1000;

area = area_2d/100;

total_flow = trapz(time,flow_current)*60;
disp(['Total flow = ' num2str(total_flow) ' ml/min']);
peak_flow = max(flow_current);
disp(['Peak flow = ' num2str(peak_flow) ' ml/s']);
ave_vel = mean(flow_current/area);
disp(['Mean velocity = ' num2str(ave_vel) ' cm/s']);
vel_flow = max(flow_current/area);
disp(['Peak velocity = ' num2str(vel_flow) ' cm/s']);


function value = findVal(fields,values,field,length)
    index = find(cellfun(@(s) strcmp(field, s), fields));
    value = cellfun(@str2num,values(index:(index+length-1)));
