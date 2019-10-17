function varargout = pressure_gui(varargin)
% PRESSURE_GUI M-file for pressure_gui.fig
%      PRESSURE_GUI, by itself, creates a new PRESSURE_GUI or raises the existing
%      singleton*.
%
%      H = PRESSURE_GUI returns the handle to a new PRESSURE_GUI or the handle to
%      the existing singleton*.
%
%      PRESSURE_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PRESSURE_GUI.M with the given input arguments.
%
%      PRESSURE_GUI('Property','Value',...) creates a new PRESSURE_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pressure_gui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pressure_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help pressure_gui

% Last Modified by GUIDE v2.5 23-Jun-2007 12:56:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pressure_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @pressure_gui_OutputFcn, ...
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


% --- Executes just before pressure_gui is made visible.
function pressure_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pressure_gui (see VARARGIN)

% Choose default command line output for pressure_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pressure_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = pressure_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function viscosity_Callback(hObject, eventdata, handles)
% hObject    handle to viscosity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of viscosity as text
%        str2double(get(hObject,'String')) returns contents of viscosity as a double


% --- Executes during object creation, after setting all properties.
function viscosity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to viscosity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function density_Callback(hObject, eventdata, handles)
% hObject    handle to density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of density as text
%        str2double(get(hObject,'String')) returns contents of density as a double


% --- Executes during object creation, after setting all properties.
function density_CreateFcn(hObject, eventdata, handles)
% hObject    handle to density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function max_iter_Callback(hObject, eventdata, handles)
% hObject    handle to max_iter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max_iter as text
%        str2double(get(hObject,'String')) returns contents of max_iter as a double


% --- Executes during object creation, after setting all properties.
function max_iter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max_iter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pressure_calc_type.
function pressure_calc_type_Callback(hObject, eventdata, handles)
% hObject    handle to pressure_calc_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns pressure_calc_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pressure_calc_type


% --- Executes during object creation, after setting all properties.
function pressure_calc_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pressure_calc_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in grad_update.
function grad_update_Callback(hObject, eventdata, handles)
% hObject    handle to grad_update (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

calc_pressure(handles);


function calc_pressure(handles)

global tframes;
%Raw Data
global MAG;
global CD;
global MASK;
global VELX;
global VELY;
global VELZ;
global VELXt;
global VELYt;
global VELZt;
global sMASK;
global sCD;
global sMAG;
global GRADx;
global GRADy;
global GRADz;
global plist;

global delX;
global delY;
global delZ;
global tres;

visc = str2double(get(handles.viscosity,'String'))/1000;
dens = str2double(get(handles.density,'String'));
max_iter = str2double(get(handles.max_iter,'String'));
alpha = str2double(get(handles.alpha,'String'));
poly_num =  str2double(get(handles.poly_num,'String'));

calc_type = get(handles.pressure_calc_type,'Value');
%%% 1 = Time Dependent
%%% 2 = Time Independent

%%%%%%%%%%%%%%%%%%%%%STEP 1: Setup For Calc%%%%%%%%%%%%%

%%Points to use
plist = find(sMASK == 1);

if( size(GRADx) ~= size(plist))
    GRADx = zeros(size(plist));
    GRADy = zeros(size(plist));
    GRADz = zeros(size(plist));
end

%%%%%%%%%%%%%%%%%%%%%STEP 2: Gradient Calc%%%%%%%%%%%%%

DIM = size(sMASK);
npts = length(plist)


%%%DERIVATIVE CONVERSIONS
conv_vel = 1/1000; %mm/s to m/s
conv_vel2= 1/1000/1000; 

conv_dx = 1/delX*1000; % mm to m
conv_dx2= 1/delX/delX*1000*1000;

conv_dy = 1/delY*1000;
conv_dy2= 1/delY/delY*1000*1000;

conv_dz = 1/delZ*1000;
conv_dz2= 1/delZ/delZ*1000*1000;

conv_dt = 1/tres*1000;


for pos = 1:50 %npts
    
    disp(['Point Number ',int2str(pos),' of ',int2str(npts)])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%   Step 2.1 Setup Derivative Space                %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%Setup Derivative Space For Derivatives
    [x0 y0 z0] = ind2sub(DIM,plist(pos));

    %%%X Points to be used
    plistX(1,:)= [x0 y0 z0];
    points = 1;
    xP = x0;
    stop_flag = 0;

    while stop_flag == 0
        xP = xP - 1;
        if xP < 1
            stop_flag = 1;
        else
            if sMASK(xP,y0,z0) == 1
                points = points + 1;
                plistX(points,:) = [xP y0 z0];
                if points > poly_num
                    stop_flag = 1;
                end
            end
        end
    end
    
    xP = x0;
    stop_flag =0;
    while stop_flag == 0
        xP = xP + 1;
        if xP > DIM(1)
            stop_flag = 1;
        else
            if sMASK(xP,y0,z0) == 1
                points = points + 1;
                plistX(points,:) = [xP y0 z0];
                if points > 2*poly_num -1
                    stop_flag = 1;
                end
            end
        end
    end
    
    %%%Y Points to be used
    plistY(1,:)= [x0 y0 z0];
    points = 1;
    yP = y0;
    stop_flag = 0;

    while stop_flag == 0
        yP = yP - 1;
        if yP < 1
            stop_flag = 1;
        else
            if sMASK(x0,yP,z0) == 1
                points = points + 1;
                plistY(points,:) = [x0 yP z0];
                if points > poly_num
                    stop_flag = 1;
                end
            end
        end
    end
    
    yP = y0;
    stop_flag =0;
    while stop_flag == 0
        yP = yP + 1;
        if yP > DIM(2)
            stop_flag = 1;
        else
            if sMASK(x0,yP,z0) == 1
                points = points + 1;
                plistY(points,:) = [x0 yP z0];
                if points > 2*poly_num -1
                    stop_flag = 1;
                end
            end
        end
    end
    
    %%%Z Points to be used
    plistZ(1,:)= [x0 y0 z0];
    points = 1;
    zP = z0;
    stop_flag = 0;

    while stop_flag == 0
        zP = zP - 1;
        if zP < 1
            stop_flag = 1;
        else
            if sMASK(x0,y0,zP) == 1
                points = points + 1;
                plistZ(points,:) = [x0 y0 zP];
                if points > poly_num
                    stop_flag = 1;
                end
            end
        end
    end
    
    zP = z0;
    stop_flag =0;
    while stop_flag == 0
        zP = zP + 1;
        if zP > DIM(3)
            stop_flag = 1;
        else
            if sMASK(x0,y0,zP) == 1
                points = points + 1;
                plistZ(points,:) = [x0 y0 zP];
                if points > 2*poly_num - 1
                    stop_flag = 1;
                end
            end
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%   Step 2.2 Get fit values                        %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    INDX = sub2ind(DIM,plistX(:,1),plistX(:,2),plistX(:,3));
    INDY = sub2ind(DIM,plistY(:,1),plistY(:,2),plistY(:,3));
    INDZ = sub2ind(DIM,plistZ(:,1),plistZ(:,2),plistZ(:,3));
    
    Xvals = plistX(:,1);
    Yvals = plistY(:,2);
    Zvals = plistZ(:,3);
    
    poly_orderx = length(Xvals) - 1;
    poly_ordery = length(Yvals) - 1;
    poly_orderz = length(Zvals) - 1;
    
    %%%%%%GET MATRIX OF ENCODING%%%%%%%
    EX = repmat(Xvals,1,poly_orderx+1).^repmat( (poly_orderx:-1:0),length(Xvals),1)
    EY = repmat(Yvals,1,poly_ordery+1).^repmat( (poly_ordery:-1:0),length(Yvals),1)
    EZ = repmat(Zvals,1,poly_orderz+1).^repmat( (poly_orderz:-1:0),length(Zvals),1)
    
    EXT = EX';
    EYT = EY';
    EZT = EZ';
    
    EXF = EXT*EX;
    EYF = EYT*EY;
    EZF = EZT*EZ;
    

    for time = 1:tframes
        
        
        time_start = time -1 ;
        time_end   = time +1 ;
        
        if time_start < 1
            time_start = 1;
            time_end   = 3;
        elseif time_end > tframes
            time_end = tframes;
            time_start= tframes -2;
        end
            
        %%%%%Fit values in Time%%%%
        vxt = VELXt(x0,y0,z0,time_start:time_end);
        vyt = VELYt(x0,y0,z0,time_start:time_end);
        vzt = VELZt(x0,y0,z0,time_start:time_end);
         
        fit_vx_dt = polyfit(time_start:time_end',vxt(:)',2);
        fit_vy_dt = polyfit(time_start:time_end',vyt(:)',2);
        fit_vz_dt = polyfit(time_start:time_end',vzt(:)',2);
    
        INDXt = INDX + (time-1)*DIM(1)*DIM(2)*DIM(3);
        INDYt = INDY + (time-1)*DIM(1)*DIM(2)*DIM(3);
        INDZt = INDZ + (time-1)*DIM(1)*DIM(2)*DIM(3);
         
         poly_orderx
         poly_ordery
         poly_orderz
        
        
        if(poly_orderx > 0)
        %%%%%Fit values in X%%%%
        vxV = VELXt(INDXt);
        vxF = EXT*vxV;
        fit_vx_dx = polyfit(Xvals,vxV);
      
        vyV = VELYt(INDXt);
        vyF = EXT*vyV;
        fit_vy_dx = linsolve(EX,vyV);
        
        vzV = VELZt(INDXt);
        vzF = EXT*vzV;
        fit_vz_dx = linsolve(EX,vzV);
        else
           fit_vx_dx = 0;
           fit_vy_dx = 0;
           fit_vz_dx = 0;
        end
        
        if(poly_ordery > 0)
        %%%%%Fit values in Y%%%%
        vxV = VELXt(INDYt);
        vxF = EYT*vxV;
        fit_vx_dy = linsolve(EY,vxV);
        
        vyV = VELYt(INDYt);
        vyF = EYT*vyV;
        fit_vy_dy = linsolve(EY,vyV);
        
        vzV = VELZt(INDYt);
        vzF = EYT*vzV;
        fit_vz_dy = linsolve(EY,vzV);
        else
           fit_vx_dy = 0;
           fit_vy_dy = 0;
           fit_vz_dy = 0;
        end
        
        if(poly_orderz > 0)
        %%%%%Fit values in X%%%%
        vxV = VELXt(INDZt);
        vxF = EZT*vxV;
        fit_vx_dz = linsolve(EZ,vxV);
      
        vyV = VELYt(INDZt);
        vyF = EZT*vyV;
        fit_vy_dz = linsolve(EZ,vyV);
        
        vzV = VELZt(INDZt);
        vzF = EZT*vzV;
        fit_vz_dz = linsolve(EZ,vzV);
        else
           fit_vx_dz = 0;
           fit_vy_dz = 0;
           fit_vz_dz = 0;
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%   Step 2.3 Derivatives                           %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%velocity Terms
        vx = VELXt(x0,y0,z0,time);
        vy = VELYt(x0,y0,z0,time);
        vz = VELZt(x0,y0,z0,time);
        
        %%First Derivatives
        dvxdx = polyval( polyder(fit_vx_dx),x0);
        dvxdy = polyval( polyder(fit_vx_dy),y0);
        dvxdz = polyval( polyder(fit_vx_dx),z0);
        
        dvydx = polyval( polyder(fit_vy_dx),x0);
        dvydy = polyval( polyder(fit_vy_dy),y0);
        dvydz = polyval( polyder(fit_vy_dx),z0);
        
        dvzdx = polyval( polyder(fit_vz_dx),x0);
        dvzdy = polyval( polyder(fit_vz_dy),y0);
        dvzdz = polyval( polyder(fit_vz_dx),z0);
        
        dvxdt = polyval( polyder(fit_vx_dt),x0);
        dvydt = polyval( polyder(fit_vy_dt),y0);
        dvzdt = polyval( polyder(fit_vz_dt),z0);
        
        %%2nd Derivatives
        dvxdx2 = polyval( polyder(polyder(fit_vx_dx)),x0);
        dvxdy2 = polyval( polyder(polyder(fit_vx_dy)),y0);
        dvxdz2 = polyval( polyder(polyder(fit_vx_dx)),z0);
        
        dvydx2 = polyval( polyder(polyder(fit_vy_dx)),x0);
        dvydy2 = polyval( polyder(polyder(fit_vy_dy)),y0);
        dvydz2 = polyval( polyder(polyder(fit_vy_dx)),z0);
        
        dvzdx2 = polyval( polyder(polyder(fit_vz_dx)),x0);
        dvzdy2 = polyval( polyder(polyder(fit_vz_dy)),y0);
        dvzdz2 = polyval( polyder(polyder(fit_vz_dx)),z0);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%   Step 2.4 Navier-Stokes                         %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        GRADx(pos,time) =-conv_vel*( -conv_dt*dens*( dvxdt ) - conv_vel*( conv_dx*dens*vx*dvxdx - conv_dy*dens*vy*dvxdy - conv_dz*dens*vz*dvxdz )+ ...
            conv_dx2*visc*dvxdx2 + conv_dy2*visc*dvxdy2 + conv_dz2*visc*dvxdz2 );  
        GRADy(pos,time) =-conv_vel*( -conv_dt*dens*( dvydt ) - conv_vel*( conv_dx*dens*vx*dvydx - conv_dy*dens*vy*dvydy - conv_dz*dens*vz*dvydz )+ ...
            conv_dx2*visc*dvydx2 + conv_dy2*visc*dvydy2 + conv_dz2*visc*dvydz2 );
        GRADz(pos,time) =-conv_vel*( -conv_dt*dens*( dvzdt ) - conv_vel*( conv_dx*dens*vx*dvzdx - conv_dy*dens*vy*dvzdy - conv_dz*dens*vz*dvzdz )+ ...
            conv_dx2*visc*dvzdx2 + conv_dy2*visc*dvzdy2 + conv_dz2*visc*dvzdz2 ); 
        
    end
    
    clear plistY
    clear plistX
    clear plistZ
end


function update_image(handles)

global sMASK;
global GRADx;
global GRADy;
global GRADz;
global PRESSURE;
global plist;
global press_axis;
global tframes;
global IMAGE;

global VELX;
global VELY;
global VELZ;
global VELXt;
global VELYt;
global VELZt;

image_object = get(handles.image_object,'Value');
%%% 1 = Gradient X
%%% 2 = Gradient Y
%%% 3 = Gradient Z
%%% 4 = Velocity X 
%%% 5 = Velocity Y 
%%% 6 = Velocity Z
%%% 7 = Pressure 

image_type = get(handles.image_type,'Value');
%%% 1 = Time Resolved
%%% 2 = Time Averaged

image_dir = get(handles.image_dir,'Value');
%%% 1 = X
%%% 2 = Y
%%% 3 = Z

num_frames = tframes;
if image_type == 2 
    num_frames = 1;
end

new_fig = ishandle(press_axis)
figure(press_axis)

daspect([1 1 1])

if(new_fig==0)
scrsz = get(0,'ScreenSize');
SIZE=[(scrsz(3)*1/2-64) scrsz(4)*1/4 scrsz(3)*1/2 scrsz(4)*1/2];
set(gcf,'Position',SIZE);
set(gcf,'Name','Pressure Window');
end


IMAGE = zeros(size(sMASK));
DIM=size(sMASK);
VENC = 600;

for t=1:num_frames

if image_object == 1
    IMAGE(plist)=GRADx(:,t);
    min_im = min(GRADx(:));
    max_im = max(GRADx(:));
elseif image_object == 2
    IMAGE(plist)=GRADy(:,t);
    min_im = min(GRADy(:));
    max_im = max(GRADy(:));
elseif image_object == 3
    IMAGE(plist)=GRADz(:,t);
    min_im = min(GRADy(:));
    max_im = max(GRADz(:));
elseif image_object == 4
    time_plist = plist + (t-1)*DIM(1)*DIM(2)*DIM(3);
    IMAGE(plist)=VELXt(time_plist);
    min_im = -VENC;
    max_im = VENC;
elseif image_object == 5
    time_plist = plist + (t-1)*DIM(1)*DIM(2)*DIM(3);
    IMAGE(plist)=VELYt(time_plist);   
    min_im = -VENC;
    max_im = VENC;
elseif image_object == 6
    time_plist = plist + (t-1)*DIM(1)*DIM(2)*DIM(3);
    IMAGE(plist)=VELZt(time_plist);     
    min_im = -VENC;
    max_im = VENC;
elseif image_object == 7
    IMAGE(plist)=0.007501*PRESSURE(:,t); 
    min_im = 0.007501*min(PRESSURE(:));
    max_im = 0.007501*max(PRESSURE(:));
end

%%%%%%%NOW DISPLAY%%%%%
if( image_dir == 1 )
    d1 = 2;
    d2 = 3;
elseif image_dir == 2
    d1 = 1;
    d2 = 3;
elseif image_dir == 3
    d1 = 1;
    d2 = 2;
end
figure(press_axis)
imagesc(reshape(sum(IMAGE,image_dir),[size(IMAGE,d1) size(IMAGE,d2)])./reshape(sum(sMASK,image_dir),[size(IMAGE,d1) size(IMAGE,d2)]),[min_im max_im])
colorbar('FontSize',22)
set(gca,'XTick',[],'YTick',[]);
daspect([1 1 1]);
pause(0.1);

end



function alpha_Callback(hObject, eventdata, handles)
% hObject    handle to alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha as text
%        str2double(get(hObject,'String')) returns contents of alpha as a double


% --- Executes during object creation, after setting all properties.
function alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
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




% --- Executes on selection change in image_object.
function image_object_Callback(hObject, eventdata, handles)
% hObject    handle to image_object (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns image_object contents as cell array
%        contents{get(hObject,'Value')} returns selected item from image_object


% --- Executes during object creation, after setting all properties.
function image_object_CreateFcn(hObject, eventdata, handles)
% hObject    handle to image_object (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in image_type.
function image_type_Callback(hObject, eventdata, handles)
% hObject    handle to image_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns image_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from image_type


% --- Executes during object creation, after setting all properties.
function image_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to image_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in image_update.
function image_update_Callback(hObject, eventdata, handles)
% hObject    handle to image_update (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update_image(handles);


% --- Executes on selection change in image_dir.
function image_dir_Callback(hObject, eventdata, handles)
% hObject    handle to image_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns image_dir contents as cell array
%        contents{get(hObject,'Value')} returns selected item from image_dir


% --- Executes during object creation, after setting all properties.
function image_dir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to image_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in pressure_update.
function pressure_update_Callback(hObject, eventdata, handles)
% hObject    handle to pressure_update (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

update_pressure(handles);



function update_pressure(handles)

global sMASK;
global GRADx;
global GRADy;
global GRADz;
global PRESSURE;
global plist;
global xlist;
global ylist;
global zlist;

global press_axis;
global tframes;
global IMAGE;

global delX;
global delY;
global delZ;

global VELX;
global VELY;
global VELZ;
global VELXt;
global VELYt;
global VELZt;

max_iter = str2double(get(handles.max_iter,'String'));
alpha = str2double(get(handles.alpha,'String'));

%%%%%%Step 1. Get Neighbors of Each Pixel %%%%%%%
points = size(plist,1);

[xlist ylist zlist] = ind2sub(size(sMASK),plist);
nbs=zeros(points,6);

for pos = 1:points
    x0 = xlist(pos);
    y0 = ylist(pos);
    z0 = zlist(pos);
    
    
    xeq = xlist == x0;
    yeq = ylist == y0;
    zeq = zlist == z0;
    
    nbs_xp = find( (xlist ==(x0+1) ).*yeq.*zeq,1 );
    if size(nbs_xp,1)~= 0
        nbs(pos,1) = nbs_xp;
    end
        
    nbs_xm = find( (xlist ==(x0-1) ).*yeq.*zeq,1 );
    if size(nbs_xm,1)~= 0
        nbs(pos,2) = nbs_xm;
    end
    
    
    nbs_yp = find( (ylist ==(y0+1) ).*xeq.*zeq,1 );
    if size(nbs_yp,1)~= 0
        nbs(pos,3) = nbs_yp;
    end
    
    nbs_ym = find( (ylist ==(y0-1) ).*xeq.*zeq,1 );
    if size(nbs_ym,1)~= 0
        nbs(pos,4) = nbs_ym;
    end
    
    
    nbs_zp = find( ( zlist==(z0+1) ).*yeq.*xeq,1 );
    if size(nbs_zp,1)~= 0
        nbs(pos,5) = nbs_zp;
    end
    
    nbs_zm = find( ( zlist==(z0-1) ).*yeq.*xeq,1 );
    if size(nbs_zm,1)~= 0
        nbs(pos,6) = nbs_zm;
    end
end


%%%%%%Step 2. Initial Setup (Simple Integration) %%%%%%%

start_pos =  floor( 1 + rand()*points );
num_counted = 1;
PRESSURE(start_pos,:)= 0*GRADx(1,:);
cnt_idx(1)=start_pos;
counted = zeros(points,1);
counted(start_pos)=1;
start_cnt = 1;

PRESSURE(:) = 0;
while num_counted < points

    num_countedp = num_counted;
    for pos = start_cnt:num_countedp

        idx = cnt_idx(pos);

        xp = nbs(idx,1);
        xm = nbs(idx,2);
        yp = nbs(idx,3);
        ym = nbs(idx,4);
        zp = nbs(idx,5);
        zm = nbs(idx,6);

        if xp ~= 0
            if counted(xp)==0
                counted(xp) = 1;
                num_counted = num_counted+1;
                cnt_idx(num_counted)=xp;
                PRESSURE(xp,:) = PRESSURE(idx,:) + GRADx(idx,:)*delX/1000;
            end
        end

        if xm ~= 0
            if counted(xm)==0
                counted(xm) = 1;
                num_counted = num_counted+1;
                cnt_idx(num_counted)=xm;
                PRESSURE(xm,:) = PRESSURE(idx,:) - GRADx(idx,:)*delX/1000;
            end
        end

        if yp ~= 0
            if counted(yp)==0
                counted(yp) = 1;
                num_counted = num_counted+1;
                cnt_idx(num_counted)=yp;
                PRESSURE(yp,:) = PRESSURE(idx,:) + GRADy(idx,:)*delY/1000;
            end
        end

        if ym ~= 0
            if counted(ym)==0
                counted(ym) = 1;
                num_counted = num_counted+1;
                cnt_idx(num_counted)=ym;
                PRESSURE(ym,:) = PRESSURE(idx,:) - GRADy(idx,:)*delY/1000;
            end
        end

        if zp ~= 0
            if counted(zp)==0
                counted(zp) = 1;
                num_counted = num_counted+1;
                cnt_idx(num_counted)=zp;
                PRESSURE(zp,:) = PRESSURE(idx,:) + GRADz(idx,:)*delZ/1000;
            end
        end

        if zm ~= 0
            if counted(zm)==0
                counted(zm) = 1;
                num_counted = num_counted+1;
                cnt_idx(num_counted)=zm;
                PRESSURE(zm) = PRESSURE(idx) - GRADz(idx)*delZ/1000;
            end
        end
        start_cnt = pos;

    end
end

        
%%%%%%%%%%%Step 3 ITERATE TO CLEAN UP
PRESSURE_old = PRESSURE;
errord = 999;
while errord > 1
    
    PRESSURE_old = PRESSURE;
    for pos = 1:points
        idx = pos;
        xp = nbs(pos,1);
        xm = nbs(pos,2);
        yp = nbs(pos,3);
        ym = nbs(pos,4);
        zp = nbs(pos,5);
        zm = nbs(pos,6);
        
        GRAD_TERM = 0;
        nb_count = 0;
        
        if xp ~= 0
             GRAD_TERM = GRAD_TERM + PRESSURE_old(xp,:) - GRADx(idx,:)*delX/1000;
             nb_count = nb_count + 1;
        end

        if xm ~= 0
             GRAD_TERM =  GRAD_TERM + PRESSURE_old(xm,:) + GRADx(idx,:)*delX/1000;
             nb_count = nb_count + 1;
        end

        if yp ~= 0
             GRAD_TERM =  GRAD_TERM + PRESSURE_old(yp,:) - GRADy(idx,:)*delY/1000;
             nb_count = nb_count + 1;
        end

        if ym ~= 0
             GRAD_TERM =  GRAD_TERM + PRESSURE_old(ym,:) + GRADy(idx,:)*delY/1000;
             nb_count = nb_count + 1;
        end
        
        if zp ~= 0
             GRAD_TERM =  GRAD_TERM + PRESSURE_old(zp,:) - GRADz(idx,:)*delZ/1000;
             nb_count = nb_count + 1;
        end

        if zm ~= 0
             GRAD_TERM =  GRAD_TERM + PRESSURE_old(zm,:) + GRADz(idx,:)*delZ/1000;
             nb_count = nb_count + 1;
        end
                
        PRESSURE(idx,:) = (1- alpha)*PRESSURE_old(idx,:) + alpha*GRAD_TERM/nb_count;
    end

    errord = sum( abs(PRESSURE(:)- PRESSURE_old(:)))/numel(PRESSURE)
    pause(0.1);
    
end

























