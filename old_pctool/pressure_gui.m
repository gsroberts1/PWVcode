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

% Last Modified by GUIDE v2.5 25-Mar-2008 13:20:02

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

global POINTS_FOUND;
POINTS_FOUND = 0;

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

global m_xlength;
global m_ylength;
global m_zlength;


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

if( length(size(GRADx)) == 2 )
    GRADx = zeros(size(plist,1),size(plist,2),tframes);
    GRADy = zeros(size(plist,1),size(plist,2),tframes);
    GRADz = zeros(size(plist,1),size(plist,2),tframes);
end

%%%%%%%%%%%%%%%%%%%%%STEP 2: Gradient Calc%%%%%%%%%%%%%

DIM = size(sMASK);
points = length(plist);


%%%DERIVATIVE CONVERSIONS
conv_vel = -1/1000; %mm/s to m/s
conv_vel2= 1/1000/1000;

conv_dx = 1/delX*1000; % mm to m
conv_dx2= 1/delX/delX*1000*1000;

conv_dy = 1/delY*1000;
conv_dy2= 1/delY/delY*1000*1000;

conv_dz = 1/delZ*1000;
conv_dz2= 1/delZ/delZ*1000*1000;

conv_dt = 1/tres*1000;

for pos = 1:points


    if 50*floor(pos/50)==pos
        disp(['Point Number ',int2str(pos),' of ',num2str(points)])
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%   Step 2.1 Setup Derivative Space                %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    %%%%Setup Derivative Space For Derivatives
    [x0 y0 z0] = ind2sub(DIM,plist(pos));

    TEST = 2;
    if TEST == 1
        GRADx(pos,:) = VELXt(x0,y0,z0,:).^3;
        GRADy(pos,:) = VELYt(x0,y0,z0,:).^3;
        GRADz(pos,:) = VELZt(x0,y0,z0,:).^3;
    elseif TEST == 2

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%   Step 2.3 Derivatives                           %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if ( (x0 < 2 ) || (y0 < 2 ) || (z0 < 2) || (x0 > m_xlength-1) || (y0 > m_ylength-1) || (z0 > m_zlength-1))
            GRADx(pos,:) = 0;
            GRADy(pos,:) = 0;
            GRADz(pos,:) = 0;
        else

        for time = 1:tframes

            time_plus1 = mod( time,tframes)+1;
            time_minus1= mod( time-2,tframes)+1;

            %%%velocity Terms
            vx = conv_vel*VELXt(x0,y0,z0,time);
            vy = conv_vel*VELYt(x0,y0,z0,time);
            vz = conv_vel*VELZt(x0,y0,z0,time);

            %%%mask values
            mxp = sMASK(x0+1,y0,z0);
            myp = sMASK(x0,y0+1,z0);
            mzp = sMASK(x0,y0,z0+1);
            mxm = sMASK(x0-1,y0,z0);
            mym = sMASK(x0,y0-1,z0);
            mzm = sMASK(x0,y0,z0-1);
            
            %%First Derivatives
            dvxdt = conv_vel*( VELXt(x0,y0,z0,time_plus1) - VELXt(x0,y0,z0,time_minus1) )/(2*tres/1000);
            dvydt = conv_vel*( VELYt(x0,y0,z0,time_plus1) - VELYt(x0,y0,z0,time_minus1) )/(2*tres/1000);
            dvzdt = conv_vel*( VELZt(x0,y0,z0,time_plus1) - VELZt(x0,y0,z0,time_minus1) )/(2*tres/1000);
            %dvxdt = conv_vel*( VELXt(x0,y0,z0,time_plus1) - VELXt(x0,y0,z0,time) )/(tres/1000);
            %dvydt = conv_vel*( VELYt(x0,y0,z0,time_plus1) - VELYt(x0,y0,z0,time) )/(tres/1000);
            %dvzdt = conv_vel*( VELZt(x0,y0,z0,time_plus1) - VELZt(x0,y0,z0,time) )/(tres/1000);

            
            dvxdx = conv_vel*( VELXt(x0+1,y0,z0,time)*mxp - VELXt(x0-1,y0,z0,time)*mxm )/(2*delX/1000);
            dvxdy = conv_vel*( VELXt(x0,y0+1,z0,time)*myp - VELXt(x0,y0-1,z0,time)*mym )/(2*delY/1000);
            dvxdz = conv_vel*( VELXt(x0,y0,z0+1,time)*mzp - VELXt(x0,y0,z0-1,time)*mzm )/(2*delZ/1000);

            dvydx = conv_vel*( VELYt(x0+1,y0,z0,time)*mxp - VELYt(x0-1,y0,z0,time)*mxm )/(2*delX/1000);
            dvydy = conv_vel*( VELYt(x0,y0+1,z0,time)*myp - VELYt(x0,y0-1,z0,time)*mym )/(2*delY/1000);
            dvydz = conv_vel*( VELYt(x0,y0,z0+1,time)*mzp - VELYt(x0,y0,z0-1,time)*mzm )/(2*delZ/1000);

            dvzdx = conv_vel*( VELZt(x0+1,y0,z0,time)*mxp - VELZt(x0-1,y0,z0,time)*mxm )/(2*delX/1000);
            dvzdy = conv_vel*( VELZt(x0,y0+1,z0,time)*myp - VELZt(x0,y0-1,z0,time)*mym )/(2*delY/1000);
            dvzdz = conv_vel*( VELZt(x0,y0,z0+1,time)*mzp - VELZt(x0,y0,z0-1,time)*mzm )/(2*delZ/1000);

            %%2nd Derivatives
            dvxdx2 = conv_vel*( VELXt(x0+1,y0,z0,time)*mxp -2*VELXt(x0,y0,z0,time) + VELXt(x0-1,y0,z0,time)*mxm )/(delX/1000)^2;
            dvxdy2 = conv_vel*( VELXt(x0,y0+1,z0,time)*myp -2*VELXt(x0,y0,z0,time) + VELXt(x0,y0-1,z0,time)*mym )/(delY/1000)^2;
            dvxdz2 = conv_vel*( VELXt(x0,y0,z0+1,time)*mzp -2*VELXt(x0,y0,z0,time) + VELXt(x0,y0,z0-1,time)*mzm )/(delZ/1000)^2;

            dvydx2 = conv_vel*( VELYt(x0+1,y0,z0,time)*mxp -2*VELYt(x0,y0,z0,time) + VELYt(x0-1,y0,z0,time)*mxm )/(delX/1000)^2;
            dvydy2 = conv_vel*( VELYt(x0,y0+1,z0,time)*myp -2*VELYt(x0,y0,z0,time) + VELYt(x0,y0-1,z0,time)*mym )/(delY/1000)^2;
            dvydz2 = conv_vel*( VELYt(x0,y0,z0+1,time)*mzp -2*VELYt(x0,y0,z0,time) + VELYt(x0,y0,z0-1,time)*mzm )/(delZ/1000)^2;

            dvzdx2 = conv_vel*( VELZt(x0+1,y0,z0,time)*mxp -2*VELZt(x0,y0,z0,time) + VELZt(x0-1,y0,z0,time)*mxm )/(delX/1000)^2;
            dvzdy2 = conv_vel*( VELZt(x0,y0+1,z0,time)*myp -2*VELZt(x0,y0,z0,time) + VELZt(x0,y0-1,z0,time)*mym )/(delY/1000)^2;
            dvzdz2 = conv_vel*( VELZt(x0,y0,z0+1,time)*mzp -2*VELZt(x0,y0,z0,time) + VELZt(x0,y0,z0-1,time)*mzm )/(delZ/1000)^2;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%   Step 2.4 Navier-Stokes                         %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            GRADx(pos,time) =   -dens*dvxdt ...
                -dens*vx*dvxdx  ...
                -dens*vy*dvxdy  ...
                -dens*vz*dvxdz ...
                +visc*dvxdx2 ...
                +visc*dvxdy2 ...
                +visc*dvxdz2;
            GRADy(pos,time) =  -dens*dvydt  ...
                -dens*vx*dvydx ...
                -dens*vy*dvydy ...
                -dens*vz*dvydz ...
                +visc*dvydx2 ...
                +visc*dvydy2 ...
                +visc*dvydz2;
            GRADz(pos,time) =  -dens*dvzdt    ...
                -dens*vx*dvzdx ...
                -dens*vy*dvzdy ...
                -dens*vz*dvzdz  ...
                +visc*dvzdx2 ...
                +visc*dvzdy2 ...
                +visc*dvzdz2;
        end
        end
    else
        %%%X Points to be used
        plistX(1,:)= [x0 y0 z0];
        lpoints = 1;
        xP = x0;
        stop_flag = 0;

        while stop_flag == 0
            xP = xP - 1;
            if xP < 1
                stop_flag = 1;
            elseif sMASK(xP,y0,z0) == 1
                lpoints = lpoints + 1;
                plistX(lpoints,:) = [xP y0 z0];
                if lpoints > poly_num
                    stop_flag = 1;
                end
            else
                stop_flag =1;
            end
        end

        xP = x0;
        stop_flag =0;
        while stop_flag == 0
            xP = xP + 1;
            if xP > DIM(1)
                stop_flag = 1;
            elseif sMASK(xP,y0,z0) == 1
                lpoints = lpoints + 1;
                plistX(lpoints,:) = [xP y0 z0];
                if lpoints > 2*poly_num -1
                    stop_flag = 1;
                end
            else
                stop_flag =1;
            end
        end

        %%%Y Points to be used
        plistY(1,:)= [x0 y0 z0];
        lpoints = 1;
        yP = y0;
        stop_flag = 0;

        while stop_flag == 0
            yP = yP - 1;
            if yP < 1
                stop_flag = 1;
            elseif sMASK(x0,yP,z0) == 1
                lpoints = lpoints + 1;
                plistY(lpoints,:) = [x0 yP z0];
                if lpoints > poly_num
                    stop_flag = 1;
                end
            else
                stop_flag =1;
            end
        end

        yP = y0;
        stop_flag =0;
        while stop_flag == 0
            yP = yP + 1;
            if yP > DIM(2)
                stop_flag = 1;
            elseif sMASK(x0,yP,z0) == 1
                lpoints = lpoints + 1;
                plistY(lpoints,:) = [x0 yP z0];
                if lpoints > 2*poly_num
                    stop_flag = 1;
                end
            else
                stop_flag =1;
            end
        end

        %%%Z Points to be used
        plistZ(1,:)= [x0 y0 z0];
        lpoints = 1;
        zP = z0;
        stop_flag = 0;

        while stop_flag == 0
            zP = zP - 1;
            if zP < 1
                stop_flag = 1;
            elseif sMASK(x0,y0,zP) == 1
                lpoints = lpoints + 1;
                plistZ(lpoints,:) = [x0 y0 zP];
                if lpoints > poly_num
                    stop_flag = 1;
                end
            else
                stop_flag =1;
            end
        end

        zP = z0;
        stop_flag =0;
        while stop_flag == 0
            zP = zP + 1;
            if zP > DIM(3)
                stop_flag = 1;
            elseif sMASK(x0,y0,zP) == 1
                lpoints = lpoints + 1;
                plistZ(lpoints,:) = [x0 y0 zP];
                if lpoints > 2*poly_num
                    stop_flag = 1;
                end
            else
                stop_flag =1;
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

        Xavg = mean(Xvals);
        Yavg = mean(Yvals);
        Zavg = mean(Zvals);

        Xvals_Scale = ( Xvals - Xavg );
        Yvals_Scale = ( Yvals - Yavg );
        Zvals_Scale = ( Zvals - Zavg );

        poly_orderx = length(Xvals) - 1;
        poly_ordery = length(Yvals) - 1;
        poly_orderz = length(Zvals) - 1;
        poly_ordert = 2;

        pow_x = (poly_orderx:-1:0);
        pow_y = (poly_ordery:-1:0);
        pow_z = (poly_orderz:-1:0);
        pow_t = (poly_ordert:-1:0);

        pow_x1 = (pow_x - 1).*( pow_x -1 >= 0);
        pow_y1 = (pow_y - 1).*( pow_y -1 >= 0);
        pow_z1 = (pow_z - 1).*( pow_z -1 >= 0);
        pow_t1 = (pow_t - 1).*( pow_t -1 >= 0);

        pow_x2 = (pow_x - 2).*( pow_x -2 >= 0);
        pow_y2 = (pow_y - 2).*( pow_y -2 >= 0);
        pow_z2 = (pow_z - 2).*( pow_z -2 >= 0);
        pow_t2 = (pow_t - 2).*( pow_t -2 >= 0);

        %%%%%%GET MATRIX OF ENCODING%%%%%%%
        EX = repmat(Xvals_Scale,1,poly_orderx+1).^repmat( (poly_orderx:-1:0),length(Xvals),1);
        EY = repmat(Yvals_Scale,1,poly_ordery+1).^repmat( (poly_ordery:-1:0),length(Yvals),1);
        EZ = repmat(Zvals_Scale,1,poly_orderz+1).^repmat( (poly_orderz:-1:0),length(Zvals),1);

        EXT = EX';
        EYT = EY';
        EZT = EZ';

        EXF = EXT*EX;
        EYF = EYT*EY;
        EZF = EZT*EZ;

        % figure
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

            %%Linsolve based is faster
            Tvals = ( time_start:time_end )';
            ET = repmat(Tvals,1,poly_ordert+1).^repmat( (poly_ordert:-1:0),length(Tvals),1);

            %%%%%Fit values in Time%%%%
            vxt = -conv_vel*VELXt(x0,y0,z0,Tvals);
            vyt = -conv_vel*VELYt(x0,y0,z0,Tvals);
            vzt = -conv_vel*VELZt(x0,y0,z0,Tvals);

            %         fit_vx_dt = polyfit(time_start:time_end',vxt(:)',2);
            %         fit_vy_dt = polyfit(time_start:time_end',vyt(:)',2);
            %         fit_vz_dt = polyfit(time_start:time_end',vzt(:)',2);

            fit_vx_dt = linsolve(ET,vxt(:));
            fit_vy_dt = linsolve(ET,vyt(:));
            fit_vz_dt = linsolve(ET,vzt(:));


            %         clf
            %         plot(Tvals(:),vxt(:),'*')
            %         hold on
            %         plot(time_start:0.01:time_end,polyval(fit_vx_dt,time_start:0.01:time_end))
            %         drawnow

            INDXt = INDX + (time-1)*DIM(1)*DIM(2)*DIM(3);
            INDYt = INDY + (time-1)*DIM(1)*DIM(2)*DIM(3);
            INDZt = INDZ + (time-1)*DIM(1)*DIM(2)*DIM(3);

            %         poly_orderx
            %         poly_ordery
            %         poly_orderz


            if(poly_orderx > 0)
                %%%%%Fit values in X%%%%
                vxV = -conv_vel*VELXt(INDXt);
                %vxF = EXT*vxV;
                fit_vx_dx = linsolve(EX,vxV);

                vyV = -conv_vel*VELYt(INDXt);
                %vyF = EXT*vyV;
                fit_vy_dx = linsolve(EX,vyV);

                vzV = -conv_vel*VELZt(INDXt);
                %vzF = EXT*vzV;
                fit_vz_dx = linsolve(EX,vzV);

            else
                fit_vx_dx = 0;
                fit_vy_dx = 0;
                fit_vz_dx = 0;
            end


            if(poly_ordery > 0)
                %%%%%Fit values in Y%%%%
                vxV = -conv_vel*VELXt(INDYt);
                %vxF = EYT*vxV;
                fit_vx_dy = linsolve(EY,vxV);

                vyV = -conv_vel*VELYt(INDYt);
                %vyF = EYT*vyV;
                fit_vy_dy = linsolve(EY,vyV);

                vzV = -conv_vel*VELZt(INDYt);
                %vzF = EYT*vzV;
                fit_vz_dy = linsolve(EY,vzV);
            else
                fit_vx_dy = 0;
                fit_vy_dy = 0;
                fit_vz_dy = 0;
            end

            if(poly_orderz > 0)
                %%%%%Fit values in X%%%%
                vxV = -conv_vel*VELXt(INDZt);
                %vxF = EZT*vxV;
                fit_vx_dz = linsolve(EZ,vxV);

                vyV = -conv_vel*VELYt(INDZt);
                %vyF = EZT*vyV;
                fit_vy_dz = linsolve(EZ,vyV);

                vzV = -conv_vel*VELZt(INDZt);
                %vzF = EZT*vzV;
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
            vx = -conv_vel*VELXt(x0,y0,z0,time);
            vy = -conv_vel*VELYt(x0,y0,z0,time);
            vz = -conv_vel*VELZt(x0,y0,z0,time);

            %%First Derivatives
            %         dvxdx = polyval( polyder(fit_vx_dx),x0-Xavg)
            %         dvxdy = polyval( polyder(fit_vx_dy),y0-Yavg)
            %         dvxdz = polyval( polyder(fit_vx_dz),z0-Zavg)
            %
            dvxdx = sum( ( (x0 - Xavg).^(pow_x1) ).*pow_x.*fit_vx_dx'.*(pow_x-1 >= 0 ) );
            dvxdy = sum( ( (y0 - Yavg).^(pow_y1) ).*pow_y.*fit_vx_dy'.*(pow_y-1 >= 0 ) );
            dvxdz = sum( ( (z0 - Zavg).^(pow_z1) ).*pow_z.*fit_vx_dz'.*(pow_z-1 >= 0 ) );

            %         dvydx = polyval( polyder(fit_vy_dx),x0-Xavg);
            %         dvydy = polyval( polyder(fit_vy_dy),y0-Yavg);
            %         dvydz = polyval( polyder(fit_vy_dz),z0-Zavg);
            %
            dvydx = sum( ( (x0 - Xavg).^(pow_x1) ).*pow_x.*fit_vy_dx'.*(pow_x-1 >= 0 ) );
            dvydy = sum( ( (y0 - Yavg).^(pow_y1) ).*pow_y.*fit_vy_dy'.*(pow_y-1 >= 0 ) );
            dvydz = sum( ( (z0 - Zavg).^(pow_z1) ).*pow_z.*fit_vy_dz'.*(pow_z-1 >= 0 ) );

            %         dvzdx = polyval( polyder(fit_vz_dx),x0-Xavg);
            %         dvzdy = polyval( polyder(fit_vz_dy),y0-Yavg);
            %         dvzdz = polyval( polyder(fit_vz_dz),z0-Zavg);
            %
            dvzdx = sum( ( (x0 - Xavg).^(pow_x1) ).*pow_x.*fit_vz_dx'.*(pow_x-1 >= 0 ) );
            dvzdy = sum( ( (y0 - Yavg).^(pow_y1) ).*pow_y.*fit_vz_dy'.*(pow_y-1 >= 0 ) );
            dvzdz = sum( ( (z0 - Zavg).^(pow_z1) ).*pow_z.*fit_vz_dz'.*(pow_z-1 >= 0 ) );

            %         dvxdt = polyval( polyder(fit_vx_dt),time);
            %         dvydt = polyval( polyder(fit_vy_dt),time);
            %         dvzdt = polyval( polyder(fit_vz_dt),time);

            dvxdt = sum( ( (time ).^(pow_t1) ).*pow_t.*fit_vx_dt'.*(pow_t-1 >= 0 ) );
            dvydt = sum( ( (time ).^(pow_t1) ).*pow_t.*fit_vy_dt'.*(pow_t-1 >= 0 ) );
            dvzdt = sum( ( (time ).^(pow_t1) ).*pow_t.*fit_vz_dt'.*(pow_t-1 >= 0 ) );


            %%2nd Derivatives
            %        dvxdx2 = polyval( polyder(polyder(fit_vx_dx)),x0-Xavg);
            %        dvxdy2 = polyval( polyder(polyder(fit_vx_dy)),y0-Yavg);
            %        dvxdz2 = polyval( polyder(polyder(fit_vx_dz)),z0-Zavg);
            %

            dvxdx2 = sum( ( (x0 - Xavg).^(pow_x2) ).*pow_x1.*pow_x.*fit_vx_dx'.*(pow_x-2 >= 0 ) );
            dvxdy2 = sum( ( (y0 - Yavg).^(pow_y2) ).*pow_y1.*pow_y.*fit_vx_dy'.*(pow_y-2 >= 0 ) );
            dvxdz2 = sum( ( (z0 - Zavg).^(pow_z2) ).*pow_z1.*pow_z.*fit_vx_dz'.*(pow_z-2 >= 0 ) );

            %         dvydx2 = polyval( polyder(polyder(fit_vy_dx)),x0-Xavg);
            %         dvydy2 = polyval( polyder(polyder(fit_vy_dy)),y0-Yavg);
            %         dvydz2 = polyval( polyder(polyder(fit_vy_dz)),z0-Zavg);

            dvydx2 = sum( ( (x0 - Xavg).^(pow_x2) ).*pow_x1.*pow_x.*fit_vy_dx'.*(pow_x-2 >= 0 ) );
            dvydy2 = sum( ( (y0 - Yavg).^(pow_y2) ).*pow_y1.*pow_y.*fit_vy_dy'.*(pow_y-2 >= 0 ) );
            dvydz2 = sum( ( (z0 - Zavg).^(pow_z2) ).*pow_z1.*pow_z.*fit_vy_dz'.*(pow_z-2 >= 0 ) );

            %         dvzdx2 = polyval( polyder(polyder(fit_vz_dx)),x0-Xavg);
            %         dvzdy2 = polyval( polyder(polyder(fit_vz_dy)),y0-Yavg);
            %         dvzdz2 = polyval( polyder(polyder(fit_vz_dz)),z0-Zavg);

            dvzdx2 = sum( ( (x0 - Xavg).^(pow_x2) ).*pow_x1.*pow_x.*fit_vz_dx'.*(pow_x-2 >= 0 ) );
            dvzdy2 = sum( ( (y0 - Yavg).^(pow_y2) ).*pow_y1.*pow_y.*fit_vz_dy'.*(pow_y-2 >= 0 ) );
            dvzdz2 = sum( ( (z0 - Zavg).^(pow_z2) ).*pow_z1.*pow_z.*fit_vz_dz'.*(pow_z-2 >= 0 ) );


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%   Step 2.4 Navier-Stokes                         %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            GRADx(pos,time) =   -conv_dt*dens*dvxdt ...
                -conv_dx*dens*vx*dvxdx  ...
                -conv_dy*dens*vy*dvxdy  ...
                -conv_dz*dens*vz*dvxdz ...
                +conv_dx2*visc*dvxdx2 ...
                +conv_dy2*visc*dvxdy2 ...
                +conv_dz2*visc*dvxdz2;
            GRADy(pos,time) =  -conv_dt*dens*dvydt  ...
                -conv_dx*dens*vx*dvydx ...
                -conv_dy*dens*vy*dvydy ...
                -conv_dz*dens*vz*dvydz ...
                +conv_dx2*visc*dvydx2 ...
                +conv_dy2*visc*dvydy2 ...
                +conv_dz2*visc*dvydz2;
            GRADz(pos,time) =  -conv_dt*dens*dvzdt    ...
                -conv_dx*dens*vx*dvzdx ...
                -conv_dy*dens*vy*dvzdy ...
                -conv_dz*dens*vz*dvzdz  ...
                +conv_dx2*visc*dvzdx2 ...
                +conv_dy2*visc*dvzdy2 ...
                +conv_dz2*visc*dvzdz2;

        end
    end
    clear plistY
    clear plistX
    clear plistZ
end


function update_image(handles,save_option)

global sMASK;
global GRADx;
global GRADy;
global GRADz;
global PRESSURE;
global plist;
global tframes;
global IMAGE;
global press_axis_plot;
global VELX;
global VELY;
global VELZ;
global VELXt;
global VELYt;
global VELZt;
global tres;

%%Reference Box info
global rbox_handle;
global rbox_exist;
global rphix;
global rphiy;
global rphiz;
global rxpos;
global rypos;
global rzpos;
global rxsize;
global rysize;
global rzsize;
global rbox_alpha;
global rbox_color;

%%Velocity Box Info
global vbox_handle;
global vbox_exist;
global vphix;
global vphiy;
global vphiz;
global vxpos;
global vypos;
global vzpos;
global vxsize;
global vysize;
global vzsize;
global vbox_alpha;
global vbox_color;

clear global slice_data;
global slice_data;

global pressure_ref;
global pressure_vel;
global velocity_ref;
global velocity_vel;

global m_xlength;
global m_ylength;
global m_zlength;

warning('off','MATLAB:divideByZero')

image_object = get(handles.image_object,'Value');
%%% 1 = Box Avg
%%% 2 = Avg in Z

image_type = get(handles.image_type,'Value');
%%% 1 = Time Resolved
%%% 2 = Time Averaged

num_frames = tframes;
if image_type == 2
    num_frames = 1;
end

new_fig = ishandle(press_axis_plot);
if(new_fig==0)
    press_axis = figure;
    scrsz = get(0,'ScreenSize');
    SIZE=[(scrsz(3)*1/2-64) scrsz(4)*1/4 scrsz(3)*1/2 scrsz(4)*1/2];
    set(gcf,'Position',SIZE);
    set(gcf,'Name','Pressure Window');
    daspect([1 1 1])
    set(gca,'color','w');
    set(gcf,'color','w');
    disp('New Figure');
else
    figure(press_axis_plot);
end

if image_object == 1
    pressure_ref =0;
    pressure_vel =0;
    velocity_ref =0;
    velocity_vel =0;

    %%Temp Storage fo Memory Managed Pressure
    for frame =1:num_frames
        PIMAGE = zeros(size(sMASK));
        if num_frames ==1
            PIMAGE(plist) = 0.007501*mean(PRESSURE,2);
        else
            PIMAGE(plist) = 0.007501*PRESSURE(:,frame);
        end

        if(size(rbox_exist,1)~=0)
            %%%%%Find Indices Within the box
            [xB,yB,zB] = meshgrid((-rxsize:rxsize),(-rysize:rysize),(-rzsize:rzsize));

            Rx = [ 1 0        0;
                0 cos(rphix) sin(rphix);
                0 -sin(rphix) cos(rphix)];
            Ry = [ cos(rphiy) 0   -sin(rphiy);
                0         1        0;
                sin(rphiy) 0    cos(rphiy)];
            Rz = [ cos(rphiz) sin(rphiz) 0;
                -sin(rphiz) cos(rphiz) 0;
                0 0 1];

            RT = Rx*Ry*Rz;

            XB = [xB(:) yB(:) zB(:)]';
            XB = ( RT*XB )' + repmat([rxpos rypos rzpos],[numel(XB)/3 1]);;
            
            %%Get a starting point
            XF = floor(XB(:,1));
            YF = floor(XB(:,2));
            ZF = floor(XB(:,3));

            disp(['Point is : ',num2str(XF(1)),'  ',num2str(YF(2)),'  ',num2str(ZF(3))]);
            
            %%find points outside matrix
            MASK = (XF < 1 ) | (YF < 1) | (ZF < 1) | (XF > m_ylength-1) | (YF > m_xlength-1) | (ZF > m_zlength-1);
            idx = find(MASK);
            XF(idx) = 1;
            YF(idx) = 1;
            ZF(idx) = 1;
            
            %%Interpolate sMASK
            idx = sub2ind(size(sMASK),YF,XF,ZF);
            BOX_MASK = sMASK(idx);
                        
            %%Interpolate Velocity
            Vtemp = VELXt(:,:,:,frame);
            BOX_VELOCITY = Vtemp(idx).^2;
            Vtemp = VELYt(:,:,:,frame);
            BOX_VELOCITY(:) = BOX_VELOCITY(:) + Vtemp(idx).^2;
            Vtemp = VELZt(:,:,:,frame);
            BOX_VELOCITY(:) = BOX_VELOCITY(:) + Vtemp(idx).^2;
            BOX_VELOCITY = sqrt(BOX_VELOCITY).*BOX_MASK;
                        
            %%Interpolate Velocity
            BOX_PRESS = PIMAGE(idx);
                                    
            pressure_ref(frame) = sum(BOX_MASK(:).*BOX_PRESS(:))/sum(BOX_MASK(:));
            velocity_ref(frame) = max(BOX_VELOCITY(:));
            
        else
            pressure_ref(frame)=0;
            velocity_ref(frame)=0;
        end%ref box

        if(size(rbox_exist,1)~=0)
            %%%%%Find Indices Within the box
            [xB,yB,zB] = meshgrid((-vxsize:vxsize),(-vysize:vysize),(-vzsize:vzsize));

            Rx = [ 1 0        0;
                0 cos(vphix) sin(vphix);
                0 -sin(vphix) cos(vphix)];
            Ry = [ cos(vphiy) 0   -sin(vphiy);
                0         1        0;
                sin(vphiy) 0    cos(vphiy)];
            Rz = [ cos(vphiz) sin(vphiz) 0;
                -sin(vphiz) cos(vphiz) 0;
                0 0 1];

            RT = Rx*Ry*Rz;

            XB = [xB(:) yB(:) zB(:)]';
            XB = ( RT*XB )' + repmat([vxpos vypos vzpos],[numel(XB)/3 1]);;

            %%Get a starting point
            XF = floor(XB(:,1));
            YF = floor(XB(:,2));
            ZF = floor(XB(:,3));

            %%find points outside matrix
            MASK = (XF < 1 ) | (YF < 1) | (ZF < 1) | (XF > m_ylength-1) | (YF > m_xlength-1) | (ZF > m_zlength-1);
            idx = find(MASK);
            sum(idx)
            XF(idx) = 1;
            YF(idx) = 1;
            ZF(idx) = 1;
            
            %%Interpolate sMASK
            idx = sub2ind(size(sMASK),YF,XF,ZF);
            BOX_MASK = sMASK(idx);
                        
            %%Interpolate Velocity
            Vtemp = VELXt(:,:,:,frame);
            BOX_VELOCITY = Vtemp(idx).^2;
            Vtemp = VELYt(:,:,:,frame);
            BOX_VELOCITY(:) = BOX_VELOCITY(:) + Vtemp(idx).^2;
            Vtemp = VELZt(:,:,:,frame);
            BOX_VELOCITY(:) = BOX_VELOCITY(:) + Vtemp(idx).^2;
            BOX_VELOCITY = sqrt(BOX_VELOCITY).*BOX_MASK;
            
             %%Interpolate Velocity
            BOX_PRESS = PIMAGE(idx);
            
            sum(BOX_MASK(:))
            
            pressure_vel(frame) = sum(BOX_MASK(:).*BOX_PRESS(:))/sum(BOX_MASK(:));
            velocity_vel(frame) = max(BOX_VELOCITY(:));
          
        else
            pressure_vel(frame)=0;
            velocity_vel(frame)=0;
        end%ref box


    end %frames
    velocity_vel
    velocity_ref
    pressure_vel
    pressure_ref
    


    subplot(121)
    plot(tres* (0:num_frames-1),pressure_vel -pressure_ref,'o','linewidth',3);
    hold on;
    plot(tres*(0:0.1:num_frames-1),spline(1:num_frames,pressure_vel -pressure_ref,1:0.1:num_frames),'-','linewidth',3);
    xlabel('Gated Time (ms)','FontSize',18);
    ylabel('\Delta Pressure (mmHg)','FontSize',18);
    set(gca,'FontSize',16);

    subplot(122)
    plot(tres* (0:num_frames-1),velocity_ref,'o','MarkerEdgeColor','r','linewidth',3);
    hold on;
    plot(tres* (0:num_frames-1),velocity_vel,'o','MarkerEdgeColor','b','linewidth',3);
    plot(tres* (0:num_frames-1),velocity_vel-velocity_ref,'o','MarkerEdgeColor','k','linewidth',3);
    legend('Loc 1','Loc 2','Difference');
    plot(tres*(0:0.1:num_frames-1),spline(1:num_frames,velocity_vel,1:0.1:num_frames),'-b','linewidth',3);
    plot(tres*(0:0.1:num_frames-1),spline(1:num_frames,velocity_ref,1:0.1:num_frames),'-r','linewidth',3);
    plot(tres*(0:0.1:num_frames-1),spline(1:num_frames,velocity_vel-velocity_ref,1:0.1:num_frames),'-k','linewidth',3);
    xlabel('Gated Time (ms)','FontSize',18);
    ylabel('Velocity (mm/s)','FontSize',18);
    set(gca,'FontSize',16);

    %%%Update Values
    set(handles.loc1_velocity,'String',num2str(max(velocity_ref)));
    set(handles.loc2_velocity,'String',num2str(max(velocity_vel)));
    set(handles.avg_delP,'String',num2str(mean(pressure_vel-pressure_ref)));
    set(handles.max_delP,'String',num2str(max(abs(pressure_vel-pressure_ref))));
elseif image_object == 2
    %%%%Mip the current slice in Z

    %%Temp Storage fo Memory Managed Pressure
    for frame =1:num_frames
        PIMAGE = zeros(size(sMASK));
        if num_frames ==1
            PIMAGE(plist) = 0.007501*mean(PRESSURE,2);
        else
            PIMAGE(plist) = 0.007501*PRESSURE(:,frame);
        end

        if(size(rbox_exist,1)~=0)
            %%%%%Find Indices Within the box
            [xB,yB,zB] = meshgrid((-rxsize:rxsize),(-rysize:rysize),(-rzsize:rzsize));

            Rx = [ 1 0        0;
                0 cos(rphix) sin(rphix);
                0 -sin(rphix) cos(rphix)];
            Ry = [ cos(rphiy) 0   -sin(rphiy);
                0         1        0;
                sin(rphiy) 0    cos(rphiy)];
            Rz = [ cos(rphiz) sin(rphiz) 0;
                -sin(rphiz) cos(rphiz) 0;
                0 0 1];

            RT = Rx*Ry*Rz;

            XB = [xB(:) yB(:) zB(:)]';
            XB = ( RT*XB )' + repmat([rxpos rypos rzpos],[numel(XB)/3 1]);;

            BOXP = zeros(size(xB));
            for pos=1:numel(xB)
                X = XB(pos,:);
                BOXP(pos) = nn3d(PIMAGE,X(2),X(1),X(3));
            end

            if frame==1
                BOXM = zeros(size(xB));
                for pos=1:numel(xB)
                    X = XB(pos,:);
                    BOXM(pos) = nn3d(sMASK,X(2),X(1),X(3));
                end
            end
            mip_t(:,:,:,frame) = BOXP;
        end


        image = mean(mip_t,4);
        imagem= sum(image,3)./sum(BOXM,3);

        clf
        imagesc(imagem)
    end




end



return


if image_object > 7
    figure(press_axis);
    clf
    hpatch = patch(isosurface(sMASK,0.5));
    reducepatch(hpatch,0.1);
    set(hpatch,'FaceColor','red','EdgeColor', 'none');
    alpha(0.2);
    set(press_axis, 'Renderer','OpenGL')
    set(press_axis, 'RendererMode','Manual');
    set(gca,'color','black');
    set(gcf,'color','black');
    daspect([1 1 1])
else


    IMAGE = zeros(size(sMASK));
    DIM=size(sMASK);
    VENC = 600;

    cmap = jet(256);
    for t=1:num_frames

        if image_object == 1
            if num_frames == 1
                IMAGE(plist) = mean(GRADx,2);
            else
                IMAGE(plist)=GRADx(:,t);
            end
            min_im = min(IMAGE(:));
            max_im = max(IMAGE(:));
        elseif image_object == 2
            if num_frames == 1
                IMAGE(plist) = mean(GRADy,2);
            else
                IMAGE(plist)=GRADy(:,t);
            end
            min_im = min(IMAGE(:));
            max_im = max(IMAGE(:));
        elseif image_object == 3
            if num_frames == 1
                IMAGE(plist) = mean(GRADz,2);
            else
                IMAGE(plist)=GRADz(:,t);
            end
            min_im = min(IMAGE(:));
            max_im = max(IMAGE(:));
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
            if num_frames == 1
                IMAGE(plist)=0.007501*mean(PRESSURE,2);
                min_im = min(IMAGE(:))-1;
                max_im = max(IMAGE(:));
            else
                IMAGE(plist)=0.007501*PRESSURE(:,t);
                min_im = 0.007501*min(PRESSURE(:))-1;
                max_im = 0.007501*max(PRESSURE(:));
            end

            cmap(1,:)=0;
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
        slice_data(:,:,t)=reshape(sum(IMAGE,image_dir),[size(IMAGE,d1) size(IMAGE,d2)])./reshape(sum(sMASK,image_dir),[size(IMAGE,d1) size(IMAGE,d2)]);
        colormap(cmap);
        imagesc(slice_data(:,:,t),[min_im max_im])
        colorbar('FontSize',22)
        set(gca,'XTick',[],'YTick',[]);
        daspect([1 1 1]);
        pause(0.1);

        if save_option == 1
            base_name = get(handles.save_name,'String')
            fname=[base_name,int2str(t),'.jpg'];
            print('-opengl','-f1','-r200','-djpeg100',fname);
            saveas(press_axis,[base_name,'.fig']);
        end


    end

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
update_image(handles,0);


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
global POINTS_FOUND;

signVx = 1;
signVy = 1;
signVz = 1;

[size_x size_y size_z] = size(sMASK);

max_iter = str2double(get(handles.max_iter,'String'));
alpha = str2double(get(handles.alpha,'String'));
max_error = str2double(get(handles.iter_error,'String'));

%%%%%%Step 1. Get Neighbors of Each Pixel %%%%%%%
%if POINTS_FOUND == 0
if 1==1
    POINTS_FOUND = 1;
    points = size(plist,1);

    [xlist ylist zlist] = ind2sub(size(sMASK),plist);
    nbs=zeros(points,6);

    sMASK = sMASK.*zeros(size(sMASK));
    sMASK(plist) = 1:points;

    for pos = 1:points

        x0 = xlist(pos);
        y0 = ylist(pos);
        z0 = zlist(pos);

        if(x0 < size_x)
            if( sMASK(x0+1,y0,z0)>0)
                nbs(pos,1) = sMASK(x0+1,y0,z0);
            end
        end

        if(x0 > 1)
            if( sMASK(x0-1,y0,z0)>0)
                nbs(pos,2) = sMASK(x0-1,y0,z0);
            end
        end

        if(y0 < size_y)
            if( sMASK(x0,y0+1,z0)>0)
                nbs(pos,3) = sMASK(x0,y0+1,z0);
            end
        end

        if(y0 > 1)
            if( sMASK(x0,y0-1,z0)>0)
                nbs(pos,4) = sMASK(x0,y0-1,z0);
            end
        end

        if(z0 < size_z)
            if( sMASK(x0,y0,z0+1)>0)
                nbs(pos,5) = sMASK(x0,y0,z0+1);
            end
        end

        if(z0 > 1)
            if( sMASK(x0,y0,z0-1)>0)
                nbs(pos,6) = sMASK(x0,y0,z0-1);
            end
        end
    end

    sMASK = sMASK > 0;
end


%%%%%%Step 2. Initial Setup (Simple Integration) %%%%%%%
num_counted = 0;
counted = zeros(points,1);
start_cnt = 1;
new_points =1;
PRESSURE = zeros(points,tframes);
regions = 0;

while num_counted < points

    start_pos =  find(counted==0,1);
    num_counted = num_counted + 1;
    counted(start_pos) = 1;
    PRESSURE(start_pos)= 0;
    cnt_idx(num_counted)=start_pos;
    new_points = 100;

    regions = regions + 1;
    region_idx(regions) = start_pos;
    disp(['Start Region Number ',int2str(regions)]);

    while new_points > 0

        num_countedp = num_counted;
        for pos = start_pos:num_counted

            idx = cnt_idx(pos);

            xp = nbs(idx,1);
            xm = nbs(idx,2);
            yp = nbs(idx,3);
            ym = nbs(idx,4);
            zp = nbs(idx,5);
            zm = nbs(idx,6);

            if xp ~= 0 && counted(xp)==0
                counted(xp) = 1;
                num_counted = num_counted+1;
                cnt_idx(num_counted)=xp;
                PRESSURE(xp,:) = PRESSURE(idx,:) +signVx*GRADx(idx,:)*delX/1000;
            end

            if xm ~= 0 && counted(xm)==0
                counted(xm) = 1;
                num_counted = num_counted+1;
                cnt_idx(num_counted)=xm;
                PRESSURE(xm,:) = PRESSURE(idx,:) - signVx*GRADx(idx,:)*delX/1000;
            end

            if yp ~= 0 &&  counted(yp)==0
                counted(yp) = 1;
                num_counted = num_counted+1;
                cnt_idx(num_counted)=yp;
                PRESSURE(yp,:) = PRESSURE(idx,:) +signVy*GRADy(idx,:)*delY/1000;
            end

            if ym ~= 0 && counted(ym)==0
                counted(ym) = 1;
                num_counted = num_counted+1;
                cnt_idx(num_counted)=ym;
                PRESSURE(ym,:) = PRESSURE(idx,:) - signVy*GRADy(idx,:)*delY/1000;
            end

            if zp ~= 0 && counted(zp)==0
                counted(zp) = 1;
                num_counted = num_counted+1;
                cnt_idx(num_counted)=zp;
                PRESSURE(zp,:) = PRESSURE(idx,:) + signVz*GRADz(idx,:)*delZ/1000;
            end

            if zm ~= 0 && counted(zm)==0
                counted(zm) = 1;
                num_counted = num_counted+1;
                cnt_idx(num_counted)=zm;
                PRESSURE(zm,:) = PRESSURE(idx,:) - signVz*GRADz(idx,:)*delZ/1000;
            end
        end


        new_points = num_counted - num_countedp;

        %           if(verbose)
        %             disp(['Points Counted ',num2str(num_counted),' of ',num2str(points)]);;
        %           end
        %
    end
end


%%%%%%%%%%%Step 3 ITERATE TO CLEAN UP
PRESSURE_old = PRESSURE;
errord = 1e99;
max_error
iter = 0;
while errord > max_error && iter < max_iter
%if 1==0
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
            GRAD_TERM = GRAD_TERM + PRESSURE_old(xp,:) - signVx*GRADx(idx,:)*delX/1000;
            nb_count = nb_count + 1;
        end

        if xm ~= 0
            GRAD_TERM =  GRAD_TERM + PRESSURE_old(xm,:) + signVx*GRADx(idx,:)*delX/1000;
            nb_count = nb_count + 1;
        end

        if yp ~= 0
            GRAD_TERM =  GRAD_TERM + PRESSURE_old(yp,:) - signVy*GRADy(idx,:)*delY/1000;
            nb_count = nb_count + 1;
        end

        if ym ~= 0
            GRAD_TERM =  GRAD_TERM + PRESSURE_old(ym,:) + signVy*GRADy(idx,:)*delY/1000;
            nb_count = nb_count + 1;
        end

        if zp ~= 0
            GRAD_TERM =  GRAD_TERM + PRESSURE_old(zp,:) - signVz*GRADz(idx,:)*delZ/1000;
            nb_count = nb_count + 1;
        end

        if zm ~= 0
            GRAD_TERM =  GRAD_TERM + PRESSURE_old(zm,:) + signVz*GRADz(idx,:)*delZ/1000;
            nb_count = nb_count + 1;
        end

        if nb_count >0
            PRESSURE(idx,:) = (1- alpha)*PRESSURE_old(idx,:) + alpha*GRAD_TERM/nb_count;
        end
    end

    iter = iter +1;
    %errord = sum( sqrt( (PRESSURE(:)- PRESSURE_old(:)).^2) )/numel(PRESSURE);
    errord = mean( abs(PRESSURE(:)- PRESSURE_old(:)));
    disp(['Iter ',num2str(iter),'  Error ',num2str(errord)])
end





function iter_error_Callback(hObject, eventdata, handles)
% hObject    handle to iter_error (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iter_error as text
%        str2double(get(hObject,'String')) returns contents of iter_error as a double


% --- Executes during object creation, after setting all properties.
function iter_error_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iter_error (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on slider movement.
function xmarker_slider_Callback(hObject, eventdata, handles)
% hObject    handle to xmarker_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
update_marker(handles);

% --- Executes during object creation, after setting all properties.
function xmarker_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xmarker_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function ymarker_slider_Callback(hObject, eventdata, handles)
% hObject    handle to ymarker_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

update_marker(handles);
% --- Executes during object creation, after setting all properties.
function ymarker_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ymarker_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function zmarker_slider_Callback(hObject, eventdata, handles)
% hObject    handle to zmarker_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
update_marker(handles);

% --- Executes during object creation, after setting all properties.
function zmarker_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zmarker_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end




function update_marker(handles)

global m_xlength;
global m_ylength;
global m_zlength;
global press_axis;
global marker_handle;
global marker_handle_val;

global xpos_ref;
global ypos_ref;
global zpos_ref;

global xpos_val;
global ypos_val;
global zpos_val;

xpos_ref = floor(1 + (get(handles.xmarker_slider,'Value')*( m_ylength-1) ));
ypos_ref = floor(1 + (get(handles.ymarker_slider,'Value')*( m_xlength-1) ));
zpos_ref = floor(1 + (get(handles.zmarker_slider,'Value')*( m_zlength-1) ));

xpos_val = floor(1 + (get(handles.xmarker_slider_val,'Value')*( m_ylength-1) ));
ypos_val = floor(1 + (get(handles.ymarker_slider_val,'Value')*( m_xlength-1) ));
zpos_val = floor(1 + (get(handles.zmarker_slider_val,'Value')*( m_zlength-1) ));

image_dir = get(handles.image_dir,'Value');
%%% 1 = X
%%% 2 = Y
%%% 3 = Z

image_object = get(handles.image_object,'Value');

figure(press_axis);
if ishandle(marker_handle)
    delete(marker_handle);
end

if ishandle(marker_handle_val)
    delete(marker_handle_val);
end

if( image_object > 7)
    hold on
    marker_handle = plot3(ypos_ref,xpos_ref,zpos_ref,'x','color','w','MarkerSize',16,'linewidth',5);
    marker_handle_val = plot3(ypos_val,xpos_val,zpos_val,'x','color','y','MarkerSize',16,'linewidth',5);
    hold off
else
    if image_dir ==1
        hold on
        marker_handle = plot(zpos_ref,ypos_ref,'x','color','w','MarkerSize',16,'linewidth',5);
        marker_handle_val = plot(zpos_val,ypos_val,'x','color','k','MarkerSize',16,'linewidth',5);
        hold off
    elseif image_dir ==2
        hold on
        marker_handle = plot(zpos_ref,xpos_ref,'x','color','w','MarkerSize',16,'linewidth',5);
        marker_handle_val = plot(zpos_val,xpos_val,'x','color','k','MarkerSize',16,'linewidth',5);
        hold off
    elseif image_dir == 3
        hold on
        marker_handle = plot(ypos_ref,xpos_ref,'x','color','w','MarkerSize',16',linewidth',5);
        marker_handle_val = plot(ypos_val,xpos_val,'x','color','k','MarkerSize',16,'linewidth',5);
        hold off
    end
end

% --- Executes on button press in press_plot.
function press_plot_Callback(hObject, eventdata, handles)
% hObject    handle to press_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

update_plot(handles,0);


function update_plot(handles,save_option)

global xpos_ref;
global ypos_ref;
global zpos_ref;
global xpos_val;
global ypos_val;
global zpos_val;
global press_axis_plot;
global slice_data;

global tframes;
global tres;

global VELX;
global VELY;
global VELZ;
global VELXt;
global VELYt;
global VELZt;

global plist;
global xlist;
global ylist;
global zlist;
global PRESSURE;
image_dir = get(handles.image_dir,'Value');
%%% 1 = X
%%% 2 = Y
%%% 3 = Z



new_fig = ishandle(press_axis_plot)
figure(press_axis_plot)

daspect([1 1 1])

if(new_fig==0)
    scrsz = get(0,'ScreenSize');
    SIZE=[(scrsz(3)*1/2-64) scrsz(4)*1/4 scrsz(3)*1/2 scrsz(4)*1/2];
    set(gcf,'Position',SIZE);
    set(gcf,'Name','Pressure Plot Window');
end

image_object = get(handles.image_object,'Value');
%%% 1 = Gradient X
%%% 2 = Gradient Y
%%% 3 = Gradient Z
%%% 4 = Velocity X
%%% 5 = Velocity Y
%%% 6 = Velocity Z
%%% 7 = Reference Box
%%% 8 = Measurement Box

if image_object == 7
    if image_dir ==1
        plot(tres* (0:tframes-1),reshape(slice_data(ypos_val,zpos_val,:)-slice_data(ypos_ref,zpos_ref,:),[tframes 1 1]),'o','linewidth',3);
        hold on
        plot(tres* (0:0.1:tframes-1),spline(tres* (0:tframes-1),reshape(slice_data(ypos_val,zpos_val,:)-slice_data(ypos_ref,zpos_ref,:),[tframes 1 1]),tres* (0:0.1:tframes-1)),'linewidth',3);
        hold off
        xlabel('Gated Time (ms)','FontSize',18);
        ylabel('Pressure Difference (mmHg)','FontSize',18);
        set(gca,'FontSize',16);
    elseif image_dir ==2
        plot(tres* (0:tframes-1),reshape(slice_data(xpos_val,zpos_val,:)-slice_data(xpos_ref,zpos_ref,:),[tframes 1 1]),'o','linewidth',3);
        hold on
        plot(tres* (0:0.1:tframes-1),spline(tres* (0:tframes-1),reshape(slice_data(xpos_val,zpos_val,:)-slice_data(xpos_ref,zpos_ref,:),[tframes 1 1]),tres* (0:0.1:tframes-1)),'linewidth',3);
        hold off
        xlabel('Gated Time (ms)','FontSize',18);
        ylabel('Pressure Difference (mmHg)','FontSize',18);
        set(gca,'FontSize',16);

    elseif image_dir == 3
        plot(tres* (0:tframes-1),reshape(slice_data(ypos_val,xpos_val,:)-slice_data(ypos_ref,xpos_ref,:),[tframes 1 1]),'o','linewidth',3);
        hold on
        plot(tres* (0:0.1:tframes-1),spline(tres* (0:tframes-1),reshape(slice_data(xpos_val,ypos_val,:)-slice_data(xpos_ref,ypos_ref,:),[tframes 1 1]),tres* (0:0.1:tframes-1)),'linewidth',3);
        hold off
        xlabel('Gated Time (ms)','FontSize',18);
        ylabel('Pressure Difference (mmHg)','FontSize',18);
        set(gca,'FontSize',16);
    end
elseif image_object < 4
    if image_dir ==1
        plot(tres* (0:tframes-1),reshape(slice_data(ypos_ref,zpos_ref,:),[tframes 1 1]),'b','linewidth',3);
        hold on
        plot(tres* (0:tframes-1),reshape(slice_data(ypos_val,zpos_val,:),[tframes 1 1]),'k','linewidth',3);
        legend('Reference Point','Eval Point','FontSize',16)
        xlabel('Gated Time (ms)','FontSize',18);
        ylabel('Gradient (Pa/m)','FontSize',18);
        set(gca,'FontSize',16);
        hold off
    elseif image_dir ==2
        plot(tres* (0:tframes-1),reshape(slice_data(xpos_ref,xpos_ref,:),[tframes 1 1]),'b','linewidth',3);
        hold on
        plot(tres* (0:tframes-1),reshape(slice_data(xpos_val,zpos_val,:),[tframes 1 1]),'k','linewidth',3);
        legend('Reference Point','Eval Point','FontSize',16)
        xlabel('Gated Time (ms)','FontSize',18);
        ylabel('Gradient (Pa/m)','FontSize',18);
        set(gca,'FontSize',16);
        hold off
    elseif image_dir == 3
        plot(tres* (0:tframes-1),reshape(slice_data(ypos_ref,xpos_ref,:),[tframes 1 1]),'b','linewidth',3);
        hold on
        plot(tres* (0:tframes-1),reshape(slice_data(ypos_val,xpos_val,:),[tframes 1 1]),'k','linewidth',3);
        legend('Reference Point','Eval Point','FontSize',16)
        xlabel('Gated Time (ms)','FontSize',18);
        ylabel('Gradient (Pa/m)','FontSize',18);
        set(gca,'FontSize',16);
        hold off
    end
elseif image_object == 8

    S=1;
    SPD =sqrt( VELXt(xpos_ref-S:xpos_ref+S,ypos_ref-S:ypos_ref+S,zpos_ref-S:zpos_ref+S,:).^2 + ...
        VELYt(xpos_ref-S:xpos_ref+S,ypos_ref-S:ypos_ref+S,zpos_ref-S:zpos_ref+S,:).^2 + ...
        VELZt(xpos_ref-S:xpos_ref+S,ypos_ref-S:ypos_ref+S,zpos_ref-S:zpos_ref+S,:).^2 );
    max_SPD = reshape(  max(max(max(SPD,[],1),[],2),[],3), [tframes 1 1 1]);

    SPD =sqrt(VELXt(xpos_val-S:xpos_val+S,ypos_val-S:ypos_val+S,zpos_val-S:zpos_val+S,:).^2 + ...
        VELYt(xpos_val-S:xpos_val+S,ypos_val-S:ypos_val+S,zpos_val-S:zpos_val+S,:).^2 + ...
        VELZt(xpos_val-S:xpos_val+S,ypos_val-S:ypos_val+S,zpos_val-S:zpos_val+S,:).^2 );
    max_SPDv = reshape(  max(max(max(SPD,[],1),[],2),[],3), [tframes 1 1 1]);


    plot(tres* (0:tframes-1),max_SPD,'o','linewidth',3);
    hold on
    plot(tres* (0:tframes-1),max_SPDv,'o','linewidth',3,'color','r');
    plot(tres* (0:tframes-1),max_SPDv-max_SPD,'o','linewidth',3,'color','k');


    legend('Reference Point','Eval Point','Difference','FontSize',16)

    plot(tres* (0:0.1:tframes-1),spline(tres* (0:tframes-1),max_SPDv,tres* (0:0.1:tframes-1)),'linewidth',3,'color','r');
    plot(tres* (0:0.1:tframes-1),spline(tres* (0:tframes-1),max_SPDv-max_SPD,tres* (0:0.1:tframes-1)),'linewidth',3,'color','k');
    plot(tres* (0:0.1:tframes-1),spline(tres* (0:tframes-1),max_SPD,tres* (0:0.1:tframes-1)),'linewidth',3);

    max_vel_ref = max( spline(tres* (0:tframes-1),max_SPD,tres* (0:0.1:tframes-1)));
    max_vel_val = max( spline(tres* (0:tframes-1),max_SPDv,tres* (0:0.1:tframes-1)));
    set(handles.max_vel_ref,'String',num2str(max_vel_ref));
    set(handles.max_vel_val,'String',num2str(max_vel_val));

    xlabel('Gated Time (ms)','FontSize',18);
    ylabel('Velocity (mm/s)','FontSize',18);
    set(gca,'FontSize',16);
    hold off
elseif image_object == 9

    idx_ref = find( (xpos_ref == xlist).*( ypos_ref == ylist).*(zpos_ref == zlist ));
    idx_val = find( (xpos_val == xlist).*( ypos_val == ylist).*(zpos_val == zlist ));

    size(idx_ref )
    size(idx_val )

    if  ( size(idx_ref,1) == 0 ) && ( size(idx_val,1) == 0 )
        title('Error Outside of Point List');
        xlabel('');
        ylabel('');
    else
        disp('plot')
        plot(tres* (0:tframes-1),0.007501*( PRESSURE(idx_val,:)-PRESSURE(idx_ref,:)),'o','linewidth',3,'color','k');
        hold on
        plot(tres* (0:0.1:tframes-1),spline(tres* (0:tframes-1),0.007501*(PRESSURE(idx_val,:)-PRESSURE(idx_ref,:)),tres* (0:0.1:tframes-1)),'linewidth',3,'color','k');

        max_del_press = max(abs( spline(tres* (0:tframes-1),0.007501*(PRESSURE(idx_val,:)-PRESSURE(idx_ref,:)),tres* (0:0.1:tframes-1))));
        avg_del_press = mean(( spline(tres* (0:tframes-1),0.007501*(PRESSURE(idx_val,:)-PRESSURE(idx_ref,:)),tres* (0:0.1:tframes-1))));
        set(handles.max_del_press,'String',num2str(max_del_press));
        set(handles.avg_del_press,'String',num2str(avg_del_press));

        xlabel('Gated Time (ms)','FontSize',18);
        ylabel('Pressure (mmHg)','FontSize',18);
        set(gca,'FontSize',16);
        hold off
    end


end

if (save_option == 1)

    base_name = get(handles.save_name,'String')
    fname=[base_name,'.jpg'];
    print('-opengl','-f2','-r200','-djpeg100',fname);
    saveas(press_axis_plot,[base_name,'.fig']);
end



% --- Executes on slider movement.
function xmarker_slider_val_Callback(hObject, eventdata, handles)
% hObject    handle to xmarker_slider_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

update_marker(handles);
% --- Executes during object creation, after setting all properties.
function xmarker_slider_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xmarker_slider_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function ymarker_slider_val_Callback(hObject, eventdata, handles)
% hObject    handle to ymarker_slider_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

update_marker(handles);

% --- Executes during object creation, after setting all properties.
function ymarker_slider_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ymarker_slider_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function zmarker_slider_val_Callback(hObject, eventdata, handles)
% hObject    handle to zmarker_slider_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
update_marker(handles);

% --- Executes during object creation, after setting all properties.
function zmarker_slider_val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zmarker_slider_val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
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


% --- Executes on button press in save_button.
function save_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

base_name = get(handles.save_name,'string');

box_name = [base_name,'.box_vals'];
mip_name = [base_name,'.mip'];

dlmwrite(hist_name,[hist_xout; hist_nout]','\t');
dlmwrite(raw_name,[1:length(box_idx); wss_out]','\t');




% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

update_image(handles,1);




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



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
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


% --- Executes on slider movement.
function box_xslice_Callback(hObject, eventdata, handles)
% hObject    handle to box_xslice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
update_box(handles);

% --- Executes during object creation, after setting all properties.
function box_xslice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to box_xslice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in box_alpha.
function box_alpha_Callback(hObject, eventdata, handles)
% hObject    handle to box_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns box_alpha contents as cell array
%        contents{get(hObject,'Value')} returns selected item from box_alpha
update_box(handles);

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
update_box(handles);

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
update_box(handles);

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
update_box(handles);

% --- Executes during object creation, after setting all properties.
function box_zrot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to box_zrot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


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
global press_axis;
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
set(handles.box_xslide,'Value',round(COMx-1)/( m_ylength-1));
set(handles.box_yslide,'Value',round(COMy-1)/( m_xlength-1));
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
global press_axis;
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


% --- Executes on selection change in box_type.
function box_type_Callback(hObject, eventdata, handles)
% hObject    handle to box_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns box_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from box_type

%%Reference Box info
global rbox_handle;
global rbox_exist;
global rphix;
global rphiy;
global rphiz;
global rxpos;
global rypos;
global rzpos;
global rxsize;
global rysize;
global rzsize;
global rbox_alpha;
global rbox_color;

%%Velocity Box Info
global vbox_handle;
global vbox_exist;
global vphix;
global vphiy;
global vphiz;
global vxpos;
global vypos;
global vzpos;
global vxsize;
global vysize;
global vzsize;
global vbox_alpha;
global vbox_color;

global m_xlength;
global m_ylength;
global m_zlength;


box_type = get(handles.box_type,'Value');
% 1 = Turn off Boxes
% 2 = Set Reference Box
% 3 = Set Measure Box

if box_type ==2 && size(rbox_exist,1)~=0
    set(handles.box_xrot,'Value',(rphix+pi/2)/(pi));
    set(handles.box_yrot,'Value',(rphiy+pi/2)/(pi));
    set(handles.box_zrot,'Value',(rphiz+pi/2)/(pi));
    set(handles.box_xslide,'Value',round(rxpos-1)/( m_ylength-1));
    set(handles.box_yslide,'Value',round(rypos-1)/( m_xlength-1));
    set(handles.box_zslide,'Value',round(rzpos-1)/( m_zlength-1));
    set(handles.box_xsize,'String',num2str(rxsize));
    set(handles.box_ysize,'String',num2str(rysize));
    set(handles.box_zsize,'String',num2str(rzsize));
    set(handles.box_alpha,'Value',(1-rbox_alpha)/0.25 + 1);
    set(handles.box_color,'Value',rbox_color);

elseif size(vbox_exist,1)~=0
    set(handles.box_xrot,'Value',(vphix+pi/2)/(pi));
    set(handles.box_yrot,'Value',(vphiy+pi/2)/(pi));
    set(handles.box_zrot,'Value',(vphiz+pi/2)/(pi));
    set(handles.box_xslide,'Value',round(vxpos-1)/( m_ylength-1));
    set(handles.box_yslide,'Value',round(vypos-1)/( m_xlength-1));
    set(handles.box_zslide,'Value',round(vzpos-1)/( m_zlength-1));
    set(handles.box_xsize,'String',num2str(vxsize));
    set(handles.box_ysize,'String',num2str(vysize));
    set(handles.box_zsize,'String',num2str(vzsize));
    set(handles.box_alpha,'Value',(1-vbox_alpha)/0.25 + 1);
    set(handles.box_color,'Value',vbox_color);

end
update_box(handles);

% --- Executes during object creation, after setting all properties.
function box_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to box_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
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




% --- Executes on slider movement.
function box_xslide_Callback(hObject, eventdata, handles)
% hObject    handle to box_xslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

update_box(handles);



function box_xsize_Callback(hObject, eventdata, handles)
% hObject    handle to box_xsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of box_xsize as text
%        str2double(get(hObject,'String')) returns contents of box_xsize as a double
update_box(handles);


function update_box(handles)

%%Anatomical Data
global sMAG;
global sCD;
global sMASK;

%%Visualization Stuff
global vis_axis;
global vis_alpha;
global press_axis;
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
pol =0;

%%Reference Box info
global rbox_handle;
global rphix;
global rphiy;
global rphiz;
global rxpos;
global rypos;
global rzpos;
global rxsize;
global rysize;
global rzsize;
global rbox_exist;
global rbox_alpha;
global rbox_color;

%%Velocity Box Info
global vbox_handle;
global vphix;
global vphiy;
global vphiz;
global vxpos;
global vypos;
global vzpos;
global vxsize;
global vysize;
global vzsize;
global vbox_exist;
global vbox_alpha;
global vbox_color;

%%%%%%Delete Object or Not %%%%%%%%%%%%%%
new_fig = ishandle(press_axis);
if new_fig==0
    press_axis = figure;
    clf
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

    hpatch = patch(isosurface(sMASK,vis_thresh));
    colormap('jet');
    reducepatch(hpatch,0.4);
    set(hpatch,'FaceColor','red','EdgeColor', 'none');
    isonormals(sCD,hpatch);

    camlight right;
    lighting gouraud
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
else
    figure(press_axis);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Determine which Box   %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

box_type = get(handles.box_type,'Value')
% 1 = Turn off Boxes
% 2 = Set Reference Box
% 3 = Set Measure Box

if box_type ==1
    if ishandle(rbox_handle)
        delete(rbox_handle);
    end
    if ishandle(vbox_handle)
        delete(vbox_handle);
    end
elseif box_type ==2
    if ishandle(rbox_handle)
        delete(rbox_handle);
    end
else
    if ishandle(vbox_handle)
        delete(vbox_handle);
    end
end

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
if box_type ~=1
    box_handle = tetramesh(Tes,XB,'FaceAlpha',box_alpha,'FaceColor',col,'EdgeColor','none');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Update the Figures in the GUI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if box_type ==2
    rbox_handle=box_handle;
    rbox_exist =1;
    rphix=phix;
    rphiy=phiy;
    rphiz=phiz;
    rxpos=xpos;
    rypos=ypos;
    rzpos=zpos;
    rxsize=xsize;
    rysize=ysize;
    rzsize=zsize;
    rbox_alpha = box_alpha;
    rbox_color = box_color;
else
    vbox_handle=box_handle;
    vbox_exist =1;
    vphix=phix;
    vphiy=phiy;
    vphiz=phiz;
    vxpos=xpos;
    vypos=ypos;
    vzpos=zpos;
    vxsize=xsize;
    vysize=ysize;
    vzsize=zsize;
    vbox_alpha = box_alpha;
    vbox_color = box_color;
end
