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

    % Last Modified by GUIDE v2.5 20-Oct-2019 00:20:29

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

    global zoomed zoomedAnat spanUD spanLR spanUDAnat spanLRAnat
    global  interpType requestInterp
    handles.anatDatasets = varargin{1};
    handles.pcDatasets = varargin{2};
    handles.magDatasets = varargin{3};
    handles.pcDatasets(1).ROI = [];
    handles.magDatasets(1).ROI = [];
    handles.anatDatasets(1).Centerline = [];
    handles.pcDatasets(1).ROIdata = [];
    handles.pcDatasets(1).ROIdataMakima = [];
    handles.pcDatasets(1).ROIdataPCHIP = [];
    handles.pcDatasets(1).ROIdataSpline = [];
    
    zoomed = 0;
    zoomedAnat = 0;
    spanUD = [0,0];
    spanLR = [0,0];
    spanUDAnat = [0,0];
    spanLRAnat = [0,0];
    interpType = 'None';
    requestInterp = 0;
    
    set(handles.PlanePopup,'String',{handles.pcDatasets.Names});
    set(handles.DatasetPopup,'String',fieldnames(handles.pcDatasets(1).Data));
    set(handles.InterpolatePopup,'String',{'None','Makima','PCHIP','Spline'});
    set(handles.AnatListbox,'String',{handles.anatDatasets.Names});
    set(handles.InterpolatePopup,'Enable','off');
    if numel(handles.magDatasets)==0
        set(handles.MAGRadio,'Enable','off')
    end 
    
    set(handles.ShowPlanesRadio,'Enable','off');
    set(handles.ZoomAnatomicalButton,'Enable','off');
    set(handles.UnzoomAnatomicalButton,'Enable','off');
    set(handles.DrawCenterlineButton,'Enable','off');
    set(handles.DeleteCenterlineButton,'Enable','off');
    set(handles.DeleteCenterlineButton,'Enable','off');
    set(handles.ComputePWVButton,'Enable','off');
    
    guidata(hObject, handles);
    updateImages(handles);
    updateAnatImages(handles);

% --- Outputs from this function are returned to the command line.
function varargout = AnalyzePWV_OutputFcn(hObject, eventdata, handles) 
    varargout{1} = handles.output;


    
    
%%%%%%%%%%%% LOAD 2DPC PLANE %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- PLANE PLOT - CREATE FUNCTION
function PlanePlot_CreateFcn(hObject, eventdata, handles)


% --- PLANE SLIDER - CALLBACK
function PlaneSlider_Callback(hObject, eventdata, handles)
    updateImages(handles);

% --- PLANE SLIDER - CREATE FUNCTION
function PlaneSlider_CreateFcn(hObject, eventdata, handles)
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end

    
% --- PC RADIO - CALLBACK
function PCRadio_Callback(hObject, eventdata, handles)
    set(handles.PlanePopup,'String',{handles.pcDatasets.Names});
    set(handles.DatasetPopup,'Enable','on');
    set(handles.DatasetPopup,'String',fieldnames(handles.pcDatasets(1).Data));
    set(handles.MAGRadio,'Value',0);
    set(handles.PlanePopup,'Value',1);
    set(handles.DatasetPopup,'Value',1);
    updateImages(handles);

    
% --- MAG RADIO - CALLBACK
function MAGRadio_Callback(hObject, eventdata, handles)
    set(handles.PlanePopup,'String',{handles.magDatasets.Names});
    set(handles.DatasetPopup,'Enable','off');
    set(handles.PCRadio,'Value',0);
    set(handles.PlanePopup,'Value',1);
    set(handles.DatasetPopup,'Value',1);
    set(handles.DatasetPopup,'String',' ');
    updateImages(handles);
    

% --- PLANE DROPDOWN - CALLBACK
function PlanePopup_Callback(hObject, eventdata, handles)
    global zoomed spanUD spanLR
    zoomed = 0; spanUD = 0; spanLR = 0;
    planeNum = get(handles.PlanePopup,'Value');
    if get(handles.PCRadio,'Value')
        set(handles.DatasetPopup,'String',fieldnames(handles.pcDatasets(planeNum).Data));
    end 
    set(handles.DatasetPopup,'Value',1);
    updateImages(handles);

% --- PLANE DROPDOWN - CREATE FUNCTION
function PlanePopup_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

    
% --- DATASET DROPDOWN - CALLBACK
function DatasetPopup_Callback(hObject, eventdata, handles)
    updateImages(handles);

% --- DATASET DROPDOWN - CREATE FUNCTION
function DatasetPopup_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    
    
% --- ZOOM BUTTON - CALLBACK
function ZoomButton_Callback(hObject, eventdata, handles)
    global zoomed spanUD spanLR
    disp('Draw an ROI to zoom in');
    rect = drawrectangle;
    positions = round(rect.Position);
    spanUD = spanUD(1)+positions(2):(spanUD(1)+positions(2)+positions(4));
    spanLR = spanLR(1)+positions(1):(spanLR(1)+positions(1)+positions(3));
    zoomed = 1;
    updateImages(handles);
 
    
% --- UNZOOM BUTTON - CALLBACK
function UnzoomButton_Callback(hObject, eventdata, handles)
    global zoomed spanUD spanLR
    zoomed = 0; spanUD = 0; spanLR = 0;
    updateImages(handles);

    
% --- DRAWROI BUTTON - CALLBACK
function DrawROIbutton_Callback(hObject, eventdata, handles)
    planeNum = get(handles.PlanePopup,'Value');
    axes(handles.PlanePlot);
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
    
    if get(handles.PCRadio,'Value')
        set(handles.MAGRadio,'Enable','off');
    else 
        set(handles.PCRadio,'Enable','off');
    end 
    set(handles.DrawROIbutton,'Enable','off');
    set(handles.PlanePopup,'Enable','off');
    set(handles.ZoomButton,'Enable','off');
    set(handles.UnzoomButton,'Enable','off');
    
    guidata(hObject,handles);
    updateImages(handles);
    else
        fprintf('An ROI has already been placed!\n');
    end 

    
% --- DELETE ROI BUTTON - CALLBACK
function DeleteROIbutton_Callback(hObject, eventdata, handles)
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
        if get(handles.PCRadio,'Value')
            set(handles.MAGRadio,'Enable','on');
        else 
            set(handles.PCRadio,'Enable','on');
        end 
    set(handles.ZoomButton,'Enable','on');
    set(handles.UnzoomButton,'Enable','on');

    
% --- LOAD ROI BUTTON - CALLBACK
function LoadROIbutton_Callback(hObject, eventdata, handles)
    global spanUD spanLR zoomed interpType requestInterp
    if ~requestInterp
        planeNum = get(handles.PlanePopup,'Value');
        if get(handles.PCRadio,'Value')
            v = handles.pcDatasets(planeNum).Data.v;
            if zoomed
                v = v(spanUD,spanLR,:);
            end 
            circle = handles.pcDatasets(planeNum).ROI;
            radius = circle.Radius; 
            center = round(circle.Center);

            [X,Y] = ndgrid(1:size(v,1),1:size(v,2));
            X = X-center(2); %shift coordinate grid
            Y = Y-center(1);
            roiMask = sqrt(X.^2+Y.^2)<=radius;

            %%% Create Uninterpolated Data
            for i=1:size(v,3)
                vTemp = v(:,:,i);
                roiDataRaw(:,i) = vTemp(roiMask);
                meanROI(i) = mean(vTemp(roiMask));
                medianROI(i) = median(vTemp(roiMask));
                maxROI(i) = max(vTemp(roiMask));
                minROI(i) = min(vTemp(roiMask));
                stdROI(i) = std(vTemp(roiMask));
                if isfield(handles.pcDatasets(planeNum).Info,'matrixx')
                    matrixx = handles.pcDatasets(planeNum).Info.matrixx;
                    fovx = handles.pcDatasets(planeNum).Info.fovx;
                    xres = fovx/matrixx;
                    timeres = handles.pcDatasets(planeNum).Info.timeres;
                else 
                    xres = handles.pcDatasets(planeNum).Info.PixelSpacing(1);
                    bpm = handles.pcDatasets(planeNum).Info.HeartRate;
                    frames = handles.pcDatasets(planeNum).Info.CardiacNumberOfImages;
                    rrInterval = (60*1000)/bpm;
                    timeres = rrInterval/frames;
                end 
                area = sum(roiMask(:))*(xres)^2;
                flowROI(i) = area.*meanROI(i);
            end 
            times = timeres.*(1:length(flowROI));
            dFlow = abs(derivative(flowROI));
            dMean = abs(derivative(meanROI));

            roiStatistics.radius = radius; roiStatistics.center = center;
            roiStatistics.roiMask = roiMask; roiStatistics.roiDataRaw = roiDataRaw;
            roiStatistics.times = times;
            roiStatistics.meanROI = meanROI; roiStatistics.dMean = dMean;
            roiStatistics.medianROI = medianROI;
            roiStatistics.maxROI = maxROI; roiStatistics.minROI = minROI;
            roiStatistics.stdROI = stdROI;
            roiStatistics.flowROI = flowROI; roiStatistics.dFlow = dFlow;

            if isstruct(handles.pcDatasets(planeNum).ROIdata)
                handles.pcDatasets(planeNum).ROIdata(end+1) = roiStatistics;
            else 
                handles.pcDatasets(planeNum).ROIdata = roiStatistics;
            end 

            %%% Create Makima Fit
            t = 1:length(flowROI);
            tq = 1:0.25:length(flowROI);
            interpTimes = timeres.*tq;
            meanROIfit = interp1(t,meanROI,tq,'makima');
            medianROIfit = interp1(t,medianROI,tq,'makima');
            maxROIfit = interp1(t,maxROI,tq,'makima');
            minROIfit = interp1(t,minROI,tq,'makima');
            stdROIfit = interp1(t,stdROI,tq,'makima');
            flowROIfit = interp1(t,flowROI,tq,'makima');
            dFlowfit = interp1(t,dFlow,tq,'makima');
            dMeanfit = interp1(t,dMean,tq,'makima');

            roiStatisticsMakima.times = interpTimes;
            roiStatisticsMakima.meanROI = meanROIfit; roiStatisticsMakima.dMean = dMeanfit;
            roiStatisticsMakima.medianROI = medianROIfit;
            roiStatisticsMakima.maxROI = maxROIfit; roiStatisticsMakima.minROI = minROIfit;
            roiStatisticsMakima.stdROI = stdROIfit;
            roiStatisticsMakima.flowROI = flowROIfit; roiStatisticsMakima.dFlow = dFlowfit;

            if isstruct(handles.pcDatasets(planeNum).ROIdataMakima)
                handles.pcDatasets(planeNum).ROIdataMakima(end+1) = roiStatisticsMakima;
            else 
                handles.pcDatasets(planeNum).ROIdataMakima = roiStatisticsMakima;
            end 

            %%% Create Piecewise Cubic Hermite Polynomial (PCHIP) Fit            
            meanROIfit = interp1(t,meanROI,tq,'pchip');
            medianROIfit = interp1(t,medianROI,tq,'pchip');
            maxROIfit = interp1(t,maxROI,tq,'pchip');
            minROIfit = interp1(t,minROI,tq,'pchip');
            stdROIfit = interp1(t,stdROI,tq,'pchip');
            flowROIfit = interp1(t,flowROI,tq,'pchip');
            dFlowfit = interp1(t,dFlow,tq,'pchip');
            dMeanfit = interp1(t,dMean,tq,'pchip');

            roiStatisticsPCHIP.times = interpTimes;
            roiStatisticsPCHIP.meanROI = meanROIfit; roiStatisticsPCHIP.dMean = dMeanfit;
            roiStatisticsPCHIP.medianROI = medianROIfit;
            roiStatisticsPCHIP.maxROI = maxROIfit; roiStatisticsPCHIP.minROI = minROIfit;
            roiStatisticsPCHIP.stdROI = stdROIfit;
            roiStatisticsPCHIP.flowROI = flowROIfit; roiStatisticsPCHIP.dFlow = dFlowfit;

            if isstruct(handles.pcDatasets(planeNum).ROIdataPCHIP)
                handles.pcDatasets(planeNum).ROIdataPCHIP(end+1) = roiStatisticsPCHIP;
            else 
                handles.pcDatasets(planeNum).ROIdataPCHIP = roiStatisticsPCHIP;
            end 

            %%% Create Cubic Spline Fit        
            meanROIfit = interp1(t,meanROI,tq,'spline');
            medianROIfit = interp1(t,medianROI,tq,'spline');
            maxROIfit = interp1(t,maxROI,tq,'spline');
            minROIfit = interp1(t,minROI,tq,'spline');
            stdROIfit = interp1(t,stdROI,tq,'spline');
            flowROIfit = interp1(t,flowROI,tq,'spline');
            dFlowfit = interp1(t,dFlow,tq,'spline');
            dMeanfit = interp1(t,dMean,tq,'spline');

            roiStatisticsSpline.times = interpTimes;
            roiStatisticsSpline.meanROI = meanROIfit; roiStatisticsSpline.dMean = dMeanfit;
            roiStatisticsSpline.medianROI = medianROIfit;
            roiStatisticsSpline.maxROI = maxROIfit; roiStatisticsSpline.minROI = minROIfit;
            roiStatisticsSpline.stdROI = stdROIfit;
            roiStatisticsSpline.flowROI = flowROIfit; roiStatisticsSpline.dFlow = dFlowfit;

            if isstruct(handles.pcDatasets(planeNum).ROIdataSpline)
                handles.pcDatasets(planeNum).ROIdataSpline(end+1) = roiStatisticsSpline;
            else 
                handles.pcDatasets(planeNum).ROIdataSpline = roiStatisticsSpline;
            end 

            clear dFlow dMean radius center roiMask area planeNum timeres v
            clear roiDataRaw meanROI medianROI maxROI minROI stdROI flowROI  
            clear dFlowfit dMeanfit circle matrixx t tq vTemp X Y xres i fovx
            clear meanROIfit medianROIfit maxROIfit minROIfit stdROIfit flowROIfit  
        end 
    end
    
    %%% Plot Curves
    handles.flow = organizeFlowInfo(handles);
    legendSet = {'Baseline',handles.flow.Name};
    count = length(legendSet)-1;
    times = handles.flow(1).Data.times;
    
    cla(handles.VelocityPlot,'reset')
    axes(handles.VelocityPlot);
    hold on; plot(times, zeros(1,length(times)) ,'Color','black','LineWidth',1.5);
    for i=1:count
        switch interpType
            case 'None'
                flow = handles.flow(i).Data.meanROI;
            case 'Makima'
                times = handles.flow(i).Makima.times;
                flow = handles.flow(i).Makima.meanROI;
            case 'PCHIP'
                times = handles.flow(i).PCHIP.times;
                flow = handles.flow(i).PCHIP.meanROI;           
            case 'Spline'
                times = handles.flow(i).Spline.times;
                flow = handles.flow(i).Spline.meanROI;    
            otherwise
        end
        
        if mean(flow)<0
            hold on; plot(times, -1*flow);
        else 
            hold on; plot(times, flow);
        end 
    end 
    legend(legendSet); hold off
    xlabel('Time (ms)'); ylabel('Mean Velocity in ROI (mm/s)');
    xlim([times(1),times(end)]);
    
    clear times interpTimes area count i 
    clear roiStatistics roiStatisticsMakima roiStatisticsPCHIP roiStatisticsSpline
    set(handles.InterpolatePopup,'Enable','on');
    
%     if numel(handles.flow)>1
%         allROIs = 1:numel(handles.flow);
%         handles.centerline.whichPlanes = [];
%         handles.centerline.whichIndices = {};
%         for i=1:numel(handles.flow)
%             ROIs2compare = allROIs ~= i;
%             ROIs2compare = nonzeros(ROIs2compare.*allROIs); 
%             for j=1:length(ROIs2compare)
%                 iterator = ROIs2compare(j);
%                 if iterator>i
%                     currName = handles.flow(i).Name;
%                     compareString = [currName ' --> ' handles.flow(iterator).Name];
%                     handles.centerline.whichPlanes(end+1) = compareString;
%                     handles.centerline.whichIndices(end+1) = [i,iterator];
%                 end 
%             end  
%         end 
%     end 

    requestInterp = 0;
    guidata(hObject,handles);
    DeleteROIbutton_Callback(hObject, eventdata, handles);
    updateImages(handles);
    axes(handles.PlanePlot);

    
% --- COMPLETE LOADING BUTTON - CALLBACK
function CompleteLoadingROI_Callback(hObject, eventdata, handles)
    set(handles.DrawROIbutton,'Enable','off');
    set(handles.DeleteROIbutton,'Enable','off');
    set(handles.LoadROIbutton,'Enable','off');

    set(handles.ShowPlanesRadio,'Enable','on');
    set(handles.ZoomAnatomicalButton,'Enable','on');
    set(handles.UnzoomAnatomicalButton,'Enable','on');
    set(handles.IdentifyPlanesButton,'Enable','on');


    
    
%%%%%%%%%%%% VELOCITY PLOT %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- VELOCITY PLOT - CREATE FUNCTION
function VelocityPlot_CreateFcn(hObject, eventdata, handles)
    

% --- INTERPOLATE POPUP - CALLBACK
function InterpolatePopup_Callback(hObject, eventdata, handles)
    global interpType requestInterp
    interp = get(handles.InterpolatePopup,'Value');
    requestInterp = 1;
    switch interp
        case 1
            interpType = 'None';
        case 2 
            interpType = 'Makima';
        case 3 
            interpType = 'PCHIP';
        case 4
            interpType = 'Spline';
        otherwise
    end 
    LoadROIbutton_Callback(hObject, eventdata, handles)

% --- INTERPOLATE POPUP - CREATE FUNCTION
function InterpolatePopup_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    
  
    
    
%%%%%%%%%%%% ANATOMICAL PLANE %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% --- ANATOMICAL PLOT - CREATE FUNCTION
function AnatomicalPlot_CreateFcn(hObject, eventdata, handles)


% --- ANATOMICAL SLIDER - CALLBACK
function SliderAnatomical_Callback(hObject, eventdata, handles)
    updateAnatImages(handles);

% --- ANATOMICAL SLIDER - CREATE FUNCTION
function SliderAnatomical_CreateFcn(hObject, eventdata, handles)
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end

 
% --- ANATOMICAL LISTBOX - CALLBACK
function AnatListbox_Callback(hObject, eventdata, handles)
    updateAnatImages(handles);

% --- ANATOMICAL LISTBOX - CREATE FUNCTION
function AnatListbox_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

    
% --- SHOW PLANE RADIO - CALLBACK
function ShowPlanesRadio_Callback(hObject, eventdata, handles) 
 
    
% --- ZOOM ANATOMICAL - CALLBACK
function ZoomAnatomicalButton_Callback(hObject, eventdata, handles)
    global zoomedAnat spanUDAnat spanLRAnat
    disp('Draw an ROI to zoom in');
    axes(handles.AnatomicalPlot)
    rect = drawrectangle;
    positions = round(rect.Position);
    spanUDAnat = spanUDAnat(1)+positions(2):(spanUDAnat(1)+positions(2)+positions(4));
    spanLRAnat = spanLRAnat(1)+positions(1):(spanLRAnat(1)+positions(1)+positions(3));
    zoomedAnat = 1;
    updateAnatImages(handles);
    
    
% --- UNZOOM ANATOMICAL - CALLBACK
function UnzoomAnatomicalButton_Callback(hObject, eventdata, handles)
    global zoomedAnat spanUDAnat spanLRAnat
    zoomedAnat = 0; spanUDAnat = 0; spanLRAnat = 0;
    updateAnatImages(handles);

    
% --- IDENTIFY PLANES - CALLBACK
function IdentifyPlanesButton_Callback(hObject, eventdata, handles)
    set(handles.DrawCenterlineButton,'Enable','on');
    set(handles.DeleteCenterlineButton,'Enable','on');
    set(handles.ComputePWVButton,'Enable','on');
    
    set(handles.ShowPlanesRadio,'Enable','off');
    set(handles.SliderAnatomical,'Enable','off');
    set(handles.AnatListbox,'Enable','off');
    set(handles.ZoomAnatomicalButton,'Enable','off');
    set(handles.UnzoomAnatomicalButton,'Enable','off');
    
    datasetNum = get(handles.AnatListbox,'Value');
    roiNames = {handles.flow.Name};
    axes(handles.AnatomicalPlot);   
    if isempty(handles.anatDatasets(datasetNum).Centerline)
        mydlg = warndlg(['Select the physical location of the ' num2str(length(roiNames)) ' chosen planes.']);
        waitfor(mydlg);
        for i=1:length(roiNames)
            handles.centerline(i).Name = roiNames{i};
            point = drawpoint;
            handles.centerline(i).Point = point;
        end 
    guidata(hObject,handles);
    updateImages(handles);
    else
        fprintf('A Centerline has already been placed!\n');
    end 

% --- DRAW CENTERLINE - CALLBACK
function DrawCenterlineButton_Callback(hObject, eventdata, handles)
    set(handles.DrawROIbutton,'Enable','off');

    datasetNum = get(handles.AnatListbox,'Value');
    axes(handles.AnatomicalPlot);   
    if isempty(handles.anatDatasets(datasetNum).Centerline)
        mydlg = warndlg('Press enter when the centerline is drawn');
        waitfor(mydlg);
        line = drawpolyline('Color','r','LineWidth',1);
        while true
            w = waitforbuttonpress; 
            switch w 
                case 1 % (keyboard press) 
                  key = get(gcf,'currentcharacter'); 
                      switch key
                          case 27 % 27 is the escape key
                              disp('User pressed the escape key. Deleting ROI.')
                              DeleteCenterlineButton_Callback(hObject, eventdata, handles)
                              break % break out of the while loop
                          case 13 % 13 is the return key 
                              disp('ROI selected');
                              line.InteractionsAllowed = 'none';
                              break
                          otherwise 
                              % Wait for a different command. 
                      end
           end
        end
        [splinePositions,lineLength] = interppolygon(line.Position,100);
        line.Position = splinePositions;
 
        if isfield(handles.anatDatasets(datasetNum).Info,'matrixx')
            matrixx = handles.anatDatasets(datasetNum).Info.matrixx;
            fovx = handles.anatDatasets(datasetNum).Info.fovx;
            xres = fovx/matrixx;
        else 
            xres = handles.anatDatasets(datasetNum).Info.PixelSpacing(1);
        end
    handles.centerline(1).Length = lineLength.*xres;
    handles.centerline(1).Centerline = line;
    
    guidata(hObject,handles);
    updateImages(handles);
    else
        fprintf('A Centerline has already been placed!\n');
    end 

% --- DELETE CENTERLINE - CALLBACK
function DeleteCenterlineButton_Callback(hObject, eventdata, handles)
    hold off
    guidata(hObject,handles);
    cla(handles.AnatomicalPlot,'reset');
    updateAnatImages(handles);
    
    set(handles.ShowPlanesRadio,'Enable','on');
    set(handles.DrawROIbutton,'Enable','on');
    set(handles.AnatListbox,'Enable','on');
    set(handles.ZoomAnatomicalButton,'Enable','on');
    set(handles.UnzoomAnatomicalButton,'Enable','on');
    

% --- COMPUTE PWV - CALLBACK
function ComputePWVButton_Callback(hObject, eventdata, handles)
    guidata(hObject,handles);
    updateImages(handles);
    set(handles.MAGRadio,'Enable','on');
    set(handles.PCRadio,'Enable','on');
    set(handles.DrawROIbutton,'Enable','on');
    set(handles.PlanePopup,'Enable','on');
    set(handles.ZoomButton,'Enable','on');
    set(handles.UnzoomButton,'Enable','on');
    set(handles.AnatListbox,'Enable','on');


    
    
%%%%%%%%%%%% PWV PLOT %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% --- PWV PLOT - CREATE FUNCTION
function TimeVsDistance_CreateFcn(hObject, eventdata, handles)


% --- Executes on button press in ExportAnalysisButton.
function ExportAnalysisButton_Callback(hObject, eventdata, handles)


% --- MESSAGE BAR - CALLBACK
function MessageBar_Callback(hObject, eventdata, handles)

% --- MESSAGE BAR - CREATE FUNCTION
function MessageBar_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end    
   

% --- TTP READOUT - CREATE FUNCTION
function ttpData_CreateFcn(hObject, eventdata, handles)


% --- TTU READOUT - CREATE FUNCTION
function ttuData_CreateFcn(hObject, eventdata, handles)


% --- TTF READOUT - CREATE FUNCTION
function ttfData_CreateFcn(hObject, eventdata, handles)


% --- XCORR READOUT - CREATE FUNCTION
function xcorrData_CreateFcn(hObject, eventdata, handles)


% --- AVERAGE READOUT - CREATE FUNCTION
function averageData_CreateFcn(hObject, eventdata, handles)




%%%%%%%%%%%% MY FUNCTIONS %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

% --- Update Images in PLANE PLOT
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
        cla(handles.PlanePlot,'reset')
        imshow(slice,[])
    end 
    
% --- Update Images in ANATOMICAL PLOT
function updateAnatImages(handles)
    global zoomedAnat spanUDAnat spanLRAnat
    datasetNum = get(handles.AnatListbox,'Value');
    images = handles.anatDatasets(datasetNum).Data;

    if ndims(images)<3
        steps = [1 1];
        set(handles.SliderAnatomical, 'SliderStep', steps);
        anatSlice = images;
    else
        dim3size = size(images,3);
        steps = [1/(dim3size-1) 10/(dim3size-1)];
        set(handles.SliderAnatomical, 'SliderStep', steps);
        sliceNum = 1+round( get(handles.SliderAnatomical,'Value').*(dim3size-1) );
        anatSlice = images(:,:,sliceNum);
    end 
    
    if zoomedAnat
        anatSlice = anatSlice(spanUDAnat,spanLRAnat);
    end

    axes(handles.AnatomicalPlot); %show image
    imshow(anatSlice,[]);

 
% --- Calculate derivative   
function fPrime = derivative(f)
    f(end+1) = f(1);
    fPrime = diff(f);
    
    
% --- TTP - time to peak calculation
function ttp = computeTTP(flow)
    numROIs = numel(flow);
    for i=1:numROIs
        if mean(flow(i).Data.flowROI)<0
            flow(i).Data.flowROI = -1*flow(i).Data.flowROI;
        end 
        [~,idx] = max(flow(i).Data.flowROI);
        ttp(i).maxFlowIdx = idx;     
    end 
    
    timeres = flow(i).Data.times(1);
    allROIs = 1:numROIs;
    for i=1:numROIs
        ROIs2compare = allROIs ~= i;
        ROIs2compare = nonzeros(ROIs2compare.*allROIs); 
        comparer = [NaN,NaN,NaN];
        for j=1:length(ROIs2compare)
            iterator = ROIs2compare(j);
            if iterator>1
                comparer(iterator) = timeres.*( ttp(iterator).maxFlowIdx - ttp(i).maxFlowIdx );
            end 
        end 
        ttp(i).deltaT = comparer;
    end 

    
% --- TTF - time to foot calculation
function ttf = computeTTF(flow)
    numROIs = numel(flow);
    for i=1:numROIs
        if mean(flow(i).Data.flowROI)<0
            flow(i).Data.flowROI = -1*flow(i).Data.flowROI;
        end 
        
    end 

    
% --- TTU - time to upstroke calculation    
function ttu = computeTTU(flow)
    numROIs = numel(flow);
    for i=1:numROIs
        if mean(flow(i).Data.flowROI)<0
            flow(i).Data.flowROI = -1*flow(i).Data.flowROI;
        end 
        
    end 

    
% --- Condense and organize flow data obtained from ROIs
function flow = organizeFlowInfo(handles)
    count = 1;
    for i=1:numel(handles.pcDatasets)
        if isstruct(handles.pcDatasets(i).ROIdata)
            if length(handles.pcDatasets(i).ROIdata)==1
                flow(count).Name = handles.pcDatasets(i).Names;
                flow(count).Data = handles.pcDatasets(i).ROIdata;
                flow(count).Info = handles.pcDatasets(i).Info;
                flow(count).Makima = handles.pcDatasets(i).ROIdataMakima;
                flow(count).PCHIP = handles.pcDatasets(i).ROIdataPCHIP;
                flow(count).Spline = handles.pcDatasets(i).ROIdataSpline;
                count = count+1;
            else 
                for j=1:length(handles.pcDatasets(i).ROIdata)
                    name = handles.pcDatasets(i).Names;
                    planeName = [name ' Plane ' num2str(j)];
                    flow(count).Name = planeName;
                    flow(count).Data = handles.pcDatasets(i).ROIdata(j);
                    flow(count).Info = handles.pcDatasets(i).Info;
                    flow(count).Makima = handles.pcDatasets(i).ROIdataMakima(j);
                    flow(count).PCHIP = handles.pcDatasets(i).ROIdataPCHIP(j);
                    flow(count).Spline = handles.pcDatasets(i).ROIdataSpline(j);
                    count = count+1;
                end 
            end 
        end 
    end 

% --- Turn PolyLine into SplineLine    
function [Y,lengthPolyLine] = interppolygon(X,N)
    if nargin < 2 || N < 2
        N = 2;
    end
    n_dim = size(X,2);
    delta_X = 0;
    
    for dim = 1:n_dim
        delta_X = delta_X + diff(X(:,dim)).^2 ;
    end
    
    lengthBetweenPoints = sqrt(delta_X);
    lengthLine = sum(lengthBetweenPoints);
    orig_metric = [0; cumsum(lengthBetweenPoints/lengthLine)];
    
    interp_metric = (0:(1/(N-1)):1)';
    Y = interp1(orig_metric,X,interp_metric,'makima');
    
    for dim = 1:n_dim
        delta_X = delta_X + diff(X(:,dim)).^2 ;
    end
    
    lengthPolyLine = sum(sqrt(delta_X));
