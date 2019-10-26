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

    % Last Modified by GUIDE v2.5 25-Oct-2019 17:48:53

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
    global  interpType showErrorBars startAnalyzing
    
    handles.anatDatasets = varargin{1};
    handles.pcDatasets = varargin{2};
    handles.magDatasets = varargin{3};
    handles.pcDatasets(1).ROI = [];
    handles.magDatasets(1).ROI = [];
    handles.anatDatasets(1).Centerline = [];
    handles.pcDatasets(1).ROIdata = [];
    handles.pcDatasets(1).ROIdataMakima = [];
    handles.pcDatasets(1).ROIdataGaussian = [];
    handles.pcDatasets(1).ROIdataSpline = [];
    zoomed = 0;
    zoomedAnat = 0;
    spanUD = [0,0];
    spanLR = [0,0];
    spanUDAnat = [0,0];
    spanLRAnat = [0,0];
    interpType = 'None';
    showErrorBars = 0;
    startAnalyzing = 0;
    
    set(handles.PlanePopup,'String',{handles.pcDatasets.Names});
    set(handles.DatasetPopup,'String',fieldnames(handles.pcDatasets(1).Data));
    set(handles.InterpolatePopup,'String',{'None','Makima','Gaussian','Spline'});
    set(handles.AnatListbox,'String',{handles.anatDatasets.Names});
    set(handles.InterpolatePopup,'Enable','off');
    if numel(handles.magDatasets)==0
        set(handles.MAGRadio,'Enable','off')
    end 
    
    set(handles.ShowPlanesRadio,'Enable','off');
    set(handles.DrawCenterlineButton,'Enable','off');
    set(handles.DeleteCenterlineButton,'Enable','off');
    set(handles.DeleteCenterlineButton,'Enable','off');
    set(handles.ComputePWVButton,'Enable','off');
    set(handles.ShowPlanesRadio,'Enable','off')
    set(handles.ttpRadio,'Value',1);
    set(handles.ttpointRadio,'Value',1);
    set(handles.ttuRadio,'Value',1);
    set(handles.ttfRadio,'Value',1);
    set(handles.xcorrRadio,'Value',1);
        
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
    clear planeNum
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
    axes(handles.PlanePlot)
    
    disp('Draw an ROI to zoom in');
    rect = drawrectangle;
    positions = round(rect.Position);
    spanUD = spanUD(1)+positions(2):(spanUD(1)+positions(2)+positions(4));
    spanLR = spanLR(1)+positions(1):(spanLR(1)+positions(1)+positions(3));
    zoomed = 1;
    
    clear rect positions
    updateImages(handles);
 
    
% --- UNZOOM BUTTON - CALLBACK
function UnzoomButton_Callback(hObject, eventdata, handles)
   global zoomed spanUD spanLR
    axes(handles.PlanePlot)
    
    zoomed = 0; spanUD = 0; spanLR = 0;
    updateImages(handles);

    
% --- DRAWROI BUTTON - CALLBACK
function DrawROIbutton_Callback(hObject, eventdata, handles)
    axes(handles.PlanePlot);
    
    if get(handles.PCRadio,'Value')
        set(handles.MAGRadio,'Enable','off');
    else 
        set(handles.PCRadio,'Enable','off');
    end 
    set(handles.DrawROIbutton,'Enable','off');
    set(handles.PlanePopup,'Enable','off');
    set(handles.ZoomButton,'Enable','off');
    set(handles.UnzoomButton,'Enable','off');
    
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
    
    clear planeNum mydlg w key
    guidata(hObject,handles);
    updateImages(handles);
    else
        fprintf('An ROI has already been placed!\n');
    end 

    
% --- DELETE ROI BUTTON - CALLBACK
function DeleteROIbutton_Callback(hObject, eventdata, handles)
   axes(handles.PlanePlot)
    hold off
    
    planeNum = get(handles.PlanePopup,'Value');
    if get(handles.PCRadio,'Value')
        handles.pcDatasets(planeNum).ROI = [];
    else 
        handles.magDatasets(planeNum).ROI = [];
    end 

    set(handles.DrawROIbutton,'Enable','on');
    set(handles.PlanePopup,'Enable','on');
        if get(handles.PCRadio,'Value')
            set(handles.MAGRadio,'Enable','on');
        else 
            set(handles.PCRadio,'Enable','on');
        end 
    set(handles.ZoomButton,'Enable','on');
    set(handles.UnzoomButton,'Enable','on');

    clear planeNum 
    guidata(hObject,handles);
    updateImages(handles);
    
% --- LOAD ROI BUTTON - CALLBACK
function LoadROIbutton_Callback(hObject, eventdata, handles)
    global spanUD spanLR zoomed
    
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

        %%% Create Linear Interpolated Data
        t = 1:size(v,3);
        tq = 1:0.25:size(v,3);
        
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
        
        for i=1:size(v,3)
            vTemp = v(:,:,i);
            roiDataRaw(:,i) = vTemp(roiMask);
            meanROI(i) = mean(vTemp(roiMask));
            medianROI(i) = median(vTemp(roiMask));
            maxROI(i) = max(vTemp(roiMask));
            minROI(i) = min(vTemp(roiMask));
            stdROI(i) = std(vTemp(roiMask));
            area = sum(roiMask(:))*(xres)^2;
            flowROI(i) = area.*meanROI(i);
        end 
        dFlow = abs(derivative(flowROI));
        dMean = abs(derivative(meanROI));
        
        times = timeres.*(1:length(flowROI));
        timesInterp = timeres.*tq;
        meanROIfit = interp1(times,meanROI,timesInterp,'linear');
        medianROIfit = interp1(times,medianROI,timesInterp,'linear');
        maxROIfit = interp1(times,maxROI,timesInterp,'linear');
        minROIfit = interp1(times,minROI,timesInterp,'linear');
        stdROIfit = interp1(times,stdROI,timesInterp,'linear');
        flowROIfit = interp1(times,flowROI,timesInterp,'linear');
        dFlowfit = interp1(times,dFlow,timesInterp,'linear');
        dMeanfit = interp1(times,dMean,timesInterp,'linear');
        

        roiStatistics.radius = radius; roiStatistics.center = center;
        roiStatistics.roiMask = roiMask; roiStatistics.roiDataRaw = roiDataRaw;
        roiStatistics.times = timesInterp;
        roiStatistics.meanROI = meanROIfit; roiStatistics.dMean = dMeanfit;
        roiStatistics.medianROI = medianROIfit;
        roiStatistics.maxROI = maxROIfit; roiStatistics.minROI = minROIfit;
        roiStatistics.stdROI = stdROIfit;
        roiStatistics.flowROI = flowROIfit; roiStatistics.dFlow = dFlowfit;

        if isstruct(handles.pcDatasets(planeNum).ROIdata)
            handles.pcDatasets(planeNum).ROIdata(end+1) = roiStatistics;
        else 
            handles.pcDatasets(planeNum).ROIdata = roiStatistics;
        end 

        %%% Create Makima Fit
        meanROIfit = interp1(times,meanROI,timesInterp,'makima');
        medianROIfit = interp1(times,medianROI,timesInterp,'makima');
        maxROIfit = interp1(times,maxROI,timesInterp,'makima');
        minROIfit = interp1(times,minROI,timesInterp,'makima');
        stdROIfit = interp1(times,stdROI,timesInterp,'makima');
        flowROIfit = interp1(times,flowROI,timesInterp,'makima');
        dFlowfit = interp1(times,dFlow,timesInterp,'makima');
        dMeanfit = interp1(times,dMean,timesInterp,'makima');

        roiStatisticsMakima.times = timesInterp;
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

        %%% Create Interpolated Curve with Gaussian Smoothing           
        meanROIfit = interp1(times,smoothdata(meanROI,'gaussian',5),timesInterp,'linear');
        medianROIfit = interp1(times,smoothdata(medianROI,'gaussian',5),timesInterp,'linear');
        maxROIfit = interp1(times,smoothdata(maxROI,'gaussian',5),timesInterp,'linear');
        minROIfit = interp1(times,smoothdata(minROI,'gaussian',5),timesInterp,'linear');
        stdROIfit = interp1(times,smoothdata(stdROI,'gaussian',5),timesInterp,'linear');
        flowROIfit = interp1(times,smoothdata(flowROI,'gaussian',5),timesInterp,'linear');
        dFlowfit = interp1(times,smoothdata(dFlow,'gaussian',5),timesInterp,'linear');
        dMeanfit = interp1(times,smoothdata(dMean,'gaussian',5),timesInterp,'linear');

        roiStatisticsGaussian.radius = radius; roiStatisticsGaussian.center = center;
        roiStatisticsGaussian.roiMask = roiMask; roiStatisticsGaussian.roiDataRaw = roiDataRaw;
        roiStatisticsGaussian.times = timesInterp;
        roiStatisticsGaussian.meanROI = meanROIfit; roiStatisticsGaussian.dMean = dMeanfit;
        roiStatisticsGaussian.medianROI = medianROIfit;
        roiStatisticsGaussian.maxROI = maxROIfit; roiStatisticsGaussian.minROI = minROIfit;
        roiStatisticsGaussian.stdROI = stdROIfit;
        roiStatisticsGaussian.flowROI = flowROIfit; roiStatisticsGaussian.dFlow = dFlowfit;

        if isstruct(handles.pcDatasets(planeNum).ROIdataGaussian)
            handles.pcDatasets(planeNum).ROIdataGaussian(end+1) = roiStatisticsGaussian;
        else 
            handles.pcDatasets(planeNum).ROIdataGaussian = roiStatisticsGaussian;
        end 

        %%% Create Cubic Spline Fit   
        meanROIfit = csaps([0 times times(end)+times(1)],[0 meanROI 0],0.0001,timesInterp);
        medianROIfit = csaps([0 times times(end)+times(1)],[0 medianROI 0],0.0001,timesInterp);
        maxROIfit = csaps([0 times times(end)+times(1)],[0 maxROI 0],0.0001,timesInterp);
        minROIfit = csaps([0 times times(end)+times(1)],[0 minROI 0],0.0001,timesInterp);
        stdROIfit = csaps([0 times times(end)+times(1)],[0 stdROI 0],0.0001,timesInterp);
        flowROIfit = csaps([0 times times(end)+times(1)],[0 flowROI 0],0.0001,timesInterp);
        dFlowfit = csaps([0 times times(end)+times(1)],[0 dFlow 0],0.0001,timesInterp);
        dMeanfit = csaps([0 times times(end)+times(1)],[0 dMean 0],0.0001,timesInterp);

        roiStatisticsSpline.times = timesInterp;
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
        clear roiStatistics roiStatisticsMakima roiStatisticsGaussian roiStatisticsSpline
        clear times timesInterp
    end 

    %%% Plot Curves
    handles.flow = organizeFlowInfo(handles);
    plotVelocity(handles);

    set(handles.InterpolatePopup,'Enable','on');
    guidata(hObject,handles);
    DeleteROIbutton_Callback(hObject, eventdata, handles);
    updateImages(handles);
    axes(handles.PlanePlot);
    
% --- COMPLETE LOADING BUTTON - CALLBACK
function CompleteLoadingROI_Callback(hObject, eventdata, handles)
    if isfield(handles,'flow') && numel(handles.flow)>1
        set(handles.DeleteROIbutton,'Enable','off');
        set(handles.LoadROIbutton,'Enable','off');

        set(handles.DrawCenterlineButton,'Enable','on');
        set(handles.ShowPlanesRadio,'Enable','on');
        set(handles.ZoomAnatomicalButton,'Enable','on');
        set(handles.UnzoomAnatomicalButton,'Enable','on');

        zLocs = EnterZslices(handles.flow);

        for i=1:numel(handles.flow)
            handles.flow(i).zLocs = zLocs(i);
        end 

        for i =1:numel(handles.anatDatasets)
            images = handles.anatDatasets(i).Data;
            dims(1) = size(images,1);
            dims(2) = size(images,2);
            dims(3) = size(images,3);
            handles.anatDatasets(i).dims = dims;

            originShift = [handles.anatDatasets(i).Info.ImagePositionPatient;1]; % origin is top left corner of image
            xres = handles.anatDatasets(i).Info.PixelSpacing(1);
            yres = handles.anatDatasets(i).Info.PixelSpacing(2);
            zres = handles.anatDatasets(i).Info.SliceThickness;

            % sometimes get extremely small values that should be 0, so round
            xVector = round( handles.anatDatasets(i).Info.ImageOrientationPatient(1:3) ,8); % what direction rows run w/r/to x
            yVector = round( handles.anatDatasets(i).Info.ImageOrientationPatient(4:6) ,8); % what direction the cols run w/r/to y
            zVector = [cross(xVector,yVector);0];
            xVector = [xVector;0];
            yVector = [yVector;0];
            rotationMatrix = [xres*xVector yres*yVector zres*zVector originShift];
            sliceDim = find(rotationMatrix(:,3));
            handles.anatDatasets(i).rotationMatrix = rotationMatrix;     
            handles.anatDatasets(i).sliceDim = sliceDim;


            topFrontLeft = rotationMatrix*[1;1;1;1];
            topFrontLeft(4) = []; % remove dummy dimension
            backBottomRight = rotationMatrix*[dims(1);dims(2);dims(3);1];
            backBottomRight(4) = [];
            spanZ(1) = topFrontLeft(3);
            spanZ(2) = backBottomRight(3);
            handles.anatDatasets(i).spanZ = spanZ;
        end   

        clear backBottomRight dims i images originShift rotationMatrix sliceDim
        clear spanZ zLocs topFrontLeft xres xVector yres yVector zres zVector
        guidata(hObject,handles);
    else
        set(handles.MessageBar,'String','You need at least 2 ROIs for PWV measurements.');
    end 
    
    
%%%%%%%%%%%% VELOCITY PLOT %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- VELOCITY PLOT - CREATE FUNCTION
function VelocityPlot_CreateFcn(hObject, eventdata, handles)
    

% --- INTERPOLATE POPUP - CALLBACK
function InterpolatePopup_Callback(hObject, eventdata, handles)
    global interpType
    
    interp = get(handles.InterpolatePopup,'Value');
    switch interp
        case 1
            interpType = 'None';
        case 2 
            interpType = 'Makima';
        case 3 
            interpType = 'Gaussian';
        case 4
            interpType = 'Spline';
        otherwise
    end 
    
    clear interp
    plotVelocity(handles)

% --- INTERPOLATE POPUP - CREATE FUNCTION
function InterpolatePopup_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    
% --- ERROR BAR RADIO - CALLBACK
function ErrorBarRadio_Callback(hObject, eventdata, handles)
    global showErrorBars
    if get(handles.ErrorBarRadio,'Value')
        showErrorBars = 1;
    else 
        showErrorBars = 0;
    end 

    plotVelocity(handles)
    
 


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
    ShowPlanesRadio_Callback(hObject, eventdata, handles)
    updateAnatImages(handles);

% --- ANATOMICAL LISTBOX - CREATE FUNCTION
function AnatListbox_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

    
% --- SHOW PLANE RADIO - CALLBACK
function ShowPlanesRadio_Callback(hObject, eventdata, handles) 
    if get(handles.ShowPlanesRadio,'Value')
        datasetNum = get(handles.AnatListbox,'Value');
        yres = handles.anatDatasets(datasetNum).Info.PixelSpacing(2);
        anatRotation = handles.anatDatasets(datasetNum).rotationMatrix;
        colsRunningDir = sign(nonzeros(anatRotation(:,3)));
        spanZ = handles.anatDatasets(datasetNum).spanZ;
        for i=1:numel(handles.flow)
            planeLinePhysical = colsRunningDir*(handles.flow(i).zLocs-spanZ(1));
            planeLineRow = round(planeLinePhysical/yres);
            handles.flow(i).planeLineRow = planeLineRow;
            guidata(hObject,handles);
        end 
    end 
    
    clear i datasetNum yres anatRotation colsRunningDir spanZ 
    clear planeLinePhysical planeLineRow
    updateAnatImages(handles);
    guidata(hObject,handles);
    
% --- ZOOM ANATOMICAL - CALLBACK
function ZoomAnatomicalButton_Callback(hObject, eventdata, handles)
    global zoomedAnat spanUDAnat spanLRAnat
    axes(handles.AnatomicalPlot)
    
    disp('Draw an ROI to zoom in');
    rect = drawrectangle;
    positions = round(rect.Position);
    spanUDAnat = spanUDAnat(1)+positions(2):(spanUDAnat(1)+positions(2)+positions(4));
    spanLRAnat = spanLRAnat(1)+positions(1):(spanLRAnat(1)+positions(1)+positions(3));
    zoomedAnat = 1;
    
    clear rect positions 
    updateAnatImages(handles);
    
    
% --- UNZOOM ANATOMICAL - CALLBACK
function UnzoomAnatomicalButton_Callback(hObject, eventdata, handles)
    global zoomedAnat spanUDAnat spanLRAnat
    axes(handles.AnatomicalPlot)
    
    zoomedAnat = 0; spanUDAnat = 0; spanLRAnat = 0;
    updateAnatImages(handles);


% --- DRAW CENTERLINE - CALLBACK
function DrawCenterlineButton_Callback(hObject, eventdata, handles)
    axes(handles.AnatomicalPlot); 
    
    set(handles.ShowPlanesRadio,'Value',1);
    set(handles.ShowPlanesRadio,'Enable','off');
    set(handles.DrawROIbutton,'Enable','off');
    set(handles.DeleteCenterlineButton,'Enable','on');
    set(handles.ComputePWVButton,'Enable','on');
    set(handles.SliderAnatomical,'Enable','off');
    set(handles.AnatListbox,'Enable','off');
    set(handles.ZoomAnatomicalButton,'Enable','off');
    set(handles.UnzoomAnatomicalButton,'Enable','off');
    set(handles.ShowPlanesRadio,'Value',1);
    ShowPlanesRadio_Callback(hObject, eventdata, handles);
    
    datasetNum = get(handles.AnatListbox,'Value');
    yres = handles.anatDatasets(datasetNum).Info.PixelSpacing(2);
    anatRotation = handles.anatDatasets(datasetNum).rotationMatrix;
    colsRunningDir = sign(nonzeros(anatRotation(:,3)));
    spanZ = handles.anatDatasets(datasetNum).spanZ;
    for i=1:numel(handles.flow)
        planeLinePhysical = colsRunningDir*(handles.flow(i).zLocs-spanZ(1));
        planeLineRow = round(planeLinePhysical/yres);
        handles.flow(i).planeLineRow = planeLineRow;
        guidata(hObject,handles);
    end 

    datasetNum = get(handles.AnatListbox,'Value');
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
        
    splinePositions = interppolygon(line.Position,100);
    line.Position = splinePositions;
    
    if isfield(handles.anatDatasets(datasetNum).Info,'matrixx')
        matrixx = handles.anatDatasets(datasetNum).Info.matrixx;
        fovx = handles.anatDatasets(datasetNum).Info.fovx;
        xres = fovx/matrixx;
        matrixy = handles.anatDatasets(datasetNum).Info.matrixy;
        fovy = handles.anatDatasets(datasetNum).Info.fovy;
        yres = fovy/matrixy;
    else 
        xres = handles.anatDatasets(datasetNum).Info.PixelSpacing(1);
        yres = handles.anatDatasets(datasetNum).Info.PixelSpacing(2);
    end

    distances = zeros(1,length(splinePositions)-1);
    for i=1:length(splinePositions)-1
        distances(i) = sqrt( xres.*(splinePositions(i,1)-splinePositions(i+1,1)).^2 + yres.*(splinePositions(i,2)-splinePositions(i+1,2)).^2 );
    end 
    distances(end+1)=0;
    handles.anatDatasets(datasetNum).Distances = distances;
    handles.anatDatasets(datasetNum).Centerline = line;
    
    clear datasetNum yres anatRotation colsRunningDir spanZ i w mydlg 
    clear key planeLinePhysical planeLineRow matrixx fovx xres matrixy 
    clear fovy yres distances splinePositions
    guidata(hObject,handles);
    updateImages(handles);
    else
        fprintf('A Centerline has already been placed!\n');
    end 

% --- DELETE CENTERLINE - CALLBACK
function DeleteCenterlineButton_Callback(hObject, eventdata, handles)
    cla(handles.AnatomicalPlot,'reset');
    
    datasetNum = get(handles.AnatListbox,'Value');
    handles.anatDatasets(datasetNum).Centerline = [];
    handles.anatDatasets(datasetNum).Length = [];

    set(handles.ShowPlanesRadio,'Enable','on');
    set(handles.DrawROIbutton,'Enable','on');
    set(handles.AnatListbox,'Enable','on');
    set(handles.ZoomAnatomicalButton,'Enable','on');
    set(handles.UnzoomAnatomicalButton,'Enable','on');
    set(handles.SliderAnatomical,'Enable','on');
    set(handles.ShowPlanesRadio,'Value',0);
    
    
    clear datasetNum
    hold off
    guidata(hObject,handles);
    updateAnatImages(handles);
   
    
% --- COMPUTE PWV - CALLBACK
function ComputePWVButton_Callback(hObject, eventdata, handles)
    global startAnalyzing
    
    startAnalyzing = 1;
    set(handles.DrawROIbutton,'Enable','off');
    set(handles.DeleteCenterlineButton,'Enable','off');
    set(handles.ShowPlanesRadio,'Enable','on');
    set(handles.AnatListbox,'Enable','on');
    set(handles.ZoomAnatomicalButton,'Enable','on');
    set(handles.UnzoomAnatomicalButton,'Enable','on');
    set(handles.SliderAnatomical,'Enable','on');

    distances = [handles.anatDatasets.Distances];
    centerline = [handles.anatDatasets.Centerline];
    planeRows = [handles.flow.planeLineRow];
    y = centerline.Position(:,2);
    
    difference = ones(length(y),1);
    for i=1:length(planeRows)
        difference = difference.*(y-planeRows(i));
    end 
    minima = islocalmin(abs(difference));
    
    if y(1)-planeRows(1)<0
        minima(1) = 1;
    end 
    
    if y(end)-planeRows(end)<0
        minima(end) = 1;
    end 
    roiIdxAlongLine = nonzeros(minima.*(1:length(y))');
    
    numROIs = numel(handles.flow);  
    allROIs = 1:numROIs;
    for i=1:numROIs
        ROIs2compare = allROIs ~= i;
        ROIs2compare = nonzeros(ROIs2compare.*allROIs); 
        distances2ROIs = [NaN,NaN,NaN];
        for j=1:length(ROIs2compare)
            iterator = ROIs2compare(j);
            if iterator>i
                idx1 = roiIdxAlongLine(i);
                idx2 = roiIdxAlongLine(iterator);
                dist = sum(distances(idx1:idx2));
                distances2ROIs(iterator) = dist;
            end 
        end 
        handles.flow(i).Distance = distances2ROIs;
    end 
    
    guidata(hObject,handles);
    TTs = computeTTs(handles.flow);
  
    distance = []; TTPeak = []; TTPoint = []; TTFoot = []; TTUpstroke = []; Xcorr = [];
    for i = 1:numel(handles.flow)
        mask = ~isnan(handles.flow(i).Distance);
        distance = [distance handles.flow(i).Distance(mask)];
        TTPeak = [TTPeak TTs(i).TTPeak(mask)];
        TTPoint = [TTPoint TTs(i).TTPoint(mask)];
        TTFoot = [TTFoot TTs(i).TTFoot(mask)];
        TTUpstroke = [TTUpstroke TTs(i).TTUpstroke(mask)];
        Xcorr = [Xcorr TTs(i).Xcorr(mask)];
    end
    
    if numel(distance)==1
        set(handles.distanceHeader1,'String',' Plane 1 --> 2');
        set(handles.distance1,'String',[num2str(round(distance(1),1)) ' mm']); 
    elseif numel(distance)==3
        set(handles.distanceHeader1,'String',' Plane 1 --> 2');
        set(handles.distance1,'String',[num2str(round(distance(1),1)) ' mm']); 
        set(handles.distanceHeader2,'String',' Plane 1 --> 3');
        set(handles.distance2,'String',[num2str(round(distance(2),1)) ' mm']); 
        set(handles.distanceHeader3,'String',' Plane 2 --> 3');
        set(handles.distance3,'String',[num2str(round(distance(3),1)) ' mm']);       
    else
        set(handles.distanceHeader1,'String',' Plane 1 --> 2');
        set(handles.distance1,'String',[num2str(round(distance(1),1)) ' mm']); 
        set(handles.distanceHeader2,'String',' Plane 1 --> 3');
        set(handles.distance2,'String',[num2str(round(distance(2),1)) ' mm']); 
        set(handles.distanceHeader3,'String',' Plane 1 --> 4');
        set(handles.distance3,'String',[num2str(round(distance(3),1)) ' mm']); 
        set(handles.distanceHeader4,'String',' Plane 2 --> 3');
        set(handles.distance4,'String',[num2str(round(distance(4),1)) ' mm']);
        set(handles.distanceHeader5,'String',' Plane 2 --> 4');
        set(handles.distance5,'String',[num2str(round(distance(5),1)) ' mm']); 
        set(handles.distanceHeader6,'String',' Plane 3 --> 4');
        set(handles.distance6,'String',[num2str(round(distance(6),1)) ' mm']);         
    end 
    
    cla(handles.TimeVsDistance,'reset');
    axes(handles.TimeVsDistance); hold on; 
    xlabel('Distance (mm)'); ylabel('Time Shift (ms)');
    sz = 30;
    
    numMethods = get(handles.ttpRadio,'Value')+get(handles.ttpointRadio,'Value')...
    +get(handles.ttfRadio,'Value')+get(handles.ttuRadio,'Value')+get(handles.xcorrRadio,'Value');
    numCompares = numel(distance);
    
    legendSet = {};
    average = zeros(1,numCompares);
    if get(handles.ttpRadio,'Value')
        scatter(distance,TTPeak,sz,'filled','MarkerFaceColor',[0 0.4470 0.7410]);
        for i=1:numCompares
            average(i) = average(i)+TTPeak(i);
        end 
        legendSet{end+1} = 'TTPeak';
    end 
    if get(handles.ttpointRadio,'Value')
        scatter(distance,TTPoint,sz,'filled','MarkerFaceColor',[0.8500 0.3250 0.0980]);
        for i=1:numCompares
            average(i) = average(i)+TTPoint(i);
        end 
        legendSet{end+1} = 'TTPoint';
    end 
    if get(handles.ttfRadio,'Value')
        scatter(distance,TTFoot,sz,'filled','MarkerFaceColor',[0.4940 0.1840 0.5560]);
        for i=1:numCompares
            average(i) = average(i)+TTFoot(i);
        end 
        legendSet{end+1} = 'TTFoot';
    end 
    if get(handles.ttuRadio,'Value')
        scatter(distance,TTUpstroke,sz,'filled','MarkerFaceColor',[0.4660 0.6740 0.1880]);
        for i=1:numCompares
            average(i) = average(i)+TTUpstroke(i);
        end 
        legendSet{end+1} = 'TTUpstroke';
    end 
    if get(handles.xcorrRadio,'Value')
        scatter(distance,Xcorr,sz,'filled','MarkerFaceColor',[0.6350 0.0780 0.1840]);
        for i=1:numCompares
            average(i) = average(i)+Xcorr(i);
        end 
        legendSet{end+1} = 'XCorr';
    end 

    average = average/numMethods;
    scatter(distance,average,40,'black');
    legendSet{end+1} = 'AVERAGE';
    

    d = min(distance):max(distance);
    linePeakFit = polyfit(distance,TTPeak,1);
    linePeak = linePeakFit(1)*d + linePeakFit(2);
    linePointFit = polyfit(distance,TTPoint,1);
    linePoint= linePointFit(1)*d + linePointFit(2);
    lineFootFit = polyfit(distance,TTFoot,1);
    lineFoot= lineFootFit(1)*d + lineFootFit(2);
    lineUpstrokeFit = polyfit(distance,TTUpstroke,1);
    lineUpstroke = lineUpstrokeFit(1)*d + lineUpstrokeFit(2);
    lineXcorrFit = polyfit(distance,Xcorr,1);
    lineXcorr = lineXcorrFit(1)*d + lineXcorrFit(2);
    lineAverageFit = polyfit(distance,average,1);
    lineAverage= lineAverageFit(1)*d + lineAverageFit(2);
    PWVpeak = 1/linePeakFit(1);
    PWVpoint = 1/linePointFit(1);
    PWVfoot = 1/lineFootFit(1);
    PWVupstroke = 1/lineUpstrokeFit(1);
    PWVxcorr = 1/lineXcorrFit(1);
    PWVaverage = 1/lineAverageFit(1);
    
    hold on; 
    if get(handles.ttpRadio,'Value')
        plot(d,linePeak,':','LineWidth',0.2,'MarkerFaceColor',[0 0.4470 0.7410]);
        set(handles.ttpData,'String',[num2str(round(PWVpeak,2)) ' m/s']);
    end 
    if get(handles.ttpointRadio,'Value')
        plot(d,linePoint,':','LineWidth',0.2,'MarkerFaceColor',[0.8500 0.3250 0.0980]);
        set(handles.ttpointData,'String',[num2str(round(PWVpoint,2)) ' m/s']);
    end 
    if get(handles.ttfRadio,'Value')
        plot(d,lineFoot,':','LineWidth',0.2,'MarkerFaceColor',[0.4940 0.1840 0.5560]);
        set(handles.ttfData,'String',[num2str(round(PWVfoot,2)) ' m/s']);
    end 
    if get(handles.ttuRadio,'Value')
        plot(d,lineUpstroke,':','LineWidth',0.2,'MarkerFaceColor',[0.4660 0.6740 0.1880]);
        set(handles.ttuData,'String',[num2str(round(PWVupstroke,2)) ' m/s']);
    end 
    if get(handles.xcorrRadio,'Value')
        plot(d,lineXcorr,':','LineWidth',0.2,'MarkerFaceColor',[0.6350 0.0780 0.1840]);
        set(handles.xcorrData,'String',[num2str(round(PWVxcorr,2)) ' m/s']);
    end 
    plot(d,lineAverage,'-k','LineWidth',0.2);
    legend(legendSet,'Location','southeast');
    hold off;
    
    
    if numMethods>0
        set(handles.averageData,'String',[num2str(round(PWVaverage,2)) ' m/s']);
    else 
        set(handles.averageData,'String','0 m/s');
    end
    
    clear allROIs average centerline d difference dist distance distances distances2ROIs
    clear i idx idx2 iterator j legendSet lineAverage lineAverageFit idx1
    clear lineFoot lineFootFit linePeak linePeakFit linePoint linePointFit
    clear lineUpstroke lineUpstrokeFit lineXcorr lineXcorrFit mask
    clear minima numCompares numMethods numROIs planeRows PWVaverage PWVfoot 
    clear PWVpeak PWVpoint PWVupstroke PWVxcorr sz startAnalyzing Xcorr
    clear TTFoot TTPeak TTPoint TTs TTUpstroke y ROIs2compare roiIdxAlongLine
    
% --- COMPLETE RESET BUTTON - CALLBACK
function resetButton_Callback(hObject, eventdata, handles)
    mydlg = warndlg('WARNING: Pressing ENTER will completely reset the ROI data. Press ESC to cancel');
    waitfor(mydlg);
    while true
        w = waitforbuttonpress; 
        switch w 
            case 1 % (keyboard press) 
              key = get(gcf,'currentcharacter'); 
                  switch key
                      case 27 % 27 is the escape key
                          set(handles.MessageBar,'String','Complete reset cancelled.');
                          break % break out of the while loop
                      case 13 % 13 is the return key 
                          set(handles.MessageBar,'String','Complete reset initiated.');
                          break
                      otherwise 
                          % Wait for a different command. 
                  end
       end
    end

    global zoomed zoomedAnat spanUD spanLR spanUDAnat spanLRAnat
    global  interpType showErrorBars startAnalyzing
    
    handles.pcDatasets(1).ROI = [];
    handles.magDatasets(1).ROI = [];
    handles.anatDatasets(1).Centerline = [];
    handles.pcDatasets(1).ROIdata = [];
    handles.pcDatasets(1).ROIdataMakima = [];
    handles.pcDatasets(1).ROIdataGaussian = [];
    handles.pcDatasets(1).ROIdataSpline = [];
    handles.flow = [];
    zoomed = 0;
    zoomedAnat = 0;
    spanUD = [0,0];
    spanLR = [0,0];
    spanUDAnat = [0,0];
    spanLRAnat = [0,0];
    interpType = 'None';
    showErrorBars = 0;
    startAnalyzing = 0;
    
    set(handles.PlanePopup,'String',{handles.pcDatasets.Names});
    set(handles.DatasetPopup,'String',fieldnames(handles.pcDatasets(1).Data));
    set(handles.InterpolatePopup,'String',{'None','Makima','Gaussian','Spline'});
    set(handles.AnatListbox,'String',{handles.anatDatasets.Names});
    set(handles.InterpolatePopup,'Enable','off');
    if numel(handles.magDatasets)==0
        set(handles.MAGRadio,'Enable','off')
    end 
    
    set(handles.ShowPlanesRadio,'Enable','off');
    set(handles.DrawCenterlineButton,'Enable','off');
    set(handles.DeleteCenterlineButton,'Enable','off');
    set(handles.DeleteCenterlineButton,'Enable','off');
    set(handles.ComputePWVButton,'Enable','off');
    set(handles.ShowPlanesRadio,'Enable','off')
    set(handles.ttpRadio,'Value',1);
    set(handles.ttpointRadio,'Value',1);
    set(handles.ttuRadio,'Value',1);
    set(handles.ttfRadio,'Value',1);
    set(handles.xcorrRadio,'Value',1);
        
    
    clear mydlg w key 
    guidata(hObject, handles);
    updateImages(handles);
    updateAnatImages(handles);



   
%%%%%%%%%%%% PWV PLOT %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% --- PWV PLOT - CREATE FUNCTION
function TimeVsDistance_CreateFcn(hObject, eventdata, handles)


% --- EXPORT ANALYSIS - CALLBACK
function ExportAnalysisButton_Callback(hObject, eventdata, handles)


% --- MESSAGE BAR - CALLBACK
function MessageBar_Callback(hObject, eventdata, handles)

% --- MESSAGE BAR - CREATE FUNCTION
function MessageBar_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end    
   

% --- TTPeak READOUT - CREATE FUNCTION
function ttpData_CreateFcn(hObject, eventdata, handles)
% --- TTPeak RADIO - CALLBACK
function ttpRadio_Callback(hObject, eventdata, handles)
    global startAnalyzing
    
    if ~get(handles.ttpRadio,'Value')
        set(handles.ttpData,'String',' ');
    end 
    if startAnalyzing
        ComputePWVButton_Callback(hObject, eventdata, handles);
    end 

% --- TTPoint READOUT - CREATE FUNCTION
function ttpointData_CreateFcn(hObject, eventdata, handles)
% --- TTPoint RADIO - CALLBACK
function ttpointRadio_Callback(hObject, eventdata, handles)
    global startAnalyzing
    
    if ~get(handles.ttpointRadio,'Value')
        set(handles.ttpointData,'String',' ');
    end 
    if startAnalyzing
        ComputePWVButton_Callback(hObject, eventdata, handles);
    end 

% --- TTUpstroke READOUT - CREATE FUNCTION
function ttuData_CreateFcn(hObject, eventdata, handles)
% --- TTUpstroke RADIO - CALLBACK
function ttuRadio_Callback(hObject, eventdata, handles)
    global startAnalyzing
    
    if ~get(handles.ttuRadio,'Value')
        set(handles.ttuData,'String',' ');
    end 
    if startAnalyzing
        ComputePWVButton_Callback(hObject, eventdata, handles);
    end 

% --- TTFoot READOUT - CREATE FUNCTION
function ttfData_CreateFcn(hObject, eventdata, handles)
% --- TTFoot RADIO - CALLBACK
function ttfRadio_Callback(hObject, eventdata, handles)
    global startAnalyzing
    
    if ~get(handles.ttfRadio,'Value')
        set(handles.ttfData,'String',' ');
    end 
    if startAnalyzing
        ComputePWVButton_Callback(hObject, eventdata, handles);
    end 

% --- Xcorr READOUT - CREATE FUNCTION
function xcorrData_CreateFcn(hObject, eventdata, handles)
% --- Xcorr RADIO - CALLBACK
function xcorrRadio_Callback(hObject, eventdata, handles)
    global startAnalyzing
    
    if ~get(handles.xcorrRadio,'Value')
        set(handles.xcorrData,'String',' ');
    end 
    if startAnalyzing
        ComputePWVButton_Callback(hObject, eventdata, handles);
    end 

% --- AVERAGE READOUT - CREATE FUNCTION
function averageData_CreateFcn(hObject, eventdata, handles)


% --- DISTANCE 1 PlANES - CREATE FUNCTION
function distanceHeader1_CreateFcn(hObject, eventdata, handles)
% --- DISTANCE 1 - CREATE FUNCTION
function distance1_CreateFcn(hObject, eventdata, handles)

% --- DISTANCE 2 PlANES - CREATE FUNCTION
function distanceHeader2_CreateFcn(hObject, eventdata, handles)
% --- DISTANCE 2 - CREATE FUNCTION
function distance2_CreateFcn(hObject, eventdata, handles)

% --- DISTANCE 3 PlANES - CREATE FUNCTION
function distanceHeader3_CreateFcn(hObject, eventdata, handles)
% --- DISTANCE 3 - CREATE FUNCTION
function distance3_CreateFcn(hObject, eventdata, handles)

% --- DISTANCE 4 PlANES - CREATE FUNCTION
function distanceHeader4_CreateFcn(hObject, eventdata, handles)
% --- DISTANCE 4 - CREATE FUNCTION
function distance4_CreateFcn(hObject, eventdata, handles)

% --- DISTANCE 5 PlANES - CREATE FUNCTION
function distanceHeader5_CreateFcn(hObject, eventdata, handles)
% --- DISTANCE 5 - CREATE FUNCTION
function distance5_CreateFcn(hObject, eventdata, handles)

% --- DISTANCE 6 PlANES - CREATE FUNCTION
function distance6_CreateFcn(hObject, eventdata, handles)
% --- DISTANCE 6 - CREATE FUNCTION
function distanceHeader6_CreateFcn(hObject, eventdata, handles)



%%%%%%%%%%%% MY FUNCTIONS %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

% --- Update Images in PLANE PLOT
function updateImages(handles)
    global spanUD spanLR zoomed 
    axes(handles.PlanePlot);
    
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
        steps = [1 500];
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

    if ~isempty(handles.pcDatasets(planeNum).ROI)
        hold on
        imshow(slice,[]);
    else 
        cla(handles.PlanePlot,'reset')
        imshow(slice,[])
    end 
    clear dataset datasetNum images imageSet isPC planeNum slice steps 
    
% --- Update Images in ANATOMICAL PLOT
function updateAnatImages(handles)
    global zoomedAnat spanUDAnat spanLRAnat 
    axes(handles.AnatomicalPlot); %show image
    
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
    
    if get(handles.ShowPlanesRadio,'Value') && handles.anatDatasets(datasetNum).sliceDim~=3
        imshow(anatSlice,[]);
        for i=1:numel(handles.flow)
            pcRow = handles.flow(i).planeLineRow;
            anatSlice(pcRow,:) = max(max(anatSlice))+100;
            hold on; plot(1:size(anatSlice,1),pcRow.*ones(1,size(anatSlice,2)),'y');
        end 
    else 
       imshow(anatSlice,[]); 
    end 
    clear anatSlice datasetNum dim3size i images pcRow sliceNum spanLRAnat spanUDAnat steps

    
% --- "Time to" calculations (TTPeak, TTPoint, TTUpstroke, TTFoot, Xcorr)
function TTs = computeTTs(flow)
    global startAnalyzing interpType
    
    numROIs = numel(flow);
    for i=1:numROIs
        if startAnalyzing
            switch interpType
                case 'None'
                    flowTemp = flow(i).Data.meanROI;
                case 'Makima'
                    flowTemp = flow(i).Makima.meanROI;
                case 'Gaussian'
                    flowTemp = flow(i).Gaussian.meanROI;  
                case 'Spline'
                    flowTemp = flow(i).Spline.meanROI; 
                otherwise
            end 
        else
            flowTemp = flow(i).Data.meanROI;
        end 
        
        if mean(flowTemp)<0
            flowTemp = -1*flowTemp;
        else 
            flowTemp = flowTemp;
        end 
        flowTemp = normalize(flowTemp,'range');
        
        times = flow(i).Data.times;
        timeres = times(2)-times(1);
        [maxPeakVel,maxPeakVelIdx] = max(flowTemp);
        upstroke = flowTemp(1:maxPeakVelIdx);
        [~,EightyPointIdx] = min(abs(upstroke-0.8*maxPeakVel));
        [~,FiftyPointIdx] = min(abs(upstroke-0.5*maxPeakVel));
        [~,TwentyPointIdx] = min(abs(upstroke-0.2*maxPeakVel));
        curvePoints(i).maxPeakVelIdx = maxPeakVelIdx;   
        curvePoints(i).maxPeakVel = maxPeakVel;        
        curvePoints(i).EightyPointIdx = EightyPointIdx;
        curvePoints(i).EightyPoint = flowTemp(EightyPointIdx);
        curvePoints(i).FiftyPointIdx = FiftyPointIdx;
        curvePoints(i).FiftyPoint = flowTemp(FiftyPointIdx);
        curvePoints(i).TwentyPointIdx = TwentyPointIdx;
        curvePoints(i).TwentyPoint = flowTemp(TwentyPointIdx);
        flows(i,:) = flowTemp;
    end 
    
    % TTPeak - time to peak calculation
    allROIs = 1:numROIs;
    for i=1:numROIs
        ROIs2compare = allROIs ~= i;
        ROIs2compare = nonzeros(ROIs2compare.*allROIs); 
        TTPeak = [NaN,NaN,NaN];
        for j=1:length(ROIs2compare)
            iterator = ROIs2compare(j);
            if iterator>i
                TTPeak(iterator) = timeres.*( curvePoints(iterator).maxPeakVelIdx - curvePoints(i).maxPeakVelIdx );
            end 
        end 
        TTs(i).TTPeak = TTPeak;
    end      
    
    % TTPoint - time to point calculation
    for i=1:numROIs
        ROIs2compare = allROIs ~= i;
        ROIs2compare = nonzeros(ROIs2compare.*allROIs); 
        TTPoint = [NaN,NaN,NaN];
        for j=1:length(ROIs2compare)
            iterator = ROIs2compare(j);
            if iterator>i
                TTPoint(iterator) = timeres.*( curvePoints(iterator).FiftyPointIdx - curvePoints(i).FiftyPointIdx );
            end 
        end 
        TTs(i).TTPoint = TTPoint;
    end    
    
    % TTUpstroke - time to upstroke calculation
    for i=1:numROIs
        ROIs2compare = allROIs ~= i;
        ROIs2compare = nonzeros(ROIs2compare.*allROIs); 
        TTUpstroke = [NaN,NaN,NaN];  
        
        for j=1:length(ROIs2compare)
            iterator = ROIs2compare(j);
            if iterator>i
                [sigmoid1,upslope1,~,tUpstroke1,~,timeresSigmoid1] = sigFit(flows(i,:));
                [sigmoid2,upslope2,~,tUpstroke2,~,timeresSigmoid2] = sigFit(flows(iterator,:));
                TTUpstroke(iterator) = timeres.*(timeresSigmoid1*tUpstroke2 - timeresSigmoid2*tUpstroke1);
            end 
        end 
        TTs(i).TTUpstroke = TTUpstroke;
    end   
    
    % TTFoot - time to foot calculation
    for i=1:numROIs
        ROIs2compare = allROIs ~= i;
        ROIs2compare = nonzeros(ROIs2compare.*allROIs); 
        TTFoot = [NaN,NaN,NaN];
        for j=1:length(ROIs2compare)
            iterator = ROIs2compare(j);
            if iterator>i
                % m = (y2-y1)/(x2-x1); y1=80%max flow and y2=20%max flow
                % x1 and x2 are the times (indices) where these values occur
                m1 = (curvePoints(i).EightyPoint-curvePoints(i).TwentyPoint)/(curvePoints(i).EightyPointIdx-curvePoints(i).TwentyPointIdx);
                % t1 is the x-intercept of this line. 
                % This can be solved analytically; b = x1-y1/m
                t1 = curvePoints(i).TwentyPointIdx - (curvePoints(i).TwentyPoint/m1);
                % Do the same for the consequent flow curve
                m2 = (curvePoints(iterator).EightyPoint-curvePoints(iterator).TwentyPoint)/(curvePoints(iterator).EightyPointIdx-curvePoints(iterator).TwentyPointIdx);
                t2 = curvePoints(iterator).TwentyPointIdx - (curvePoints(iterator).TwentyPoint/m2);
                TTFoot(iterator) = timeres.*(t2-t1);
            end 
        end 
        TTs(i).TTFoot = TTFoot;
    end 
    
    % XCorr - cross correlation calculation
    for i=1:numROIs
        ROIs2compare = allROIs ~= i;
        ROIs2compare = nonzeros(ROIs2compare.*allROIs); 
        Xcorr = [NaN,NaN,NaN];
        for j=1:length(ROIs2compare)
            flows(i) = normalize(flows(i),'range');
            iterator = ROIs2compare(j);
            if iterator>i
                XcorrPlot = xcorr(flows(i,:),flows(iterator,:));
                [~,maxXcorrIdx] = max(XcorrPlot);
                shift = -1*(maxXcorrIdx - length(flows(i,:)));
                Xcorr(iterator) = shift.*timeres;
            end 
        end 
        TTs(i).Xcorr = Xcorr;
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
                flow(count).Gaussian = handles.pcDatasets(i).ROIdataGaussian;
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
                    flow(count).Gaussian = handles.pcDatasets(i).ROIdataGaussian(j);
                    flow(count).Spline = handles.pcDatasets(i).ROIdataSpline(j);
                    count = count+1;
                end 
            end 
        end 
    end 

    
% --- Turn PolyLine into SplineLine    
function Y = interppolygon(X,N)
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
    
    
% --- Plot Velocities    
function plotVelocity(handles)
    global interpType showErrorBars
    cla(handles.VelocityPlot,'reset')
    axes(handles.VelocityPlot);
    
    legendSet = {'Baseline',handles.flow.Name};
    count = length(legendSet)-1;
    times = handles.flow(1).Data.times;
    
    plot(times, zeros(1,length(times)) ,'Color','black','LineWidth',1.5);
    for i=1:count
        switch interpType
            case 'None'
                flow = handles.flow(i).Data.meanROI;
                stdev = handles.flow(i).Data.stdROI;
            case 'Makima'
                times = handles.flow(i).Makima.times;
                flow = handles.flow(i).Makima.meanROI;
                stdev = handles.flow(i).Makima.stdROI;
            case 'Gaussian'
                times = handles.flow(i).Gaussian.times;
                flow = handles.flow(i).Gaussian.meanROI;  
                stdev = handles.flow(i).Gaussian.stdROI;
            case 'Spline'
                times = handles.flow(i).Spline.times;
                flow = handles.flow(i).Spline.meanROI; 
                stdev = handles.flow(i).Spline.stdROI;
            otherwise
        end
        
        if mean(flow)<0
            if showErrorBars
                hold on; errorbar(times,-1*flow,stdev);
            else
                hold on; plot(times,-1*flow);
            end 
        else 
            if showErrorBars
                hold on; errorbar(times,flow,stdev);
            else
                hold on; plot(times,flow);
            end
        end 
    end 
    legend(legendSet); hold off
    xlabel('Time (ms)'); ylabel('Mean Velocity in ROI (mm/s)');
    xlim([times(1),times(end)]);
    

% --- Sigmoid Fit Function
function [sigmoid,upslope,curvature,t1,t2,timeres] = sigFit(meanROI)
    [~,peak] = max(meanROI);
    upslope = meanROI(1:peak);
    t = 1:peak;
    times = 1:0.1:peak;
    t2 = length(times);
    upslope = interp1(t,upslope,times);
    upslope = normalize(upslope,'range');

    timeres = times(2)-times(1);

    %c1 = b, c2 = a, c3 = x0, c4 = dx
    sigmoidModel = @(c) c(1) + ( (c(2)-c(1)) ) ./ ( 1+exp((times-c(3))./c(4)) ) - upslope;
    c0 = [max(upslope),min(upslope),peak/2,timeres/2];
    opts = optimset('Display', 'off');
    c = lsqnonlin(sigmoidModel,c0,[],[],opts);
    sigmoid = c(1) + ( (c(2)-c(1)) ) ./ ( 1+exp((times-c(3))./c(4)) );

    for i=2:numel(upslope)-1
        dt = 0.5*(times(i+1)-times(i-1));
        ddt = 4*(times(i+1)-2*times(i)+times(i-1));
        dy = 0.5*(sigmoid(i+1)-sigmoid(i-1));
        ddy = 4*(sigmoid(i+1)-2*sigmoid(i)+sigmoid(i-1));
        curvature(i) = (ddy*dt - ddt*dy)./((dt^2 + dy^2).^(3/2));
    end 

    [~,t1] = max(curvature);
    
% --- Calculate derivative   
function fPrime = derivative(f)
    f(end+1) = f(1);
    fPrime = diff(f);
   
