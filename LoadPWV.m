function varargout = LoadPWV(varargin)
    % LOADPWV MATLAB code for LoadPWV.fig
    %      LOADPWV, by itself, creates a new LOADPWV or raises the existing
    %      singleton*.
    %
    %      H = LOADPWV returns the handle to a new LOADPWV or the handle to
    %      the existing singleton*.
    %
    %      LOADPWV('CALLBACK',hObject,eventData,handles,...) calls the local
    %      function named CALLBACK in LOADPWV.M with the given input arguments.
    %
    %      LOADPWV('Property','Value',...) creates a new LOADPWV or raises the
    %      existing singleton*.  Starting from the left, property value pairs are
    %      applied to the GUI before LoadPWV_OpeningFcn gets called.  An
    %      unrecognized property name or invalid value makes property application
    %      stop.  All inputs are passed to LoadPWV_OpeningFcn via varargin.
    %
    %      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
    %      instance to run (singleton)".
    %
    % See also: GUIDE, GUIDATA, GUIHANDLES

    % Edit the above text to modify the response to help LoadPWV

    % Last Modified by GUIDE v2.5 26-Oct-2019 10:49:51

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @LoadPWV_OpeningFcn, ...
                       'gui_OutputFcn',  @LoadPWV_OutputFcn, ...
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


% --- Executes just before LoadPWV is made visible.
function LoadPWV_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to LoadPWV (see VARARGIN)

    % Choose default command line output for LoadPWV
    handles.output = hObject;

    % Update handles structure
    guidata(hObject, handles);
    
    handles.global.preloaded = 0;
    handles.anatDatasets.Names = '';
    handles.anatDatasets.Data = 0;
    handles.global.anatDataIter = 1;    
    handles.pcDatasets.Names = '';
    handles.pcDatasets.Data = 0;
    handles.global.pcDataIter = 1;
    handles.magDatasets.Names = '';
    handles.magDatasets.Data = 0;
    handles.global.magDataIter = 1;
    handles.pcviprDatasets.Names = '';
    handles.pcviprDatasets.Data = 0;
    handles.global.pcviprDataIter = 1;

    guidata(hObject, handles);
    
% --- Outputs from this function are returned to the command line.
function varargout = LoadPWV_OutputFcn(hObject, eventdata, handles) 
    varargout{1} = handles.output;


    

    
%%%%%%%%%%%%%% PRELOAD CASE PANEL %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- PRELOAD CASE CALLBACK
function PreloadCaseButton_Callback(hObject, eventdata, handles)
    [casefile, casefileDir] = uigetfile({'*.mat;','Useable Files (*.mat)';
       '*.mat','MAT-files (*.mat)'; ...
       '*.*',  'All Files (*.*)'}, 'Select the case file (NAME_pwvcasefile_DATE.mat)');
    load([casefileDir '\' casefile]);
    
    handles.anatDatasets = anatDatasets;
    handles.pcDatasets = pcDatasets;
    handles.magDatasets = magDatasets;
    handles.pcviprDatasets = pcviprDatasets;
    
    set(handles.ListboxAnatomical,'String',{handles.anatDatasets.Names});
    set(handles.Listbox2DPC,'String',{handles.pcDatasets.Names});
    set(handles.Listbox2DMAG,'String',{handles.magDatasets.Names});
    set(handles.Listbox4DFlow,'String',{handles.pcviprDatasets.Names});

    set(handles.ErrorMessageBar,'String','Case file loaded successfully')

    handles.global.anatDataIter = numel(anatDatasets) + 1;
    handles.global.pcDataIter = numel(pcDatasets) + 1;
    handles.global.magDataIter = numel(magDatasets) + 1;
    handles.global.pcviprDataIter = numel(pcviprDatasets) + 1;
    
    clear casefile casefileDir    
    guidata(hObject, handles);


    

    
%%%%%%%%%%%%%% ANATOMICAL PANEL %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- DATA NAME ANATOMICAL CALLBACK
function DataNameAnatomical_Callback(hObject, eventdata, handles)

% --- DATA NAME ANATOMICAL CREATE FUNCTION
function DataNameAnatomical_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

    
% --- BUTTON LOAD ANATOMICAL CALLBACK
function ButtonLoadAnatomical_Callback(hObject, eventdata, handles)
    anatDataIter = handles.global.anatDataIter;
    set(handles.ErrorMessageBar,'String','')
    if isempty(get(handles.DataNameAnatomical,'String'))
        set(handles.ErrorMessageBar,'String','Name is required to load data! Type a name into the text box and press the Load Data button.')
    else 
        name = get(handles.DataNameAnatomical,'String');
        [anatomicalFile, anatomicalDir] = uigetfile({'*.dcm;*.mat;','Useable Files (*.dcm,*.mat)';
       '*.dcm',  'DICOM files (*.dcm)'; ...
       '*.mat','MAT-files (*.mat)'; ...
       '*.*',  'All Files (*.*)'}, 'Select ONE anatomical file in the dataset');
        [~,~,extension] = fileparts(anatomicalFile);
        dirInfo = dir(fullfile(anatomicalDir,['*' extension]));
        if isequal(extension,'.dcm')
            handles.anatDatasets(anatDataIter).Info = dicominfo(fullfile(anatomicalDir,dirInfo(1).name));
            for i=1:length(dirInfo)
                data(:,:,i) = single(dicomread(fullfile(anatomicalDir,dirInfo(i).name)));
            end        
        else
            temp = load(fullfile(anatomicalDir,dirInfo(1).name));
            data = single(struct2array(temp));
        end

        if sum(strcmp({handles.anatDatasets.Names},name))>0
            set(handles.ErrorMessageBar,'String','This name has already been used, try loading again with a different name')
            set(handles.DataNameAnatomical,'String','')
        else 
            set(handles.ErrorMessageBar,'String',['"' name '" has been loaded successfully. Add another dataset OR select "Loading Complete" if finished loading data']);
            handles.anatDatasets(anatDataIter).Names = name;
            handles.anatDatasets(anatDataIter).Data = data;
            handles.global.anatDataIter = anatDataIter + 1;
            set(handles.DataNameAnatomical,'String','')
            set(handles.ListboxAnatomical,'String',{handles.anatDatasets.Names});
        end 
    end 
    
    clear name temp data dirInfo extension anatomicalFile anatomicalDir i
    guidata(hObject, handles);

% --- BUTTON LOAD ANATOMICAL CREATE FUNCTION
function ButtonLoadAnatomical_CreateFcn(hObject, eventdata, handles)


% --- LISTBOX ANATOMICAL CALLBACK
function ListboxAnatomical_Callback(hObject, eventdata, handles)

% --- LISTBOX ANATOMICAL CREATE FUNCTION
function ListboxAnatomical_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


% --- REMOVE ANATOMICAL CALLBACK
function RemoveAnatomical_Callback(hObject, eventdata, handles)
    index = get(handles.ListboxAnatomical, 'value');
    handles.anatDatasets(index) = [];
    set(handles.ListboxAnatomical,'Value',1);
    set(handles.ListboxAnatomical,'String',{handles.anatDatasets.Names});
    handles.global.anatDataIter = handles.global.anatDataIter - 1;
    
    clear index
    guidata(hObject, handles);

% --- REMOVE ANATOMICAL CREATE FUNCTION
function RemoveAnatomical_CreateFcn(hObject, eventdata, handles)





%%%%%%%%%%%%%% 2DPC PANEL %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- 2D PC DATA NAME CALLBACK
function DataName2DPC_Callback(hObject, eventdata, handles)

% --- 2D PC DATA NAME CREATE FUNCTION
function DataName2DPC_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    
    
% --- BUTTON LOAD 2DPC CALLBACK
function ButtonLoad2DPC_Callback(hObject, eventdata, handles)
    pcDataIter = handles.global.pcDataIter;
    set(handles.ErrorMessageBar,'String','')
    if isempty(get(handles.DataName2DPC,'String'))
        set(handles.ErrorMessageBar,'String','Name is required to load data! Type a name into the text box and press the Load Data button.')
    elseif sum(strcmp({handles.pcDatasets.Names},get(handles.DataName2DPC,'String')))>0
        set(handles.ErrorMessageBar,'String','This name has already been used, try loading again with a different name')
        set(handles.DataName2DPC,'String','')
    else 
        name = get(handles.DataName2DPC,'String');
        [pcFile, pcDir] = uigetfile({'*.dcm;*.dat;','Useable Files (*.dcm,*.dat)';
       '*.dcm', 'DICOM files (*.dcm)'; ...
       '*.dat', 'Data files (*.dat)'; ...
       '*.*',  'All Files (*.*)'}, 'Select ONE 2D phase contrast file in the dataset');
        [~,~,extension] = fileparts(pcFile);
        dirInfo = dir(fullfile(pcDir,['*' extension]));
        if isequal(extension,'.dcm')
            handles.pcDatasets(pcDataIter).Info = dicominfo(fullfile(pcDir,dirInfo(1).name));
            handles.pcDatasets(pcDataIter).isDICOM = 1;
            for i=1:length(dirInfo)
                temp(:,:,i) = single(dicomread(fullfile(pcDir,dirInfo(i).name)));
            end   
            data.v = temp(:,:,1:floor(length(dirInfo)/2)); %%maybe change to velocity
            data.mag = temp(:,:,floor(length(dirInfo)/2)+1:end);
            VMEAN = BGPhaseCorrect(data.v);
        else
            fid = fopen([pcDir '\pcvipr_header.txt'], 'r');
            dataArray = textscan(fid, '%s%s%[^\n\r]', 'Delimiter', ' ', 'MultipleDelimsAsOne', true, 'ReturnOnError', false);
            fclose(fid);
            dataArray{1,2} = cellfun(@str2num,dataArray{1,2}(:), 'UniformOutput', false);
            pcviprHeader = cell2struct(dataArray{1,2}(:), dataArray{1,1}(:), 1);
            handles.pcDatasets(pcDataIter).Info = pcviprHeader;
            handles.pcDatasets(pcDataIter).isDICOM = 0;
            resx = pcviprHeader.matrixx;  
            resy = pcviprHeader.matrixy;  
            nframes = pcviprHeader.frames;      
            MAG = load_dat(fullfile(pcDir,'MAG.dat'),[resx resy]);
            CD = load_dat(fullfile(pcDir,'CD.dat'),[resx resy]);
            VMEAN = load_dat(fullfile(pcDir,'comp_vd_3.dat'),[resx resy]);
            v = zeros(resx,resy,nframes);
            mag = zeros(resx,resy,nframes);
            cd = zeros(resx,resy,nframes);
            %%% Find the non-zero velocity direction
            for j = 1:nframes   
                v(:,:,j) = load_dat(fullfile(pcDir, ['\ph_' num2str(j-1,'%03i') '_vd_3.dat']),[resx resy]);
                mag(:,:,j) = load_dat(fullfile(pcDir, ['\ph_' num2str(j-1,'%03i') '_mag.dat']),[resx resy]);
                cd(:,:,j) = load_dat(fullfile(pcDir, ['\ph_' num2str(j-1,'%03i') '_cd.dat']),[resx resy]);
            end 
            data.cd = cd;
            data.mag = mag;
            data.v = v;
            data.CD = CD;
            data.MAG = MAG;
            data.VMEAN = VMEAN;
        end
    set(handles.ErrorMessageBar,'String',['"' name '" has been loaded successfully. Add another dataset OR select "Loading Complete" if finished loading data']);
    handles.pcDatasets(pcDataIter).Names = name;
    handles.pcDatasets(pcDataIter).Data = data;
    handles.global.pcDataIter = pcDataIter + 1;
    set(handles.DataName2DPC,'String','')
    set(handles.Listbox2DPC,'String',{handles.pcDatasets.Names});
    end 
    
    clear name pcFile pcDir cd mag v CD MAG VMEAN resx resy pcviprHeader 
    clear dataArray nframes dirInfo temp i j fid ans extension data
    guidata(hObject, handles);

% --- BUTTON LOAD 2DPC CREATE FUNCTION
function ButtonLoad2DPC_CreateFcn(hObject, eventdata, handles)


% --- LISTBOX 2DPC CALLBACK
function Listbox2DPC_Callback(hObject, eventdata, handles)

% --- LISTBOX 2DPC CREATE FUNCTION
function Listbox2DPC_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


% --- REMOVE 2DPC CALLBACK
function Remove2DPC_Callback(hObject, eventdata, handles)
    index = get(handles.Listbox2DPC, 'value');
    handles.global.pcDatasets(index) = [];
    set(handles.Listbox2DPC,'Value',1);
    set(handles.Listbox2DPC,'String',{handles.global.pcDatasets.Names});
    handles.global.pcDataIter = handles.global.pcDataIter - 1;
    
    clear index
    guidata(hObject, handles);
        
% --- REMOVE 2DPC CREATE FUNCTION
function Remove2DPC_CreateFcn(hObject, eventdata, handles)





%%%%%%%%%%%%%% 2DMAG PANEL %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- DATA NAME 2DMAG CALLBACK
function DataName2DMAG_Callback(hObject, eventdata, handles)

% --- DATA NAME 2DMAG CREATE FUNCTION
function DataName2DMAG_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

    
% --- BUTTON LOAD 2DMAG CALLBACK
function ButtonLoad2DMAG_Callback(hObject, eventdata, handles)
    magDataIter = handles.global.magDataIter;
    set(handles.ErrorMessageBar,'String','')

    if isempty(get(handles.DataName2DMAG,'String'))
        set(handles.ErrorMessageBar,'String','Name is required to load data! Type a name into the text box and press the Load Data button.')
    else 
        name = get(handles.DataName2DMAG,'String');
        [magFile, magDir] = uigetfile({'*.dcm;*.mat;','Useable Files (*.dcm,*.mat)';
       '*.dcm',  'DICOM files (*.dcm)'; ...
       '*.mat','MAT-files (*.mat)'; ...
       '*.*',  'All Files (*.*)'}, 'Select ONE 2D magnitude file in the dataset');
        [~,~,extension] = fileparts(magFile);
        dirInfo = dir(fullfile(magDir,['*' extension]));
        
        if isequal(extension,'.dcm')
            handles.magDatasets(magDataIter).Info = dicominfo(fullfile(magDir,dirInfo(1).name));
            for i=1:length(dirInfo)
                data(:,:,i) = single(dicomread(fullfile(magDir,dirInfo(i).name)));
            end        
        else
            temp = load(fullfile(magDir,dirInfo(1).name));
            data = single(struct2array(temp));
        end
        
        if sum(strcmp({handles.magDatasets.Names},name))>0
            set(handles.ErrorMessageBar,'String','This name has already been used, try loading again with a different name')
            set(handles.DataName2DMAG,'String','')
        else 
            set(handles.ErrorMessageBar,'String',['"' name '" has been loaded successfully. Add another dataset OR select "Loading Complete" if finished loading data']);
            handles.magDatasets(magDataIter).Names = name;
            handles.magDatasets(magDataIter).Data = data;
            handles.global.magDataIter = magDataIter + 1;
            set(handles.DataName2DMAG,'String','')
            set(handles.Listbox2DMAG,'String',{handles.magDatasets.Names});
        end 
    end 
    
    clear name data temp dirInfo extension magFile magDir i
    guidata(hObject, handles);

% --- BUTTON LOAD 2DMAG CREATE FUNCTION
function ButtonLoad2DMAG_CreateFcn(hObject, eventdata, handles)


% --- LISTBOX 2DMAG CALLBACK
function Listbox2DMAG_Callback(hObject, eventdata, handles)

% --- LISTBOX 2DMAG CREATE FUNCTION
function Listbox2DMAG_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


% --- REMOVE 2DMAG CALLBACK
function Remove2DMAG_Callback(hObject, eventdata, handles)
    index = get(handles.Listbox2DMAG, 'value');
    handles.global.magDatasets(index) = [];
    set(handles.Listbox2DMAG,'Value',1);
    set(handles.Listbox2DMAG,'String',{handles.global.magDatasets.Names});
    handles.global.magDataIter = handles.global.magDataIter - 1;
    
    clear index 
    guidata(hObject, handles);

% --- REMOVE 2DMAG CREATE FUNCTION
function Remove2DMAG_CreateFcn(hObject, eventdata, handles)





%%%%%%%%%%%%%% 4D FLOW PANEL %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- DATA NAME 4D FLOW CALLBACK
function DataName4DFlow_Callback(hObject, eventdata, handles)

% --- DATA NAME 4D FLOW CREATE FUNCTION
function DataName4DFlow_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

    
% --- BUTTON LOAD 4D FLOW CALLBACK
function ButtonLoad4DFlow_Callback(hObject, eventdata, handles)
    pcviprDataIter = handles.global.pcviprDataIter;
    set(handles.ErrorMessageBar,'String','')

    if isempty(get(handles.DataName4DFlow,'String'))
        set(handles.ErrorMessageBar,'String','Name is required to load data! Type a name into the text box and press the Load Data button.')
    else 
        name = get(handles.DataName4DFlow,'String');
        [pcviprFile, pcviprDir] = uigetfile({'*.dcm;*.mat;','Useable Files (*.dcm,*.mat)';
       '*.dcm',  'DICOM files (*.dcm)'; ...
       '*.mat','MAT-files (*.mat)'; ...
       '*.*',  'All Files (*.*)'}, 'Select ONE 4D flow file in the dataset');
        [~,~,extension] = fileparts(pcviprFile);
        dirInfo = dir(fullfile(pcviprDir,['*' extension]));
        
        if isequal(extension,'.dcm')
            handles.pcviprDatasets(pcviprDataIter).Info = dicominfo(fullfile(pcDir,dirInfo(1).name));
            for i=1:length(dirInfo)
                data(:,:,i) = single(dicomread(fullfile(pcviprDir,dirInfo(i).name)));
            end        
        else
            temp = load(fullfile(pcviprDir,dirInfo(1).name));
            data = single(struct2array(temp));
        end
        
        if sum(strcmp({handles.pcviprDatasets.Names},name))>0
            set(handles.ErrorMessageBar,'String','This name has already been used, try loading again with a different name')
            set(handles.DataName4DFlow,'String','')
        else 
            set(handles.ErrorMessageBar,'String',['"' name '" has been loaded successfully. Add another dataset OR select "Loading Complete" if finished loading data']);
            handles.pcviprDatasets(pcviprDataIter).Names = name;
            handles.pcviprDatasets(pcviprDataIter).Data = data;
            handles.global.pcviprDataIter = pcviprDataIter + 1;
            set(handles.DataName4DFlow,'String','')
            set(handles.Listbox4DFlow,'String',{handles.pcviprDatasets.Names});
        end 
    end 
    
    clear name pcviprFile pcviprDir extension dirInfo temp data i 
    guidata(hObject, handles);

% --- BUTTON LOAD 4D FLOW CREATE FUNCTION
function ButtonLoad4DFlow_CreateFcn(hObject, eventdata, handles)


% --- LISTBOX 4D FLOW CALLBACK
function Listbox4DFlow_Callback(hObject, eventdata, handles)

% ---LISTBOX 4D FLOW CREATE FUNCTION
function Listbox4DFlow_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


% --- REMOVE 4D FLOW CALLBACK
function Remove4DFlow_Callback(hObject, eventdata, handles)
    index = get(handles.Listbox4DFlow, 'value');
    handles.global.pcviprDatasets(index) = [];
    set(handles.Listbox4DFlow,'Value',1);
    set(handles.Listbox4DFlow,'String',{handles.global.pcviprDatasets.Names});
    handles.global.pcviprDataIter = handles.global.pcviprDataIter - 1;
    
    clear index 
    guidata(hObject, handles);

% --- REMOVE 4D FLOW CREATE FUNCTION
function Remove4DFlow_CreateFcn(hObject, eventdata, handles)





%%%%%%%%%%%%%% POST LOADING %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- CASENAME CALLBACK
function CaseName_Callback(hObject, eventdata, handles)

% --- CASENAME CREATE FUNCTION
function CaseName_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

    
% --- SAVE CASE BUTTON CALLBACK
function SaveCaseButton_Callback(hObject, eventdata, handles)
    anatDatasets = handles.anatDatasets;
    pcDatasets = handles.pcDatasets;
    magDatasets = handles.magDatasets;
    pcviprDatasets = handles.pcviprDatasets;
    
    if isempty(get(handles.CaseName,'String'))
        set(handles.ErrorMessageBar,'String','Name is required to save data! Type a name into the text box and press the Save Case File button.')
    else
        set(handles.ErrorMessageBar, 'String', 'Select the directory where you will save the Case File.');
        casefileDir = uigetdir('C:\','Select the directory where you will save the Case File.');
        casename = get(handles.CaseName,'String');
        date = datestr(now);
        chopDate = [date(1:2) '-' date(4:6) '-' date(10:11) '-' date(13:14) date(16:17)];
        filename = [casename '_pwvcasefile_' chopDate '.mat'];
        set(handles.ErrorMessageBar,'String',['Saving Data as "' filename '"']);
        filenameWithPath = [casefileDir '\' filename];
        save(filenameWithPath,'anatDatasets','pcDatasets','magDatasets','pcviprDatasets');
        set(handles.ErrorMessageBar,'String','Data saved successfully');
        set(handles.CaseName,'');
    end 
    
    clear casefileDir casename date chopData filename filenameWithPath
    guidata(hObject, handles);

    
% --- ERROR MESSAGE BAR CALLBACK
function ErrorMessageBar_Callback(hObject, eventdata, handles)

% --- ERROR MESSAGE BAR CREATE FUNCTION
function ErrorMessageBar_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


% --- COMPLETE LOADING BUTTON CALLBACK
function ButtonCompleteLoading_Callback(hObject, eventdata, handles)
    if handles.global.anatDataIter==1
        set(handles.ErrorMessageBar,'String','Anatomical datasets required to complete loading process.');
    elseif handles.global.pcDataIter==1
        set(handles.ErrorMessageBar,'String','2D PC datasets required to complete loading process.');
    else
        set(handles.ErrorMessageBar,'String','Loading Process Complete...');
        AnalyzePWV(handles.anatDatasets,handles.pcDatasets,handles.magDatasets)
    end 

% --- COMPLETE LOADING BUTTON CREATE FUNCTION
function ButtonCompleteLoading_CreateFcn(hObject, eventdata, handles)





%%%%%%%%%%%%%% MY FUNCTIONS %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Dat files
function v = load_dat(name, res)

    [fid,errmsg]= fopen(name,'r');
    if fid < 0  % If name does not exist in directory
        set(handles.ErrorMessageBar,'String',['Error Opening Data : ',errmsg]);
    end

    % Reads in as short, reshapes by image res.
    v = reshape(fread(fid,'short=>single'),res);
    fclose(fid);
    clear fid errmsg 

% Perform 2D background phase correction
function vmean = BGPhaseCorrect(v)
vmean = mean(abs(v),3);
vmax = max(abs(v),[],3);

cdThresh = 0.10*max(vmax(:));
mask = vmax<cdThresh;
[rows,cols] = find(mask);
vmean = double(vmean.*mask);
vmeanVector = nonzeros(vmean(:));

f = fit( [rows, cols], vmeanVector, 'poly33' );
c = coeffvalues(f);

x0 = 1:size(mask,1);
y0 = 1:size(mask,2);
[y,x] = meshgrid(x0,y0);

fits = c(1) + c(2).*x + c(3).*y + c(4).*(x.^2) + c(5).*x.*y + c(6).*(y.^2) ...
     + c(7).*(x.^3) + c(8).*(x.^2).*y + c(9).*x.*(y.^2) + c(10).*(y.^3); 
vmean = vmean-fits;





%%%%%%%%%%%%%% UNUSED %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- MAIN PANEL CREATE FUNCTION
function MainPanel_CreateFcn(hObject, eventdata, handles)
% --- MAIN PANEL DELETE FUNCTION
function MainPanel_DeleteFcn(hObject, eventdata, handles)
% --- HEADER MAIN CREATE FUNCTION
function HeaderMain_CreateFcn(hObject, eventdata, handles)
% --- HEADER ANATOMICAL CREATEFUNCTION
function HeaderAnatomical_CreateFcn(hObject, eventdata, handles)
% --- HEADER ANATOMICAL LISTBOX CREATEFUNCTION
function HeaderAnatomicalListing_CreateFcn(hObject, eventdata, handles)
% --- LOAD ANATOMICAL TEXT CREATEFUNCTION
function LoadAnatomicalText_CreateFcn(hObject, eventdata, handles)
% --- 2D PC HEADER CREATE FUNCTION
function Header2DPC_CreateFcn(hObject, eventdata, handles)
% --- LOAD 2DPC TEXT CREATE FUNCTION
function Load2DPCText_CreateFcn(hObject, eventdata, handles)
% --- HEADER 2DPC LISTBOX CREATE FUNCTION
function Header2DPCListing_CreateFcn(hObject, eventdata, handles)
% --- HEADER 2DMAG CREATE FUNCTION
function Header2DMAG_CreateFcn(hObject, eventdata, handles)
% --- LOAD 2DMAG TEXT CREATE FUNCTION
function Load2DMAGText_CreateFcn(hObject, eventdata, handles)
% --- HEADER 2DMAG LISTING CREATE FUNCTION
function Header2DMAGListing_CreateFcn(hObject, eventdata, handles)
% --- HEADER 4DFLOW CREATE FUNCTION
function Header4DFlow_CreateFcn(hObject, eventdata, handles)
% ---LOAD 4D FLOW TEXT CREATE FUNCTION
function Load4DFlowText_CreateFcn(hObject, eventdata, handles)
% --- HEADER 4D FLOW LISTBOX
function Header4DFlowListing_CreateFcn(hObject, eventdata, handles)
