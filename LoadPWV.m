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

% Last Modified by GUIDE v2.5 10-Oct-2019 15:08:11

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

    % UIWAIT makes LoadPWV wait for user response (see UIRESUME)
    % uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = LoadPWV_OutputFcn(hObject, eventdata, handles) 
    % varargout  cell array for returning output args (see VARARGOUT);
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure
    varargout{1} = handles.output;


%%%%%%%%%%%%%% PRELOAD CASE PANEL %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- PRELOAD CASE CALLBACK
function PreloadCaseButton_Callback(hObject, eventdata, handles)
    global anatDataIter pcDataIter magDataIter pcviprDataIter preloaded

    [casefile, casefileDir] = uigetfile({'*.mat;','Useable Files (*.mat)';
       '*.mat','MAT-files (*.mat)'; ...
       '*.*',  'All Files (*.*)'}, 'Select the case file (NAME_pwvcasefile_DATE.mat)');
    load([casefileDir '\' casefile]);
    set(handles.ListboxAnatomical,'String',{anatDatasets.Names});
    set(handles.Listbox2DPC,'String',{pcDatasets.Names});
    set(handles.Listbox2DMAG,'String',{magDatasets.Names});
    set(handles.Listbox4DFlow,'String',{pcviprDatasets.Names});

    disp('Case file loaded successfully');
    set(handles.ErrorMessageBar,'String','Case file loaded successfully')

    anatDataIter = numel(anatDatasets);
    pcDataIter = numel(pcDatasets);
    magDataIter = numel(magDatasets);
    pcviprDataIter = numel(pcviprDatasets);
    preloaded=1;


    
    
%%%%%%%%%%%%%% ANATOMICAL PANEL %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- DATA NAME ANATOMICAL CALLBACK
function DataNameAnatomical_Callback(hObject, eventdata, handles)

% --- DATA NAME ANATOMICAL CREATE FUNCTION
function DataNameAnatomical_CreateFcn(hObject, eventdata, handles)
    global anatDatasets anatDataIter preloaded
    
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

    if isempty(preloaded)
        anatDatasets.Names = '';
        anatDatasets.Data = 0;
        anatDataIter = 1;
    end 

    
% --- BUTTON LOAD ANATOMICAL CALLBACK
function ButtonLoadAnatomical_Callback(hObject, eventdata, handles)
    global anatDatasets anatDataIter

    set(handles.ErrorMessageBar,'String','')
    if isempty(get(handles.DataNameAnatomical,'String'))
        disp('Name is required to load data!');
        disp('Type a name into the text box and press the Load Data button.');
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
            anatDatasets(anatDataIter).Info = dicominfo(fullfile(anatomicalDir,dirInfo(1).name));
            for i=1:length(dirInfo)
                data(:,:,i) = single(dicomread(fullfile(anatomicalDir,dirInfo(i).name)));
            end        
        else
            temp = load(fullfile(anatomicalDir,dirInfo(1).name));
            data = single(struct2array(temp));
        end

        if sum(strcmp({anatDatasets.Names},name))>0
            set(handles.ErrorMessageBar,'String','This name has already been used, try loading again with a different name')
            set(handles.DataNameAnatomical,'String','')
        else 
            set(handles.ErrorMessageBar,'String',['"' name '" has been loaded successfully. Add another dataset OR select "Loading Complete" if finished loading data']);
            anatDatasets(anatDataIter).Names = name;
            anatDatasets(anatDataIter).Data = data;
            anatDataIter = anatDataIter+1;
            set(handles.DataNameAnatomical,'String','')
            set(handles.ListboxAnatomical,'String',{anatDatasets.Names});
        end 
    end 

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
    global anatDatasets anatDataIter
    index = get(handles.ListboxAnatomical, 'value');
    anatDatasets(index) = [];
    set(handles.ListboxAnatomical,'Value',1);
    set(handles.ListboxAnatomical,'String',{anatDatasets.Names});
    anatDataIter = anatDataIter-1;

% --- REMOVE ANATOMICAL CREATE FUNCTION
function RemoveAnatomical_CreateFcn(hObject, eventdata, handles)




%%%%%%%%%%%%%% 2DPC PANEL %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- 2D PC DATA NAME CALLBACK
function DataName2DPC_Callback(hObject, eventdata, handles)

% --- 2D PC DATA NAME CREATE FUNCTION
function DataName2DPC_CreateFcn(hObject, eventdata, handles)
    global pcDatasets pcDataIter preloaded
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

    if isempty(preloaded)
        pcDatasets.Names = '';
        pcDatasets.Data = 0;
        pcDataIter = 1;
    end
    
    
% --- BUTTON LOAD 2DPC CALLBACK
function ButtonLoad2DPC_Callback(hObject, eventdata, handles)
    global pcDatasets pcDataIter
    set(handles.ErrorMessageBar,'String','')
    if isempty(get(handles.DataName2DPC,'String'))
        disp('Name is required to load data!');
        disp('Type a name into the text box and press the Load Data button.');
        set(handles.ErrorMessageBar,'String','Name is required to load data! Type a name into the text box and press the Load Data button.')
    elseif sum(strcmp({pcDatasets.Names},get(handles.DataName2DPC,'String')))>0
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
            pcDatasets(pcDataIter).Info = dicominfo(fullfile(pcDir,dirInfo(1).name));
            for i=1:length(dirInfo)
                temp(:,:,i) = single(dicomread(fullfile(pcDir,dirInfo(i).name)));
            end   
            data.v = temp(:,:,1:floor(length(dirInfo)/2)); %%maybe change to velocity
            data.mag = temp(:,:,floor(length(dirInfo)/2)+1:end);
        else
            fid = fopen([pcDir '\pcvipr_header.txt'], 'r');
            dataArray = textscan(fid, '%s%s%[^\n\r]', 'Delimiter', ' ', 'MultipleDelimsAsOne', true, 'ReturnOnError', false);
            fclose(fid);
            dataArray{1,2} = cellfun(@str2num,dataArray{1,2}(:), 'UniformOutput', false);
            pcviprHeader = cell2struct(dataArray{1,2}(:), dataArray{1,1}(:), 1);
            pcDatasets(pcDataIter).Info = pcviprHeader;
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
    pcDatasets(pcDataIter).Names = name;
    pcDatasets(pcDataIter).Data = data;
    pcDataIter = pcDataIter+1;
    set(handles.DataName2DPC,'String','')
    set(handles.Listbox2DPC,'String',{pcDatasets.Names});
    end 

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
    global pcDatasets pcDataIter
    index = get(handles.Listbox2DPC, 'value');
    pcDatasets(index) = [];
    set(handles.Listbox2DPC,'Value',1);
    set(handles.Listbox2DPC,'String',{pcDatasets.Names});
    pcDataIter = pcDataIter-1;

% --- REMOVE 2DPC CREATE FUNCTION
function Remove2DPC_CreateFcn(hObject, eventdata, handles)




%%%%%%%%%%%%%% 2DMAG PANEL %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- DATA NAME 2DMAG CALLBACK
function DataName2DMAG_Callback(hObject, eventdata, handles)

% --- DATA NAME 2DMAG CREATE FUNCTION
function DataName2DMAG_CreateFcn(hObject, eventdata, handles)
    global magDatasets magDataIter preloaded
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

    if isempty(preloaded)
        magDatasets.Names = '';
        magDatasets.Data = 0;
        magDataIter = 1;   
    end 

    
% --- BUTTON LOAD 2DMAG CALLBACK
function ButtonLoad2DMAG_Callback(hObject, eventdata, handles)
    global magDatasets magDataIter
    set(handles.ErrorMessageBar,'String','')

    if isempty(get(handles.DataName2DMAG,'String'))
        disp('Name is required to load data!');
        disp('Type a name into the text box and press the Load Data button.');
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
            magDatasets(magDataIter).Info = dicominfo(fullfile(magDir,dirInfo(1).name));
            for i=1:length(dirInfo)
                data(:,:,i) = single(dicomread(fullfile(magDir,dirInfo(i).name)));
            end        
        else
            temp = load(fullfile(magDir,dirInfo(1).name));
            data = single(struct2array(temp));
        end
        
        if sum(strcmp({magDatasets.Names},name))>0
            set(handles.ErrorMessageBar,'String','This name has already been used, try loading again with a different name')
            set(handles.DataName2DMAG,'String','')
        else 
            set(handles.ErrorMessageBar,'String',['"' name '" has been loaded successfully. Add another dataset OR select "Loading Complete" if finished loading data']);
            magDatasets(magDataIter).Names = name;
            magDatasets(magDataIter).Data = data;
            magDataIter = magDataIter+1;
            set(handles.DataName2DMAG,'String','')
            set(handles.Listbox2DMAG,'String',{magDatasets.Names});
        end 
        
    end 

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
    global magDatasets magDataIter
    index = get(handles.Listbox2DMAG, 'value');
    magDatasets(index) = [];
    set(handles.Listbox2DMAG,'Value',1);
    set(handles.Listbox2DMAG,'String',{magDatasets.Names});
    magDataIter = magDataIter-1;

% --- REMOVE 2DMAG CREATE FUNCTION
function Remove2DMAG_CreateFcn(hObject, eventdata, handles)




%%%%%%%%%%%%%% 4D FLOW PANEL %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- DATA NAME 4D FLOW CALLBACK
function DataName4DFlow_Callback(hObject, eventdata, handles)

% --- DATA NAME 4D FLOW CREATE FUNCTION
function DataName4DFlow_CreateFcn(hObject, eventdata, handles)
    global pcviprDatasets pcviprDataIter preloaded
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

    if isempty(preloaded)
        pcviprDatasets.Names = '';
        pcviprDatasets.Data = 0;
        pcviprDataIter = 1;
    end

    
% --- BUTTON LOAD 4D FLOW CALLBACK
function ButtonLoad4DFlow_Callback(hObject, eventdata, handles)
    global pcviprDatasets pcviprDataIter
    set(handles.ErrorMessageBar,'String','')

    if isempty(get(handles.DataName4DFlow,'String'))
        disp('Name is required to load data!');
        disp('Type a name into the text box and press the Load Data button.');
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
            pcviprDatasets(pcviprDataIter).Info = dicominfo(fullfile(pcDir,dirInfo(1).name));
            for i=1:length(dirInfo)
                data(:,:,i) = single(dicomread(fullfile(pcviprDir,dirInfo(i).name)));
            end        
        else
            temp = load(fullfile(pcviprDir,dirInfo(1).name));
            data = single(struct2array(temp));
        end
        
        if sum(strcmp({pcviprDatasets.Names},name))>0
            set(handles.ErrorMessageBar,'String','This name has already been used, try loading again with a different name')
            set(handles.DataName4DFlow,'String','')
        else 
            set(handles.ErrorMessageBar,'String',['"' name '" has been loaded successfully. Add another dataset OR select "Loading Complete" if finished loading data']);
            pcviprDatasets(pcviprDataIter).Names = name;
            pcviprDatasets(pcviprDataIter).Data = data;
            pcviprDataIter = pcviprDataIter+1;
            set(handles.DataName4DFlow,'String','')
            set(handles.Listbox4DFlow,'String',{pcviprDatasets.Names});
        end 
        
    end 

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
    global pcviprDatasets pcviprDataIter
    index = get(handles.Listbox4DFlow, 'value');
    pcviprDatasets(index) = [];
    set(handles.Listbox4DFlow,'Value',1);
    set(handles.Listbox4DFlow,'String',{pcviprDatasets.Names});
    pcviprDataIter = pcviprDataIter-1;

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
    global anatDatasets pcDatasets magDatasets pcviprDatasets
    if isempty(get(handles.CaseName,'String'))
        disp('Name is required to save data!');
        disp('Type a name into the text box and press the Save Case File button.');
        set(handles.ErrorMessageBar,'String','Name is required to save data! Type a name into the text box and press the Save Case File button.')
    else
        disp('Select the directory where you will save the Case File.')
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
    end 

    
% --- ERROR MESSAGE BAR CALLBACK
function ErrorMessageBar_Callback(hObject, eventdata, handles)

% --- ERROR MESSAGE BAR CREATE FUNCTION
function ErrorMessageBar_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


% --- COMPLETE LOADING BUTTON CALLBACK
function ButtonCompleteLoading_Callback(hObject, eventdata, handles)
    global anatDataIter pcDataIter
    global anatDatasets magDatasets pcDatasets
    if anatDataIter==1
        disp('Anatomical datasets required to complete loading process.');
        set(handles.ErrorMessageBar,'String','Anatomical datasets required to complete loading process.');
    elseif pcDataIter==1
        disp('2D PC datasets required to complete loading process.');
        set(handles.ErrorMessageBar,'String','2D PC datasets required to complete loading process.');
    else
        disp('Loading Process Complete...');
        set(handles.ErrorMessageBar,'String','Loading Process Complete...');
        AnalyzePWV(anatDatasets,pcDatasets,magDatasets)
        %AnalyzePWV_lite(anatDatasets,pcDatasets,magDatasets)
    end 

% --- COMPLETE LOADING BUTTON CREATE FUNCTION
function ButtonCompleteLoading_CreateFcn(hObject, eventdata, handles)




%%%%%%%%%%%%%% MY FUNCTIONS %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = load_dat(name, res)
% Load Dat
% Loads in dat files in current directory.
    [fid,errmsg]= fopen(name,'r');
    if fid < 0  % If name does not exist in directory
        disp(['Error Opening Data : ',errmsg]);
    end

    % Reads in as short, reshapes by image res.
    v = reshape(fread(fid,'short=>single'),res);
    v = imrotate(v,180);
    fclose(fid);


    
    
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