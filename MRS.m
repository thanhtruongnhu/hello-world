function varargout = MRS(varargin)
% MRS MATLAB code for MRS.fig
%      MRS, by itself, creates a new MRS or raises the existing
%      singleton*.
%
%      H = MRS returns the handle to a new MRS or the handle to
%      the existing singleton*.
%
%      MRS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MRS.M with the given input arguments.
%
%      MRS('Property','Value',...) creates a new MRS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MRS_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MRS_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MRS

% Last Modified by GUIDE v2.5 16-Aug-2019 16:40:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MRS_OpeningFcn, ...
                   'gui_OutputFcn',  @MRS_OutputFcn, ...
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


% --- Executes just before MRS is made visible.
function MRS_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reservded - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MRS (see VARARGIN)
clear global
clc;
% Choose default command line output for MRS
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MRS wait for user response (see UIRESUME)
% uiwait(handles.figure1);

global idobjname pxobjname pyobjname orobjname blobjname editother_obj editrunstep_obj phi_obj idc_obj 
global fid numberofRobots fid_data fr td_obj tp_obj

matlab_client('connect')
fid=[];                                    % handles file
fr=[];
fid_data=[];
numberofRobots=0; 

td_obj=handles.Td;
tp_obj=handles.Tp;

% get edit object names to indicate id of robots
id1_obj=handles.id1;
id2_obj=handles.id2;
id3_obj=handles.id3;
id4_obj=handles.id4;
id5_obj=handles.id5;
id6_obj=handles.id6;
id7_obj=handles.id7;
id8_obj=handles.id8;
id9_obj=handles.id9;
id10_obj=handles.id10;
id11_obj=handles.id11;
id12_obj=handles.id12;
id13_obj=handles.id13;
id14_obj=handles.id14;
id15_obj=handles.id15;
id16_obj=handles.id16;
id17_obj=handles.id17;
id18_obj=handles.id18;
id19_obj=handles.id19;
id20_obj=handles.id20;

id21_obj=handles.id21;
id22_obj=handles.id22;
id23_obj=handles.id23;
id24_obj=handles.id24;
id25_obj=handles.id25;
id26_obj=handles.id26;
id27_obj=handles.id27;
id28_obj=handles.id28;
id29_obj=handles.id29;
id30_obj=handles.id30;
id31_obj=handles.id31;
id32_obj=handles.id32;
id33_obj=handles.id33;
id34_obj=handles.id34;
id35_obj=handles.id35;
id36_obj=handles.id36;
id37_obj=handles.id37;
id38_obj=handles.id38;
id39_obj=handles.id39;
id40_obj=handles.id40;


idobjname=[id1_obj,id2_obj,id3_obj,id4_obj,id5_obj,id6_obj,id7_obj,id8_obj,id9_obj,id10_obj,id11_obj,id12_obj,id13_obj,id14_obj,id15_obj,id16_obj,id17_obj,id18_obj,id19_obj,id20_obj,...
           id21_obj,id22_obj,id23_obj,id24_obj,id25_obj,id26_obj,id27_obj,id28_obj,id29_obj,id30_obj,id31_obj,id32_obj,id33_obj,id34_obj,id35_obj,id36_obj,id37_obj,id38_obj,id39_obj,id40_obj];

% get edit object names to indicate position x of robots
px1_obj=handles.px1;
px2_obj=handles.px2;
px3_obj=handles.px3;
px4_obj=handles.px4;
px5_obj=handles.px5;
px6_obj=handles.px6;
px7_obj=handles.px7;
px8_obj=handles.px8;
px9_obj=handles.px9;
px10_obj=handles.px10;
px11_obj=handles.px11;
px12_obj=handles.px12;
px13_obj=handles.px13;
px14_obj=handles.px14;
px15_obj=handles.px15;
px16_obj=handles.px16;
px17_obj=handles.px17;
px18_obj=handles.px18;
px19_obj=handles.px19;
px20_obj=handles.px20;

px21_obj=handles.px21;
px22_obj=handles.px22;
px23_obj=handles.px23;
px24_obj=handles.px24;
px25_obj=handles.px25;
px26_obj=handles.px26;
px27_obj=handles.px27;
px28_obj=handles.px28;
px29_obj=handles.px29;
px30_obj=handles.px30;
px31_obj=handles.px31;
px32_obj=handles.px32;
px33_obj=handles.px33;
px34_obj=handles.px34;
px35_obj=handles.px35;
px36_obj=handles.px36;
px37_obj=handles.px37;
px38_obj=handles.px38;
px39_obj=handles.px39;
px40_obj=handles.px40;


pxobjname=[px1_obj,px2_obj,px3_obj,px4_obj,px5_obj,px6_obj,px7_obj,px8_obj,px9_obj,px10_obj,px11_obj,px12_obj,px13_obj,px14_obj,px15_obj,px16_obj,px17_obj,px18_obj,px19_obj,px20_obj,...
           px21_obj,px22_obj,px23_obj,px24_obj,px25_obj,px26_obj,px27_obj,px28_obj,px29_obj,px30_obj,px31_obj,px32_obj,px33_obj,px34_obj,px35_obj,px36_obj,px37_obj,px38_obj,px39_obj,px40_obj];


% get edit object names to indicate position y of robots
py1_obj=handles.py1;
py2_obj=handles.py2;
py3_obj=handles.py3;
py4_obj=handles.py4;
py5_obj=handles.py5;
py6_obj=handles.py6;
py7_obj=handles.py7;
py8_obj=handles.py8;
py9_obj=handles.py9;
py10_obj=handles.py10;
py11_obj=handles.py11;
py12_obj=handles.py12;
py13_obj=handles.py13;
py14_obj=handles.py14;
py15_obj=handles.py15;
py16_obj=handles.py16;
py17_obj=handles.py17;
py18_obj=handles.py18;
py19_obj=handles.py19;
py20_obj=handles.py20;

py21_obj=handles.py21;
py22_obj=handles.py22;
py23_obj=handles.py23;
py24_obj=handles.py24;
py25_obj=handles.py25;
py26_obj=handles.py26;
py27_obj=handles.py27;
py28_obj=handles.py28;
py29_obj=handles.py29;
py30_obj=handles.py30;
py31_obj=handles.py31;
py32_obj=handles.py32;
py33_obj=handles.py33;
py34_obj=handles.py34;
py35_obj=handles.py35;
py36_obj=handles.py36;
py37_obj=handles.py37;
py38_obj=handles.py38;
py39_obj=handles.py39;
py40_obj=handles.py40;

pyobjname=[py1_obj,py2_obj,py3_obj,py4_obj,py5_obj,py6_obj,py7_obj,py8_obj,py9_obj,py10_obj,py11_obj,py12_obj,py13_obj,py14_obj,py15_obj,py16_obj,py17_obj,py18_obj,py19_obj,py20_obj,...
           py21_obj,py22_obj,py23_obj,py24_obj,py25_obj,py26_obj,py27_obj,py28_obj,py29_obj,py30_obj,py31_obj,py32_obj,py33_obj,py34_obj,py35_obj,py36_obj,py37_obj,py38_obj,py39_obj,py40_obj];

% get edit object names to indicate orientation of robots
or1_obj=handles.or1;
or2_obj=handles.or2;
or3_obj=handles.or3;
or4_obj=handles.or4;
or5_obj=handles.or5;
or6_obj=handles.or6;
or7_obj=handles.or7;
or8_obj=handles.or8;
or9_obj=handles.or9;
or10_obj=handles.or10;
or11_obj=handles.or11;
or12_obj=handles.or12;
or13_obj=handles.or13;
or14_obj=handles.or14;
or15_obj=handles.or15;
or16_obj=handles.or16;
or17_obj=handles.or17;
or18_obj=handles.or18;
or19_obj=handles.or19;
or20_obj=handles.or20;

or21_obj=handles.or21;
or22_obj=handles.or22;
or23_obj=handles.or23;
or24_obj=handles.or24;
or25_obj=handles.or25;
or26_obj=handles.or26;
or27_obj=handles.or27;
or28_obj=handles.or28;
or29_obj=handles.or29;
or30_obj=handles.or30;
or31_obj=handles.or31;
or32_obj=handles.or32;
or33_obj=handles.or33;
or34_obj=handles.or34;
or35_obj=handles.or35;
or36_obj=handles.or36;
or37_obj=handles.or37;
or38_obj=handles.or38;
or39_obj=handles.or39;
or40_obj=handles.or40;

orobjname=[or1_obj,or2_obj,or3_obj,or4_obj,or5_obj,or6_obj,or7_obj,or8_obj,or9_obj,or10_obj,or11_obj,or12_obj,or13_obj,or14_obj,or15_obj,or16_obj,or17_obj,or18_obj,or19_obj,or20_obj,...
           or21_obj,or22_obj,or23_obj,or24_obj,or25_obj,or26_obj,or27_obj,or28_obj,or29_obj,or30_obj,or31_obj,or32_obj,or33_obj,or34_obj,or35_obj,or36_obj,or37_obj,or38_obj,or39_obj,or40_obj];


% get edit object names to indicate battery level of robots
bl1_obj=handles.bl1;
bl2_obj=handles.bl2;
bl3_obj=handles.bl3;
bl4_obj=handles.bl4;
bl5_obj=handles.bl5;
bl6_obj=handles.bl6;
bl7_obj=handles.bl7;
bl8_obj=handles.bl8;
bl9_obj=handles.bl9;
bl10_obj=handles.bl10;
bl11_obj=handles.bl11;
bl12_obj=handles.bl12;
bl13_obj=handles.bl13;
bl14_obj=handles.bl14;
bl15_obj=handles.bl15;
bl16_obj=handles.bl16;
bl17_obj=handles.bl17;
bl18_obj=handles.bl18;
bl19_obj=handles.bl19;
bl20_obj=handles.bl20;

bl21_obj=handles.bl21;
bl22_obj=handles.bl22;
bl23_obj=handles.bl23;
bl24_obj=handles.bl24;
bl25_obj=handles.bl25;
bl26_obj=handles.bl26;
bl27_obj=handles.bl27;
bl28_obj=handles.bl28;
bl29_obj=handles.bl29;
bl30_obj=handles.bl30;
bl31_obj=handles.bl31;
bl32_obj=handles.bl32;
bl33_obj=handles.bl33;
bl34_obj=handles.bl34;
bl35_obj=handles.bl35;
bl36_obj=handles.bl36;
bl37_obj=handles.bl37;
bl38_obj=handles.bl38;
bl39_obj=handles.bl39;
bl40_obj=handles.bl40;


blobjname=[bl1_obj,bl2_obj,bl3_obj,bl4_obj,bl5_obj,bl6_obj,bl7_obj,bl8_obj,bl9_obj,bl10_obj,bl11_obj,bl12_obj,bl13_obj,bl14_obj,bl15_obj,bl16_obj,bl17_obj,bl18_obj,bl19_obj,bl20_obj,...
           bl21_obj,bl22_obj,bl23_obj,bl24_obj,bl25_obj,bl26_obj,bl27_obj,bl28_obj,bl29_obj,bl30_obj,bl31_obj,bl32_obj,bl33_obj,bl34_obj,bl35_obj,bl36_obj,bl37_obj,bl38_obj,bl39_obj,bl40_obj];

editother_obj=handles.editother;
editrunstep_obj=handles.editrunstep;
phi_obj=handles.phi;
idc_obj=handles.idc;


% --- Executes on button press in sinit_button.
function sinit_button_Callback(hObject, eventdata, handles)
% hObject    handle to sinit_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hf w corner g

hf=handles.axes1;
hold on;
grid off;
cla(hf);
clc;

% % Add code for showing walls here.
[FileName, PathName] = uigetfile('*.txt');
Name = fullfile(PathName,FileName);
if PathName==0,
    return; 
end
w=load(Name);

corner=[];
for i=1:size(w,1),
    if i==1
        A=[w(i,1),w(i,2)];
        B=[w(i,3),w(i,4)];
        corner=[corner;A;B];
    else
        A=[w(i,1),w(i,2)];
        B=[w(i,3),w(i,4)];
        statusA=1;
        for j=1:size(corner,1),
            if norm(A-corner(j,:))==0
                statusA=0;
                break;
            end
        end
        if statusA==1
            corner=[corner;A];
        end
        statusB=1;
        for j=1:size(corner,1),
            if norm(B-corner(j,:))==0
                statusB=0;
                break;
            end
        end
        if statusB==1
            corner=[corner;B];
        end
    end
    plot([w(i,1) w(i,3)],[w(i,2),w(i,4)],'-b');
end
%xv=[1 -2;2 2];

%% Add code for showing gates here.
[FileName, PathName] = uigetfile('*.txt');
Name = fullfile(PathName,FileName);
if PathName==0,
    return; 
end
g=load(Name);
for i=1:size(g,1),
    plot([g(i,1) g(i,3)],[g(i,2),g(i,4)],'-.k');
    text(double((g(i,1)+g(i,3))/2),double((g(i,2)+g(i,4))/2),num2str(i),'HorizontalAlignment','center','Fontsize',8,'color','r');
end

% --- Executes on button press in RobotInit.
function RobotInit_Callback(hObject, eventdata, handles)
% hObject    handle to RobotInit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hf mode comobj natnetclient 
global rc ra rh epsilon beta phi Tc L rw Lb vmax wmax Pr
global numberofRobots xv numofvic numofHop
global G Robot IDSet editother_obj fid_data count fid Hop swarm_obj
global x3 vmobj to

% Receiving control Parameters
rc=str2num(get(handles.rc,'string'));
ra=str2num(get(handles.ra,'string'));
rh=str2num(get(handles.rh,'string'));
epsilon=str2num(get(handles.epsilon,'string'));
beta=str2num(get(handles.beta,'string'));
phi=str2num(get(handles.phi,'string'));
Tc=str2num(get(handles.tc,'string'));



% Body parameter of robot
L=0.1;      % distance between two wheel
rw=0.02;    % wheel radius
Lb=0.085;   % body lenght
vmobj=handles.vm;   
vmax=str2num(get(vmobj,'string'));%1.25;  % Average rate of PWM = 62.5 when theta=0
wmax=59;   % corresponding to dtheta/dt=pi/4
Pr=70;
numofHop=0;
count=0;

fid=[];                                    % handles file
fid_data=[];

xv=[];      % set of victim points

mode=get(handles.modeselection,'value');   % mode for experiment: "1"=real experiment; "0"=simulation
x3=[];
switch mode
    case 1
        C = matlab_client('getsquares');
        IDSet=[];
        for i=1:C.dat(C.n),
            k=8*(i-1);
            ID=C.dat(k+1);
            if ID<=41
                if ID<=40   % for children robots
                    y=C.dat(k+2);    % x
                    z=C.dat(k+3);    % y
                    x=C.dat(k+4);    % z
                    q=[C.dat(k+5),C.dat(k+6),C.dat(k+7),C.dat(k+8)]; %[qx,qy,qz,qw]
                    IDSet=[IDSet,ID];
                    [yaw,roll,pitch]=quat2angle(q,'yxz');%'yxz'
                    j=size(IDSet,2);
%                     Robot(j).x=[x,y];
%                     Robot(j).angle=yaw;                    
%                     Robot(j).body=plot(Robot(j).x(1),Robot(j).x(2),'or','Markersize',3);
%                     Robot(j).LED=plot(Robot(j).x(1),Robot(j).x(2),'o','MarkerFaceColor','w');
%                     Robot(j).idtext=text(double(Robot(j).x(1)-0.1),double(Robot(j).x(2)-0.1),num2str(IDSet(j)),'HorizontalAlignment','center','Fontsize',8,'color','r');
%                     initparameters(j);
                    RobotInit(j,x,y,yaw); % Robot i is assigned with ID=IDSet(i).
                else    % for mother robot
                    delta=[0.3675,0.00289]; %[0.3313,0.00289]; 
                    y=C.dat(k+2);    % x
                    z=C.dat(k+3);    % y
                    x=C.dat(k+4);    % z
                    q=[C.dat(k+5),C.dat(k+6),C.dat(k+7),C.dat(k+8)]; %[qx,qy,qz,qw]
                    IDSet=[IDSet,ID];
                    [yaw,roll,pitch]=quat2angle(q,'yxz');%'yxz'
                    j=size(IDSet,2);
                    R=[cos(yaw),sin(yaw);-sin(yaw),cos(yaw)]; % Counterclockwise Rotation Matrix
                    idm=j;
                    Pos=[x,y]+delta*R;
                    MRobotInit(j,Pos(1,1),Pos(1,2),yaw); % Robot i is assigned with ID=IDSet(i).
                    hm=plot(hf,Robot(idm).x(1),Robot(idm).x(2),'.r');                    
                end
            else
                if ID<=50   % for victims
                    y=C.dat(k+2);    % x
                    x=C.dat(k+4);    % z
                    xv=[xv;[x,y]];
                else        % stick for human operator
                    y=C.dat(k+2);    % x
                    z=C.dat(k+3);    % y
                    x=C.dat(k+4);    % z
                    q=[C.dat(k+5),C.dat(k+6),C.dat(k+7),C.dat(k+8)]; %[qx,qy,qz,qw]
                    [yaw,roll,pitch]=quat2angle(q,'yxz');%'yxz'
                    Hop.x=[x,y,z];
                    Hop.angle=[yaw,roll,pitch];
                    numofHop=1;
                end
            end
        end

        comobj = instrfind('Type', 'serial', 'Port', 'COM5', 'Tag', '');
        if isempty(comobj)
            comobj = serial('COM5');
        else
            fclose(comobj);
            comobj = comobj(1);
        end        
        
        comobj.OutputBufferSize = 1024;
        comobj.BaudRate = 38400;
        comobj.Timeout=5; %2
        fopen(comobj);
        

    case 2
        [FileName, PathName] = uigetfile('*.txt');
        Name = fullfile(PathName,FileName);
        if PathName==0,
            return; 
        end
        xr=load(Name);
        IDSet=[1:1:size(xr,1)];
        for i=1:size(xr,1),
            x=xr(i,1);
            y=xr(i,2);
            angle=randn(1)*pi;%xr(i,3);%randn(1)*pi;
            RobotInit(i,x,y,angle);
            x3=[x3,Robot(i).angle];
        end
        xv=[];%[1 -2.5;2 2];
        if ~isempty(xv)            plot(xv(:,1),xv(:,2),'*m');
        end
        swarm_obj=handles.swarm;
        b4s=get(swarm_obj,'value')
        switch b4s
            case 1
                Hop.x=[0,0,1.6];
                Hop.angle=[0,0,0];                
            case 2
                Hop.x=[0,0,1.6];
                Hop.angle=[-pi/2.5,0,0];                    
            case 3
                Hop.x=[0,0,1.6];
                Hop.angle=[pi/2.5,0,0];                  
            case 4
                Hop.x=[0,0,1.6];
                Hop.angle=[2*pi/3,0,0];  
            case 5
                Hop.x=[0,0,1];
                Hop.angle=[0,0,pi/2];                 
        end

        numofHop=1;        
end

numberofRobots=size(IDSet,2);
numofvic=size(xv,1);
set(handles.NumofRobots,'string',num2str(numberofRobots));
set(handles.victim,'string',num2str(numofvic));
set(handles.Hop,'string',num2str(numofHop));
G = generateGraph();    % Generate graph representing MRS
% Updating position and orientation for robots from Motion Tracking System with cycle Td

datafolder=pwd;
path=[datafolder,'\Data'];
k='_1';
filename=fullfile(path,['datatime',num2str(numberofRobots),k,'.txt']);
fid_data = fopen(filename,'w'); 
count=0;
datafolder=pwd;
path=[datafolder,'\Data'];
filename=fullfile(path,['victim','.txt']);
fv = fopen(filename,'w'); 
for i=1:size(xv,1)
    fprintf(fv,'%1.4f %1.4f\n',xv(i,1),xv(i,2));
end
fclose(fv);

to=clock;
Td=0.07;
t1=timerfind('Name','timer1');
if isempty(t1)
    t1=timer('Name','timer1','TimerFcn',@PO_RobotUpdate,'Period',Td,'ExecutionMode','fixedSpacing');
else
    stop(t1);
end
start(t1);


% % function addedRobotbyclickmouse(obj,evt) % Insert robots by click mouse
% % global numberofRobots NumofRobotsobj xr
% % global Robot G IDSet maxcreatedRobots
% % 
% % if numberofRobots==maxcreatedRobots
% %     return;
% % end
% % numberofRobots=numberofRobots+1;
% % id=numberofRobots;
% % IDSet=[IDSet,id];
% % pos=get(obj,'CurrentPoint');
% % Robot(id).x=[pos(1,1),pos(1,2)];
% % Robot(id).angle=randn(1)*pi;
% % RobotInit(id,Robot(id).x(1),Robot(id).x(2),Robot(id).angle);
% % G = generateGraph();    % Generate graph representing MRS
% % set(NumofRobotsobj,'string',num2str(numberofRobots));
% % xr=[xr;Robot(id).x];



function PO_RobotUpdate(mTimer,~)
% This function is used to update position and orientation for robots from
% motion tracking system with cycle Td
%global editother_obj idc_obj comobj blobjname
global natnetclient
global hf mode Robot numberofRobots IDSet Hop
global idobjname pxobjname pyobjname orobjname swarm_obj fid_data count 
global editother_obj to td_obj  

% t=tic;

switch mode
    case 1
        C = matlab_client('getsquares');
        for i=1:C.dat(C.n),
            k=8*(i-1);
            id=C.dat(k+1);
            if id<=41
                if id<=40   % for children robot
                    x=C.dat(k+4);    % z
                    y=C.dat(k+2);    % x
                    z=C.dat(k+3);    % y
                    q=[C.dat(k+5),C.dat(k+6),C.dat(k+7),C.dat(k+8)]; %[qx,qy,qz,qw]
                    [c,Ia,Ib]=intersect(id,IDSet);  
                    Robot(Ib).xold=Robot(Ib).x
                    Robot(Ib).x(1)=x;
                    Robot(Ib).x(2)=y;
                    [yaw,roll,pitch]=quat2angle(q,'yxz');%'yxz'
                    Robot(Ib).angle=yaw;
                    %set(Robot(Ib).body,'XData',Robot(Ib).x(1),'YData',Robot(Ib).x(2));
                    %set(Robot(Ib).idtext,'position',[Robot(Ib).x(1)-0.1,Robot(Ib).x(2)-0.1]);
                    RobotUpdate(Ib);
                    Robot(Ib).ds=norm(Robot(Ib).x-Robot(Ib).x);
                else        % for mother robot
                    delta=[0.3675,0.00289];%[0.3313,0.00289]; 
                    y=C.dat(k+2);    % x
                    z=C.dat(k+3);    % y
                    x=C.dat(k+4);    % z
                    q=[C.dat(k+5),C.dat(k+6),C.dat(k+7),C.dat(k+8)]; %[qx,qy,qz,qw]
                    [yaw,roll,pitch]=quat2angle(q,'yxz');%'yxz'
                    Robot(idm).angle=yaw;
                    R=[cos(Robot(idm).angle),sin(Robot(idm).angle); -sin(Robot(idm).angle), cos(Robot(idm).angle)]; % Counterclockwise Rotation Matrix
                    Robot(idm).x=[x,y]+delta*R;
                    MRobotUpdate(idm); % Update robot's body at a new pos.
                    plot(hf,Robot(idm).x(1),Robot(idm).x(2),'.r');
                end
            elseif id>50
                    y=C.dat(k+2);    % x
                    z=C.dat(k+3);    % y
                    x=C.dat(k+4);    % z
                    q=[C.dat(k+5),C.dat(k+6),C.dat(k+7),C.dat(k+8)]; %[qx,qy,qz,qw]
                    [yaw,roll,pitch]=quat2angle(q,'yxz');%'yxz'
                    Hop.x=[x,y,z];
                    Hop.angle=[yaw,roll,pitch];
            end
        end
    case 2
       for i=1:numberofRobots, 
            set(idobjname(i),'string',num2str(i));             % Indicate  ID of robots
            set(pxobjname(i),'string',num2str(Robot(i).x(1),'%1.2f')); % Indicate value x of robot
            set(pyobjname(i),'string',num2str(Robot(i).x(2),'%1.2f')); % Indicate value y of robot
            set(orobjname(i),'string',num2str(Robot(i).angle*180/pi,'%3.2f'));% Indicate orientation value of robot
            RobotUpdate(i); % Update robot's body at a new pos.
            b4s=get(swarm_obj,'value');
            switch b4s
                case 1
                    Hop.x=[0,0,1.6];
                    Hop.angle=[0,0,0];                
                case 2
                    Hop.x=[0,0,1.6];
                    Hop.angle=[-pi/2.5,0,0];                    
                case 3
                    Hop.x=[0,0,1.6];
                    Hop.angle=[pi/2.5,0,0];                  
                case 4
                    Hop.x=[0,0,1.6];
                    Hop.angle=[2*pi/3,0,0]; 
                case 5
                    Hop.x=[0,0,1];
                    Hop.angle=[0,0,pi/2];                     
            end            
       end
end
drawnow;
Tx=etime(clock,to);
set(td_obj,'string',num2str(Tx,'%2.4f'));
to=clock;
% datatime=toc(t);
% count=count+1;
% if ~isempty(fid_data)
%     fprintf(fid_data,'%1.4f\n',datatime);
%     set(editother_obj,'string',num2str(count));
%     if count>=300 
%        fclose(fid_data);
%        fid_data=[];
%        set(editother_obj,'string','Done');
%     end
% end

% --- Executes on button press in Targets.
function Targets_Callback(hObject, eventdata, handles)
% hObject    handle to Targets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global xd FreeTargetIDSet

[FileName, PathName] = uigetfile('*.txt');
Name = fullfile(PathName,FileName);
if PathName==0,
    return; 
end
xr=load(Name);
xd=[];
for i=1:size(xr,1),
    xd=[xd;xr(i,:)];
end

FreeTargetIDSet=[1:1:size(xd,1)];
for i=1: size(xd,1),
    plot(xd(i,1),xd(i,2),'+k');
end

% --- Executes on button press in start_button.
function start_button_Callback(hObject, eventdata, handles)
% hObject    handle to start_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Tc numberofRobots st routeLED IDSet rc epsilon mode vid recordvideoStatus
global step rms dest_index
global sigma   
global LandmarkIDSet OccupiedRobotIDSet FreeTargetIDSet Penalty fid xd cloud group xs
global G Robot pro hf
global status destg detectedgate staticnodes robotset xv bs fireman Backupnodes robotsinroom idn
global numofUVT detectedvictim stime fr  
global Ln r I2 alpha xo xp ts idc
rms=0.1;
step=0;
staticnodes=[];
Backupnodes=[];
st=1;
sigma=0;
count=1;
routeLED=[];
status=0;
pro=get(handles.protable_popmen,'value');   % Selecting a program for starting
datafolder=pwd;
path=[datafolder,'\Data'];
filename=fullfile(path,['Pos',num2str(numberofRobots),'_',num2str(mode),'.txt']);
fid = fopen(filename,'w'); 

filename=fullfile(path,['Route',num2str(numberofRobots),'_',num2str(mode),'.txt']);
fr = fopen(filename,'w'); 

%% record video
recordvideoStatus=get(handles.video,'value');
if recordvideoStatus==1
    namevideo=['MRS',num2str(numberofRobots),'.avi'];
    namevideo=fullfile(path,namevideo);
    vid = VideoWriter(namevideo);  
    vid.Quality = 100;
    vid.FrameRate = 15; %default 15; R5.19-> 9.4; R6.19->9.0
    open(vid);
end

if ~isempty(xp)
    cloud=xp;
    xo=cloud(1,:);
else
    return;
end

switch pro
    case 2      % Swarm Movement: (1) Load Scenario; (2) Robot Init; (3) Load Target; (4) Start
        robotset=[1:1:numberofRobots];
        sigma=0;
        idc=1;
    case 3      % Multi-Target Tracking
        OccupiedRobotIDSet=[];
        LandmarkIDSet=[];
    case 4  % Coverage
        robotset=[1:1:numberofRobots];
        xo=[-1 -0.5];
        dist=[];
        for i=robotset,
            dist=[dist,norm(Robot(i).x-xo)];
        end
        [value,idm]=min(dist);
        bs=robotset(idm);        
        Nbs=[];
        for i=robotset,
           if i~=bs & norm(Robot(bs).x-Robot(i).x)<=0.5
               Nbs=[Nbs,i];               
           end
           Robot(i).xo=Robot(i).x;
        end
        if size(Nbs,2)>=2
            dist=[];
            for i=Nbs,
                dist=[dist,norm(xo-Robot(i).x)];
            end
            [value,index]=sort(dist);
            fireman=[Nbs(index(1)),Nbs(index(2))];
        end  
        fireman=[];
        if ~isempty(fireman)
            for id=robotset,
               if ~any(id==fireman)
                   LED(id,'G','off'); 
                   LED(id,'R','off'); 
                   LED(id,'B','off'); 
               else
                   LED(id,'G','on'); 
               end
               Robot(id).route=[];
               Robot(id).idc=0;
               Robot(id).target=[];
               Robot(id).status=0;           
               Robot(id).xo=Robot(id).x;
            end  
        end

        
        Penalty=[];
        %cloud=[6 -0.5];
        xd=[];
        dist=[];
        robotset=[];
        for id=1:G.NodesN,
            if ~any(id==[bs,fireman])
                dist=[dist,norm(Robot(id).x-cloud)];
                robotset=[robotset,id];
            end
        end
        [value,idm]=sort(dist);
        robotset=robotset(idm);
        id1=robotset(1);
        xd=[xd;Robot(id1).x];
        OccupiedRobotIDSet=[id1];
        idn=1;
        group=[id1];
        Robot(id1).status=2;
        Robot(id1).target=1;
        VTGfull(id1);
        dist=[];
        nt=[];
        for i=2:size(xd,1)
            dist=[dist,norm(cloud-xd(i,:))];
            nt=[nt,i];
        end
        [val,idm]=sort(dist);
        nt=nt(idm);
        FreeTargetIDSet=[nt(1),nt(2)];
        LandmarkIDSet=[id1]; 
        for i=1:size(xd,1),
            if any(i==FreeTargetIDSet)
                plot(hf,xd(i,1),xd(i,2),'.r');  
            else
                plot(hf,xd(i,1),xd(i,2),'.b');  
            end
        end
        numofUVT=size(xd,1)-(size(FreeTargetIDSet,2)+1)
        status=0;
        detectedvictim=[];
        stime=0;
    case 5 % Following Route
        robotset=[1:1:numberofRobots];
        for id=robotset,
           LED(id,'G','off'); 
           LED(id,'R','off'); 
           LED(id,'B','off'); 
           Robot(id).route=[];
           Robot(id).idc=0;
           Robot(id).target=[];
           Robot(id).status=0;
        end         
        %fireman=[IvIDSet(1),IvIDSet(2)];
        fireman=[4 5];
        for id=fireman
            LED(id,'G','on');
            %Robot(id).route=[IvIDSet(3),IvIDSet(4),IvIDSet(5),IvIDSet(6),IvIDSet(7),IvIDSet(8)];
            Robot(id).route=[1,2,3];
            Robot(id).idc=2;
            Robot(id).xdtemp=[];
        end
                
    case 6 % Swarm movement and coverage
        %xo=[-1 -0.5];
        %xo=[0.15 1.5];
        OccupiedRobotIDSet=[];
        status=0;
        %cloud=[6.3 -0.5]; %5 rooms
        %cloud=[2.4 1;2.6 0;3.2 -1.7;5.8 -1.7]; %3 rooms
        dest_index=1;
        destg=[];
        detectedgate=[];
        robotset=[1:1:numberofRobots];
        dist=[];
        for i=robotset,
            dist=[dist,norm(Robot(i).x-xo)];
        end
        [value,idm]=min(dist);
        bs=robotset(idm);        
        Nbs=[];
        for i=robotset,
           if i~=bs & norm(Robot(bs).x-Robot(i).x)<=0.5
               Nbs=[Nbs,i];               
           end
        end
        if size(Nbs,2)>=2
            dist=[];
            for i=Nbs,
                dist=[dist,norm(xo-Robot(i).x)];
            end
            [value,index]=sort(dist);
            fireman=[Nbs(index(1)),Nbs(index(2))];
        end
        staticnodes=[bs];
        Backupnodes=[];
        robotsinroom=[];
        fireman=[];
        for id=robotset,
           if ~any(id==fireman)
               LED(id,'G','off'); 
               LED(id,'R','off'); 
               LED(id,'B','off'); 
           end
           Robot(id).route=[];
           Robot(id).idc=0;
           Robot(id).target=[];
           Robot(id).status=0;
           Robot(id).xo=Robot(id).x;
        end
        xs=[Robot(bs).x];
%         u=cloud(dest_index,:)-Robot(bs).x;
%         xs=Robot(bs).x+(rc-epsilon)*u/norm(u);   

    case 7 % Human operator
        status=0;
        xd=[];
        robotset=[1:1:numberofRobots];
        for id=robotset,
            Robot(id).xo=Robot(id).x;
        end
    case 8 % cyclic pursuit
       robotset=[1,2,3]; 
       status=0;
        
    case 9 % Circular Formation
        n=numberofRobots;
        w=zeros(1,n-2);
        Ln=circulant([1 -1 w],1);
        r=1;
        alpha=[2*pi/2 4*pi/2]';
        %alpha=[2*pi/4 4*pi/4 6*pi/4 0]';
        %alpha=[2*pi/8 4*pi/8 6*pi/8 0]';
        %alpha=[2*pi/3 4*pi/3 0]'; % 3robots
        %alpha=[2*pi/6 4*pi/6 0]';
        %alpha=[2*pi/7 4*pi/7 6*pi/7 8*pi/7 10*pi/7 12*pi/7 0]'; %7 robots
        %alpha=[2*pi/14 4*pi/14 6*pi/14 8*pi/14 10*pi/14 12*pi/14 0]';
        I2=eye(2,2); 
        xo=[];
        for i=1:n,
            xo=[xo;Robot(i).x];
        end
        xo
    case 10 % Concentric formation 
        
        
end

ts=clock;
t3=timerfind('Name','timer3');
if isempty(t3)
    t3=timer('Name','timer3','TimerFcn',@coreprogram,'Period',Tc,'ExecutionMode','fixedSpacing');
else
    stop(t3);
end
start(t3);

function coreprogram(mTimer,~)
%% Core program for Multi-Robot System
% global G step editrunstep_obj pro st
global G Robot editrunstep_obj step pro vid recordvideoStatus
global OccupiedRobotIDSet LandmarkIDSet  FreeTargetIDSet 
global xd hf Penalty fid 
global status rms robotsinroom robotoutroom  idt xv bs idv fireman Backupnodes ido
global detectedgate staticnodes robotset GateP group cloud dest_index st idn ad rh rc epsilon xs
global firemanstatus stime Istaticnodes
global numofUVT detectedvictim 
global Hop numofvic idg g fr xc numberofRobots mode L ra
global Ln r I2 alpha Tc vmax x3 xo ts tp_obj idc xp
    step=step+1
    set(editrunstep_obj,'string',step);
    G = generateGraph();
    firemanonwork=[4 5];
    if ~isempty(fireman)& ~any(status==firemanonwork)
        for i=1:G.NodesN,
            Robot(i).N=setdiff(Robot(i).N,fireman);
            Robot(i).Ne=setdiff(Robot(i).Ne,fireman);
            Robot(i).Na=setdiff(Robot(i).Na,fireman);            
        end
    end
    
    switch pro
        case 2  % Swarm movement
            xdk=xp(idc,:);
            ReachTarget=SM(xdk,robotset);
            if ReachTarget==1
                if idc<size(xp,1),
                    idc=idc+1;
                else
                    idc=1;
                end
            end
        case 3  % Multi-Target Tracking
            MTT();
            for id=1:G.NodesN,
               dist=[];
               for i=1:size(Robot(id).N,2),
                   j=Robot(id).N(i);
                   dist=[dist,norm(Robot(id).x-Robot(j).x)];
               end
               [valuei,idmin]=min(dist);
               j=Robot(id).N(idmin);
            end
            if size(OccupiedRobotIDSet,2)==G.NodesN,
                fclose(fid);
                id=100;
                CommandSend(id,0,0,0,0); % Stop all Robots
                t3=timerfind('Name','timer3');
                if ~isempty(t3)
                    stop(t3);% stop timer 3 for robotcontrol
                end
            end
        case 4  % Coverage  
            status
            switch status
                case 0
                    a=idn;
                    b=min(size(robotset,2),size(xd,1)-numofUVT); % robotset does not includes [bs,fireman]
                    if a<b
                        for i=a+1:b,
                                id=robotset(i);
                                group=[group,id];
                                idn=idn+1;
                        end
                    end
                    group
                    LandmarkIDSet
                    FreeTargetIDSet
                    OccupiedRobotIDSet
                    xd
                    ng=size(group,2)
                    nb=size(robotset,2)
                    nt=size(xd,1)
                    size(OccupiedRobotIDSet,2)
                    numofUVT

                    
                    detectedvictim
                   for i=group,
                       if Robot(i).status==0 & isempty(intersect(Robot(i).N,OccupiedRobotIDSet)),
                           xdk=xd(1,:);
                           HDC(i,xdk);
                       end
                   end
                   if idn<size(robotset,2)-numofUVT,
                       for i=1:numofUVT
                           id=robotset(idn+i);
                           if isempty(intersect(Robot(id).N,OccupiedRobotIDSet))
                               va=[0 0];
                               Neighbours=intersect(Robot(id).N,group);
                               if ~isempty(Neighbours)
                                   for j=Neighbours,
                                       va=va+Robot(j).va;
                                   end
                               end
                               if norm(va)>0
                                   xdk=Robot(id).x+va;
                                   HDC(id,xdk);
                               end
                           else
                               RobotStop(id); 
                           end
                       end
                   end
                   DCA(group);
                   if ((size(OccupiedRobotIDSet,2)+1==G.NodesN)|(size(OccupiedRobotIDSet,2)==size(xd,1)-numofUVT))                    
                        % firemans move toward an identified victim
                        %% identify victims
                        bs
                        route=Robot(bs).route
                        if isempty(Robot(bs).route)
                            c=1;
                            for i=group,                       
                                if ~isempty(Robot(i).VictiminSen),
                                    c=0;
                                    if isempty(detectedvictim)
                                        victim=Robot(i).VictiminSen;
                                    else
                                        victim=setdiff(Robot(i).VictiminSen,detectedvictim);
                                    end
                                    if ~isempty(victim)
                                        idv=i;
                                        detectedvictim=[detectedvictim,victim];    
                                        [cdist,route]=shortestpathFindSMC(bs,idv,[OccupiedRobotIDSet,idv]);
                                        Robot(bs).route=route;
                                        Robot(bs).idc=1;
                                        for j=fireman,
                                            Robot(j).route=route; 
                                            Robot(j).idc=2;
                                        end
                                        
                                        for k=1:10
                                            for i=1:size(Robot(bs).route,2),
                                               id=Robot(bs).route(i);
                                               LED(id,'R','on');
                                            end
                                            for i=fireman,
                                                LED(i,'G','on');
                                            end  
                                        end   
                                        firemanstatus=0;
                                        if ~isempty(Robot(bs).route)
                                            for k=1:size(Robot(bs).route,2),
                                                fprintf(fr,'%d ',Robot(bs).route(k));
                                            end
                                            fprintf(fr,'\n');
                                        end                                        
                                        break;
                                    end
                                end
                            end
                            if c==1,
                                for id=1:G.NodesN,
                                    RobotStop(id); % Stop all Robots  
                                end
                                 
                                idn=size(OccupiedRobotIDSet,2);
                                group=[];
                                status=1;
                                for i=OccupiedRobotIDSet,
                                    Robot(i).route=[];
                                    Robot(i).idc=1;
                                    Robot(i).status=0;
                                end                                 
                            end
                            
                        else        %~isempty(Robot(bs).route)
                            firemanstatus    
                            switch firemanstatus
                                case 0 
                                   idc=Robot(fireman(1)).idc;
                                   id_dest=Robot(fireman(1)).route(idc);
                                   xt=Robot(id_dest).x;
                                   [idr,ReachTarget]=SM4fireman(xt,fireman);
                                   if ReachTarget==1
                                       if Robot(fireman(1)).idc<size(Robot(fireman(1)).route,2)
                                           id=setdiff(fireman,idr);
                                           if norm(Robot(id).x-xt)>0.4
                                               HDC(id,xt);
                                           else
                                               Robot(fireman(1)).idc=Robot(fireman(1)).idc+1;
                                           end                                               
                                       else
                                           id=setdiff(fireman,idr);
                                           if norm(Robot(id).x-xt)>0.4
                                               HDC(id,xt);
                                           else
                                                firemanstatus=1;
                                           end
                                       end               
                                   end
                                case 1
                                    idc=Robot(fireman(1)).idc;
                                    if  Robot(fireman(1)).route(idc)==idv                                    
                                        xdk=xv(Robot(idv).VictiminSen(1),:);
                                        s=0;
                                        if norm(Robot(fireman(1)).x-xdk)>0.2
                                            HDC(fireman(1),xdk);
                                        else
                                            RobotStop(fireman(1));
                                            s=1;
                                        end
                                        
                                        s1=0;
                                        if norm(Robot(fireman(2)).x-xdk)>0.2
                                            HDC(fireman(2),xdk);
                                        else
                                            RobotStop(fireman(2));
                                            s1=1;
                                        end
                                        
                                        if (s==1)|(s1==1)
                                            stime=0;
                                            firemanstatus=2;
                                        end
                                    end
                                    if Robot(fireman(1)).route(idc)==bs,
                                        s=0;
                                        if norm(Robot(fireman(1)).x-Robot(fireman(1)).xo)>rms
                                            HDC(fireman(1),Robot(fireman(1)).xo);
                                        else
                                            RobotStop(fireman(1));
                                            s=1;
                                        end
                                        
                                        s1=0;
                                        if norm(Robot(fireman(2)).x-Robot(fireman(2)).xo)>rms
                                            HDC(fireman(2),Robot(fireman(2)).xo);
                                        else
                                            RobotStop(fireman(2));
                                            s1=1;
                                        end
                                        
                                        if s==1&s1==1
                                            for k=1:10,
                                                for i=Robot(bs).route,
                                                   LED(i,'R','off');
                                                end
                                                for id=fireman,
                                                    LED(id,'G','off'); 
                                                end
                                            end 
                                            stime=stime+1;
                                            if stime>10
                                                stime=0;
                                                firemanstatus=0;
                                                Robot(bs).route=[];
                                                Robot(fireman(1)).route=[];
                                                Robot(fireman(2)).route=[];
                                                Robot(fireman(1)).idc=1;
                                                Robot(fireman(2)).idc=1;
                                                if size(detectedvictim,2)==numofvic,
                                                    for id=1:G.NodesN,
                                                        RobotStop(id); % Stop all Robots   
                                                    end
                                                    idn=size(OccupiedRobotIDSet,2);
                                                    group=[];
                                                    status=1;
                                                    for i=OccupiedRobotIDSet,
                                                        Robot(i).route=[];
                                                        Robot(i).idc=1;
                                                        Robot(i).status=0;
                                                    end                                           
                                                end   
                                            end
                                        end
                                    end
                                case 2
                                    stime=stime+1;
                                    for id=fireman,
                                        RobotStop(id);
                                    end
                                    if stime>10
                                        firemanstatus=0;
                                        index=[size(Robot(bs).route,2):-1:1];
                                        Robot(fireman(1)).route=Robot(bs).route(index);
                                        Robot(fireman(1)).idc=1;
                                        Robot(fireman(2)).route=Robot(bs).route(index);
                                        Robot(fireman(2)).idc=1;
                                        stime=0;
                                    end                                    
                            end
                        end
                   end
%                     %% exploration is done
%                     if isempty(Robot(bs).route)&((size(OccupiedRobotIDSet,2)+3==G.NodesN)|(size(OccupiedRobotIDSet,2)==size(xd,1)-numofUVT))
%                         RobotStop(100); % Stop all Robots   
%                         stime=stime+1;
%                         if stime>20
%                             idn=size(OccupiedRobotIDSet,2);
%                             group=[];
%                             status=1;
%                             for i=OccupiedRobotIDSet,
%                                 Robot(i).route=[];
%                                 Robot(i).idc=1;
%                                 Robot(i).status=0;
%                             end
%                         end
%                     end
                case 1 % Thu hoi ve
                        if mod(step,30)==0
                            id=OccupiedRobotIDSet(idn);
                            group=[group,id];
                            [cdist,route]=shortestpathFindSMC(id,OccupiedRobotIDSet(1),setdiff(OccupiedRobotIDSet,group));
                            Robot(id).route=route;
                            Robot(id).idc=1;
                            if idn>2
                                idn=idn-1;
                            end
                        end
                        c=0;
                        for id=group,
                            switch Robot(id).status
                                case 0
                                    reach=InNetworkMovement(id);
                                    if reach==1
                                        Robot(id).status=1;
                                    end
                                case 1
                                    xdk=Robot(id).xo;
                                    if norm(Robot(id).x-Robot(id).xo)>rms
                                        HDC(id,xdk);
                                    else
                                        RobotStop(id);
                                        c=c+1;
                                    end
                            end
                        end
                        
                        if c==size(OccupiedRobotIDSet,2)-1,
                            t3=timerfind('Name','timer3');
                            if ~isempty(t3)
                                stop(t3);% stop timer 3 for robotcontrol
                            end 
                        end
            end
            
        case 5  % following a route
           switch mod(step,10)
               case 0
                   for i=Robot(fireman(1)).route,
                       LED(i,'R','on');
                       buzzer(i);
                   end
                   for i=Robot(fireman(2)).route,
                       LED(i,'R','on');
                       buzzer(i);
                   end                   
                       
               case 4
                   for i=Robot(fireman(1)).route,
                       LED(i,'R','off');
                   end
                   for i=Robot(fireman(2)).route,
                       LED(i,'R','off');
                   end                     
           end
           switch status
               case 0
                   idc=Robot(fireman(1)).idc;
                   id_dest=Robot(fireman(1)).route(idc);
                   xt=Robot(id_dest).x;
                   [idr,ReachTarget]=SM4fireman(xt,fireman);
                   if ReachTarget==1
                       if Robot(fireman(1)).idc<size(Robot(fireman(1)).route,2)
                           id=setdiff(fireman,idr); 
                           if norm(Robot(id).x-xt)>0.4
                               HDC(id,xt);
                           else
                               Robot(fireman(1)).idc=Robot(fireman(1)).idc+1;
                           end                            
                       else
                           id=setdiff(fireman,idr);
                           if norm(Robot(id).x-xt)>0.4
                               HDC(id,xt);
                           else
                                status=1;
                           end
                       end               
                   end
               case 1
                   for i=Robot(fireman(1)).route,
                       LED(i,'R','off');
                   end
                   for i=Robot(fireman(2)).route,
                       LED(i,'R','off');
                   end
                   
                   for i=fireman,
                       LED(i,'G','off');
                   end
                   for id=1:G.NodesN
                        RobotStop(id); % Stop all Robots
                   end
                   t3=timerfind('Name','timer3');
                   if ~isempty(t3)
                       stop(t3);% stop timer 3 for robotcontrol
                   end
           end
        case 6 % Swarm movement and coverage
            
            if status<=7           
                group=setdiff(robotset,[staticnodes,fireman,Backupnodes]);    % Group of activate robots for swarm movement
            end
            if any(status==[0 3])
                % establish a backbone on corridor
                ids=staticnodes(size(staticnodes,2));
                Nc=[];
                for i=group
                    if status==0
                        if ~isempty(Robot(i).Nc)& any(ids==Robot(i).Nc)
                            Nc=[Nc,i]; 
                        end
                    else
                        if ~isempty(Robot(i).Nc)& any(ids==Robot(i).Nc)& i~=ido
                            Nc=[Nc,i]; 
                        end
                    end
                end
                if ~isempty(Nc)
                    dist=[];
                    for i=Nc,
                        dist=[dist,norm(Robot(ids).x-Robot(i).x)];
                    end
                    [dm,idm]=min(dist);
                    staticnodes=[staticnodes,Nc(idm)];
                    group=setdiff(group,Nc(idm));
                    s1=xs(size(staticnodes,2)-1,:);
                    u=cloud(dest_index,:)-s1;
                    s2=s1+(rc-epsilon)*u/norm(u);
                    xs=[xs;s2];                    
                end            
            end 
            if size(staticnodes,2)>1 & status<7
                for i=2:size(staticnodes,2),
                    xdk=xs(i,:);
                    id=staticnodes(i);
                    if norm(Robot(id).x-xdk)>rms
                        HDC(id,xdk);
                    else
                        RobotStop(id)   %stop robot id
                    end
                end            
            end            
            
           switch status
                case 0  % swarm movement
                    reachtarget=SM(cloud(dest_index,:),group)               % Swarm movement                      
                    if reachtarget==1 
                        if dest_index<size(cloud,1)
                            dest_index=dest_index+1;
                        else
                            idt=size(staticnodes,2);
                            Istaticnodes=[];
                            status=7;
                        end
                    end                    
                    
                    for id=group,                           % Gate Detection
                       s=0;
                       if ~isempty(Robot(id).gate),
                           if isempty(detectedgate),
                               detectedgate=[detectedgate;Robot(id).gate(1,:)];
                               s=1;                               
                           else
                               for j=1:size(Robot(id).gate,1),
                                   if isempty(intersect(Robot(id).gate(j,:),detectedgate,'rows'))
                                       detectedgate=[detectedgate;Robot(id).gate(j,:)];
                                        s=1;
                                       break;
                                   end
                               end
                           end
                       end
                       if s==1
                            for id=1:G.NodesN,
                                RobotStop(id); % Stop all Robots
                            end
                            
                            idg=size(detectedgate,1);
                            A=[detectedgate(idg,1), detectedgate(idg,2)];
                            B=[detectedgate(idg,3), detectedgate(idg,4)];
                            C=[detectedgate(idg,5), detectedgate(idg,6)];
                             
                            I=(A+B)/2;
                            AB=B-A;
                            IC=C-I;
                            n=[-AB(1,2),AB(1,1)];
                            P1=I+0.45*rh*n/norm(n);%0.5*rh*n/norm(n);
                            IP1=P1-I;
                            if dot(IC,IP1)/(norm(IC)*norm(IP1))>0   % P1 va C cung phia voi AB (P1 trong phong)
                                P2=I-0.5*rh*n/norm(n);
                                GateP=[I;P1;P2];                    % GateP=[Giua, trong, ngoai]    
                            else                            
                                P1=I-0.45*rh*n/norm(n);
                                P2=I+0.5*rh*n/norm(n);
                                GateP=[I;P1;P2];
                            end
                            plot(hf, GateP(:,1),GateP(:,2),'or');
                            status=1;
                            break;
                       end
                    end
                case 1  % swarm moves to the detected gate
                    SM(GateP(3,:),group);       % Swarm movement to gate
                    s=0;
                    idg=size(detectedgate,1);
                    for id=group,      
                        dist=norm(Robot(id).x-GateP(3,:));%distP2Seg(Robot(id).x,detectedgate(idg,:));
                        if dist<=0.3 %~isempty(dist)&
                            RobotStop(100); % Stop all Robots
                            s=1;
                            break;
                        end
                    end
                    
                    if s==1   % Enough near Gate
                        RobotStop(100); % Stop all Robots
                        dist=[];
                        for id=group,
                            dist=[dist,norm(Robot(id).x-GateP(1,:))];
                        end
                        [value,Index]=sort(dist);
                        robotoutroom=group(Index)
                        if size(robotoutroom,2)>=2
                            status=2;
                        end
                    end
                case 2  % First robot is deployed into room                        
                    id1=robotoutroom(1);
                    id2=robotoutroom(2);
                    subgroup=setdiff(group,[id1,id2]);
                    for i=subgroup,
                        RobotStop(i); % Stop robot i
                    end
                    xd1=GateP(2,:); % Gate-Point inside room
                    xd2=GateP(3,:); % Gate_Point Outside room;

                    if norm(Robot(id2).x-xd2)>rms
                        HDC(id2,xd2);
                    else
                        RobotStop(id2); % Stop robot id2
                    end
                    
                    if norm(Robot(id1).x-xd1)>rms
                        HDC(id1,xd1);
                    else
                        RobotStop(id1); % Stop robot id1
                    end
                    if norm(Robot(id2).x-xd2)<=rms & norm(Robot(id1).x-xd1)<=rms
                        xd=xd1;
                        Penalty=[];
                        FreeTargetIDSet=[];
                        Robot(id1).status=2;
                        Robot(id1).target=1;  
                        nv=GateP(2,:)-GateP(1,:);
                        u=[1,0];
                        phi_i0=acos(dot(nv,u)/norm(nv));
                        idg=size(detectedgate,1);
                        VTG(id1,phi_i0,detectedgate(idg,:));
                        plot(hf,xd(:,1),xd(:,2),'.b');
                        robotsinroom=[id1,id2];
                        robotoutroom=setdiff(robotoutroom,[id1,id2]);
                        LandmarkIDSet=id1;
                        OccupiedRobotIDSet=id1;
                        status=3;
                        idn=[];
                        ad=0;
                    end
                    
                case 3 % deploy into the dectected room
                   idg=size(detectedgate,1);
                   DCAinRoom(robotsinroom,detectedgate(idg,:)); 
                   robotoutroom=setdiff(robotoutroom,staticnodes);
                   
                   s=1;
                   if isempty(robotoutroom)
                       id1=robotsinroom(size(robotsinroom,2));
                       if ad==0
                           idn=addingRobot(id1);
                           if idn~=0    % khong co robot de add
                               ad=1;
                           end
                       end
                       if ad==1
                           s=0;
                           reach1=InNetworkMovement(idn);
                           if reach1==1,
                               s=1;
                               robotoutroom=[robotoutroom,idn]; 
                               ad=0;
                           end                       
                       end
                   else
                       dist=[];                       
                       for id=robotoutroom,
                           dist=[dist,norm(Robot(id).x-GateP(3,:))];
                       end
                       [value,Index]=sort(dist);
                       robotoutroom=robotoutroom(Index);
                       id1=robotoutroom(1);
                       xd1=GateP(3,:);
                       if norm(Robot(id1).x-xd1)>rms
                           HDC(id1,xd1);
                       else
                           RobotStop(id1); % Stop robot id1
                       end
                        
                       Nc=intersect(Robot(id1).Nc,robotoutroom);
                       if ~isempty(Nc)  % Lam the nao de dung khi ko con la critical nodes?
                            for i=Nc,
                               HDC(i,Robot(id1).x); 
                            end
                       else
                           Ni=setdiff(Robot(id1).N,[robotsinroom,staticnodes,Backupnodes,fireman]);
                           if ~isempty(Ni)
                               for i=Ni,
                                  RobotStop(i); % Stop robot i 
                               end
                           end
                       end                       
                       
                       Ni=setdiff(Robot(id1).N,[robotsinroom,staticnodes,Backupnodes,fireman]);
                       if ~isempty(Ni)
                          dist=[];
                          for i=Ni,
                              dist=[dist,G.R(id1,i)];                              
                          end                        
                          [value,idm]=min(dist);
                          j=Ni(idm);
                          if value>0.3
                              HDC(j,Robot(id1).x);
                          else
                              RobotStop(j);; % Stop robot j                           
                          end
                       end
                                              
                       if (norm(Robot(id1).x-xd1)<=0.1 ) & (size(robotsinroom,2) <size(xd,1))
                           if size(robotoutroom,2)>1
                               robotsinroom=[robotsinroom,id1]
                           end
                           robotoutroom=setdiff(robotoutroom,id1)
                       end
                       ido=id1;                       
                   end
                    if ((size(OccupiedRobotIDSet,2)==size(robotsinroom,2)) & s==1)|(size(robotsinroom,2)==2 & size(xd,1)==1)
                        RobotStop(100);; % Stop all Robots
                        for id=robotsinroom,
                            Robot(id).status=0;
                            Robot(id).target=[];
                            Robot(id).route=[];
                            Robot(id).xdtemp=[];
                        end
                        if size(xd,1)==1,
                            ido=robotsinroom(2);
                            robotsinroom=robotsinroom(1);
                        end
                        LandmarkIDSet=[ido];
                        FreeTargetIDSet=[1:1:size(xd,1)];
                        OccupiedRobotIDSet=[ido,robotsinroom];
                        Robot(ido).status=2;
                        idg=size(detectedgate,1);
                        A=[detectedgate(idg,1), detectedgate(idg,2)];
                        B=[detectedgate(idg,3), detectedgate(idg,4)];
                        AB=B-A;
                        Po=GateP(3,:);
                        Pi=GateP(2,:);
                        n=Po-Pi;
                        xd=[];
                        k=0;
                        for i=0:5,
                            if mod(i,2)==0
                                row=[-2 -1 0 1];
                            else
                                row=[1 0 -1 -2];
                            end                            
                            
                            for j=row,
                                P=Po+i*0.3*n/norm(n)+j*0.3*AB/norm(AB);
                                if k<size(robotsinroom,2)
                                    robotoutroom=setdiff(group,robotsinroom);
                                    check=1;
                                    for m=robotoutroom
                                        if norm(P-Robot(m).x)<=0.2,
                                            check=0;
                                        end
                                    end
                                    if check==1,
                                        xd=[xd;P];
                                        k=k+1;
                                    end
                                end
                            end
                        end
                        plot(hf,xd(:,1),xd(:,2),'*b');
                        robotoutroom=[ido,robotsinroom(size(robotsinroom,2))]; %Tap robot thu hien di chuyen ra ngoai room, robot id1 vai tro lam cot moc dan duong.
                        idt=size(xd,1);
                        
                        idv=0; 
                        for id=robotsinroom,
                            if ~isempty(Robot(id).VictiminSen)
                                idv=id;
                                break;
                            end
                        end
%                         victim=[];
%                         subrobot=[];
%                         for id=robotsinroom,
%                             if ~isempty(Robot(id).VictiminSen)
%                                 subrobot=[subrobot,id];
%                                 VictiminSeni=setdiff(Robot(id).VictiminSen,victim)
%                                 if ~isempty(VictiminSeni)
%                                     victim=[victim,VictiminSeni];
%                                 end
%                             end
%                         end
%                         dist=[];
%                         for id=subrobot
%                             dist=[dist,norm(Robot(id).x-xv(victim(1),:))];
%                         end
%                         [val,idm]=min(dist);
%                         idv=subrobot(idm);



                        if idv==0   % dont have any victim
                            status=6;
                            stime=0;
                        else
                            [cdist1,route1]=shortestpathFindSMC(bs,ido,robotset)%setdiff(robotset,fireman));
                            [cdist2,route2]=shortestpathFindSMC(ido,idv,robotset)%setdiff(robotset,fireman));

                            if ~isempty(fireman)
                                Robot(fireman(1)).route=route1;
                                Robot(fireman(1)).idc=2;
                                index=[2:1:size(route2,2)];
                                Robot(fireman(2)).route=route2(index);
                                Robot(fireman(2)).idc=1;
                                Robot(bs).route=[Robot(fireman(1)).route,Robot(fireman(2)).route];
                                Robot(bs).idc=1;  
                                for k=1:10
                                    for i=1:size(Robot(bs).route,2),
                                       id=Robot(bs).route(i);
                                       LED(id,'R','on');
                                    end
                                    for i=fireman,
                                        LED(i,'G','on');
                                    end  
                                end
                                if ~isempty(Robot(bs).route)
                                    for k=1:size(Robot(bs).route,2),
                                        fprintf(fr,'%d ',Robot(bs).route(k));
                                    end
                                    fprintf(fr,'\n');
                                end                                
                                status=4;
                                firemanstatus=0;
                                stime=0;
                            else
                                stime=0;
                                status=6;
                            end
                        end
                        %st=0;
                    end
                case 4  % Can viet mot ham chuyen dong cho 2 fireman; di song hanh voi nhau va bam theo route
                    
%                    switch mod(step,10)
%                        case 0
%                            for i=Robot(bs).route
%                                LED(i,'R','on');
%                                buzzer(i);
%                            end                 
%                            for i=fireman,
%                                LED(i,'G','on');
%                            end
%                        case 4
%                            for i=Robot(bs).route
%                                LED(i,'R','off');
%                            end
%                    end
                    
                    switch firemanstatus
                        case 0 
                           idc=Robot(fireman(1)).idc;
                           id_dest=Robot(fireman(1)).route(idc);
                           xt=Robot(id_dest).x;
                           [idr,ReachTarget]=SM4fireman(xt,fireman);
                           if ReachTarget==1
                               if Robot(fireman(1)).idc<size(Robot(fireman(1)).route,2)
                                    Robot(fireman(1)).idc=Robot(fireman(1)).idc+1;
                               else
                                   id=setdiff(fireman,idr);
                                   if norm(Robot(id).x-xt)>0.35
                                       HDC(id,xt);
                                   else
                                        firemanstatus=1;
                                   end
                               end               
                           end
                        case 1
                           route=Robot(fireman(2)).route
                           idc=Robot(fireman(2)).route(1);
                           xdk=Robot(idc).x;
                           id1=fireman(1);
                           id2=fireman(2);
                           s1=0;
                           s2=0;
                           d1=norm(Robot(id1).x-xdk);
                           d2=norm(Robot(id2).x-xdk);
                           if d1>d2
                               r1=0.4;
                               r2=0.3;
                           else
                               r1=0.3;
                               r2=0.4;
                           end
                           if d1>r1
                               HDC(id1,xdk);
                           else
                               RobotStop(id1);
                               s1=1;
                           end
                           if d2>r2
                               HDC(id2,xdk);
                           else
                               RobotStop(id2);
                               s2=1;
                           end                   
                           if s1==1 & s2==1
                               firemanstatus=2;
                           end                            
                        case 2
                           idc=Robot(fireman(2)).idc;
                           id_dest=Robot(fireman(2)).route(idc);
                           xt=Robot(id_dest).x;
                           [idr,ReachTarget]=SM4fireman(xt,fireman);
                           if ReachTarget==1
                               if Robot(fireman(2)).idc<size(Robot(fireman(2)).route,2)
                                    Robot(fireman(2)).idc=Robot(fireman(2)).idc+1;
                               else
                                   firemanstatus=3;
                               end
                           end
                        case 3
                            xdk=xv(Robot(idv).VictiminSen(1),:);
                            %va=xdk-Robot(idv).x;
                            s=0;
                            id1=fireman(1);
                            if norm(Robot(id1).x-xdk)>0.2
                                HDC(id1,xdk);
                            else
                                RobotStop(id1);
                                s=1;
                            end
                            
                            s1=0;
                            id2=fireman(2);
                            if norm(Robot(id2).x-xdk)>0.2
                                HDC(id2,xdk);
                            else
                                RobotStop(id2);
                                s1=1;
                            end
                                
                            if (s==1)&(s1==1)
                                firemanstatus=4;
                            end
                        case 4
                            for id=fireman,
                                RobotStop(id);
                            end
                            stime=stime+1;
                            id=fireman(1);
                            if stime>10
                                status=5;
                                firemanstatus=0;
                                index1=[size(Robot(fireman(1)).route,2):-1:1];
                                route1=Robot(fireman(1)).route(index1);                                
                                index2=[size(Robot(fireman(2)).route,2):-1:1];
                                route2=Robot(fireman(2)).route(index2);                                
                                Robot(fireman(1)).route=route2;
                                Robot(fireman(2)).route=route1;
                                Robot(fireman(1)).idc=1;
                                Robot(fireman(2)).idc=1;
                                stime=0;
                            end
                    end
                case 5                 
%                    switch mod(step,10)
%                        case 0
%                            for i=Robot(bs).route,
%                                LED(i,'R','on');
%                                buzzer(i);
%                            end               
%                            for i=fireman,
%                                LED(i,'G','on');
%                            end
%                        case 4
%                            for i=Robot(bs).route
%                                LED(i,'R','off');
%                            end
%                    end  
                    
                    switch firemanstatus
                        case 0 
                           idc=Robot(fireman(1)).idc;
                           id_dest=Robot(fireman(1)).route(idc);
                           xt=Robot(id_dest).x;
                           [idr,ReachTarget]=SM4fireman(xt,fireman);
                           if ReachTarget==1
                               if Robot(fireman(1)).idc<size(Robot(fireman(1)).route,2)-1
                                    Robot(fireman(1)).idc=Robot(fireman(1)).idc+1;
                               else
                                   id=setdiff(fireman,idr);
                                   if norm(Robot(id).x-xt)>0.35
                                       HDC(id,xt);
                                   else
                                        firemanstatus=1;
                                   end
                               end               
                           end
                        case 1
                           idc=Robot(fireman(2)).route(1);
                           xdk=Robot(idc).x;
                           id1=fireman(1);
                           id2=fireman(2);
                           s1=0;
                           s2=0;
                           d1=norm(Robot(id1).x-xdk);
                           d2=norm(Robot(id2).x-xdk);
                           if d1>d2
                               r1=0.4;
                               r2=0.3;
                           else
                               r1=0.3;
                               r2=0.4;
                           end
                           if d1>r1
                               HDC(id1,xdk);
                           else
                               RobotStop(id1);
                               s1=1;
                           end
                           if d2>r2
                               HDC(id2,xdk);
                           else
                               RobotStop(id2);
                               s2=1;
                           end                   
                           if s1==1 & s2==1
                               firemanstatus=2;
                           end 
                        case 2
                           idc=Robot(fireman(2)).idc;
                           id_dest=Robot(fireman(2)).route(idc);
                           xt=Robot(id_dest).x;
                           [idr,ReachTarget]=SM4fireman(xt,fireman);
                           if ReachTarget==1
                               if Robot(fireman(2)).idc<size(Robot(fireman(2)).route,2)
                                    Robot(fireman(2)).idc=Robot(fireman(2)).idc+1;
                               else
                                   firemanstatus=3;
                               end
                           end
                        case 3
                            s=0;
                            id1=fireman(1);
                            if norm(Robot(id1).x-Robot(id1).xo)>rms
                                HDC(id1,Robot(id1).xo);
                            else 
                                RobotStop(id1);
                                s=1;
                            end
                            
                            s1=0;
                            id2=fireman(2);
                            if norm(Robot(id2).x-Robot(id2).xo)>rms
                                HDC(id2,Robot(id2).xo);
                            else
                                RobotStop(id2);
                                s1=1;
                            end
                                
                            if s==1&s1==1
                                firemanstatus=4;
                            end                            
                            
                        case 4
                           for id=fireman
                                RobotStop(id);
                                Robot(id).route=[];
                                Robot(id).idc=[];
                           end
                           for k=1:10,
                                for i=Robot(bs).route,
                                   LED(i,'R','off');
                                end
                                for id=fireman,
                                    LED(id,'G','off'); 
                                end
                            end                            
                            Robot(bs).route=[];
                            stime=55;
                            firemanstatus=0;
                            status=6;                             
                    end
                case 6
                    stime=stime+1;
                    if stime>20,
                        id=robotoutroom(size(robotoutroom,2));
                        MTTC(robotoutroom); 
                        idg=size(detectedgate,1);
                        idc=size(robotsinroom,2)-(size(robotoutroom,2)-1);
                        if distP2Seg(Robot(id).x,detectedgate(idg,:))<=0.4 & idc>0                       
                            robotoutroom=[robotoutroom,robotsinroom(idc)];
                        end
                        if size(OccupiedRobotIDSet,2)-size(robotsinroom,2)==size(xd,1)+1
                            RobotStop(100); % Stop all Robots
                            for id=robotoutroom,
                                Robot(id).status=0;
                                Robot(id).target=[];
                                Robot(id).route=[];
                                Robot(id).xdtemp=[];
                            end
                           status=0;
                        end
                    end
                case 7 % thu hoi robot
                    if ~isempty(Istaticnodes)
                        group=[group,Istaticnodes];
                    end
                    xdk=Robot(staticnodes(idt)).x;
                    reachtarget=SM(xdk,group);
                    if reachtarget
                        Istaticnodes=[Istaticnodes,staticnodes(idt)];
                        if idt>2,
                            idt=idt-1;
                        else
                            group=[group,staticnodes(idt)];
                            xd=[];
                            for i=0:7
                                for j=2:-1:-2,
                                    if size(xd,1)<size(group,2),
                                        %P=[0.15+0.3*i,-0.4+0.3*j];
                                        P=[xo(1,1)+0.3*i,xo(1,2)+0.3*j];
                                        check=1;
                                        for k=[bs,fireman],
                                            if norm(P-Robot(k).x)<=0.2
                                                check=0;
                                            end
                                        end
                                        if check==1,
                                            xd=[xd;P];
                                        end
                                    end
                                end
                            end
                            dist=[];
                            for id=group,
                                dist=[dist,norm(Robot(bs).x-Robot(id).x)];
                            end
                            [val,idm]=sort(dist);
                            group=group(idm);
                            for i=1:size(group,2),
                                id=group(i);
                                Robot(id).xo=xd(i,:);
                            end
                            status=8;
                        end
                    end
                case 8
                    c=0;
                    for id=group,
                        xdk=Robot(id).xo;
                        if norm (Robot(id).x-xdk)>rms,
                            HDC(id,xdk);
                        else
                            c=c+1;
                            RobotStop(id);
                        end
                    end
                    if c==size(group,2),
                        st=0;
                    end
                    
           end 
            
        case 7 % Human Operator
            % Nhan lenh khi hanh dong da thuc hien xong
%             id=1;
%             xd
%             idt
            status
%             xx=Hop.x(3)>0.8
%             aa1=norm(Hop.angle(1))>pi/2
            switch status
                case 0  
                    if Hop.x(3)>0.8  
                        if (size(xd,1)==0) 
                            if norm(Hop.angle(1))<=pi/6 % Forward (swarm2F)? 
                                xd=[1.5,-0.5;6,-0.5];
                                idt=2;
                                status=1;
                            else
                                stime=0;
                                if (Hop.angle(1)<-pi/4)& (Hop.angle(1)>=-pi/2)    % Turn right (swarm2LR)
                                    idg=5;
                                    idt=2;
                                    xd=[1.5,-0.5;5.8,-0.5];                                    
                                    status=2;
                                end
                                if (Hop.angle(1)> pi/4)&(Hop.angle(1)<= pi/2) % Turn left (swarm2RR)
                                    idt=2;
                                    idg=3;
                                    xd=[1.5,-0.5;5,-0.5];                                    
                                    status=2;
                                end
                                A=[g(idg,1), g(idg,2)];
                                B=[g(idg,3), g(idg,4)];
                                C=[g(idg,5), g(idg,6)];

                                I=(A+B)/2;
                                AB=B-A;
                                IC=C-I;
                                n=[-AB(1,2),AB(1,1)];
                                P1=I+0.45*rh*n/norm(n);
                                IP1=P1-I;
                                if dot(IC,IP1)/(norm(IC)*norm(IP1))>0   % P1 va C cung phia voi AB (P1 trong phong)
                                    P2=I-0.5*rh*n/norm(n);
                                    GateP=[I;P1;P2];                    % GateP=[Giua, trong, ngoai]    
                                else                            
                                    P1=I-0.45*rh*n/norm(n);
                                    P2=I+0.5*rh*n/norm(n);
                                    GateP=[I;P1;P2];
                                end
                                plot(hf, GateP(:,1),GateP(:,2),'or');                                
                            end
                        end
                        if norm(Hop.angle(1))> pi/2                   % Home (swarm2H)
                            if size(xd,1)==0 
                                xt=[1,-0.5];
                                s=1;
                                for i=robotset,
                                    if norm(Robot(i).x-xt)<0.5
                                        s=0;
                                    end
                                end
                                if s==1
                                    xd=xt;
                                    idt=1;
                                end
                            else
                                idt=size(xd,1);
                                status=4;
                            end
                            stime=0;
                        end
                    end
                case 1  % task 
                    % Swarm Movement
                    xdk=xd(idt,:);
                    reachtarget=SM(xdk,robotset) 
                    if reachtarget==1
                        if idt<size(xd,1)
                            idt=idt+1;
                        else
                            for i=1:5,
                                RobotStop(100);
                            end
                            status=0;
                        end
                    end
                case 2
                    switch stime
                        case 0
                            % Swarm Movement
                            xdk=xd(idt,:);
                            reachtarget=SM(xdk,robotset) 
                            if reachtarget==1
                                if idt<size(xd,1)
                                    idt=idt+1;
                                else
                                    for i=1:5,
                                        RobotStop(100);
                                    end
                                    stime=1;
                                end
                            end                            
                        case 1  % swarm moves to the detected gate
                            SM(GateP(3,:),robotset);       % Swarm movement to gate
                            s=0;
                            for id=robotset,      
                                dist=norm(Robot(id).x-GateP(3,:)); 
                                if dist<=0.3  
                                    RobotStop(100); % Stop all Robots
                                    s=1;
                                    break;
                                end
                            end                            
                            if s==1   % Enough near Gate
                                RobotStop(100); % Stop all Robots
                                dist=[];
                                for id=robotset,
                                    dist=[dist,norm(Robot(id).x-GateP(1,:))];
                                end
                                [value,Index]=sort(dist);
                                robotoutroom=robotset(Index);
                                stime=2;
                            end                            
                        case 2 % First robot is deployed into room  
                            id1=robotoutroom(1);
                            id2=robotoutroom(2);
                            subgroup=setdiff(robotset,[id1,id2]);
                            for i=subgroup,
                                RobotStop(i); % Stop robot i
                            end
                            xd1=GateP(2,:); % Gate-Point inside room
                            xd2=GateP(3,:); % Gate_Point Outside room;

                            if norm(Robot(id2).x-xd2)>rms
                                HDC(id2,xd2);
                            else
                                RobotStop(id2); % Stop robot id2
                            end

                            if norm(Robot(id1).x-xd1)>rms
                                HDC(id1,xd1);
                            else
                                RobotStop(id1); % Stop robot id1
                            end
                            if norm(Robot(id2).x-xd2)<=rms & norm(Robot(id1).x-xd1)<=rms
                                xd=xd1;
                                Penalty=[];
                                FreeTargetIDSet=[];
                                Robot(id1).status=2;
                                Robot(id1).target=1;  
                                nv=GateP(2,:)-GateP(1,:);
                                u=[1,0];
                                phi_i0=acos(dot(nv,u)/norm(nv));
                                VTG(id1,phi_i0,g(idg,:));
                                plot(hf,xd(:,1),xd(:,2),'.b');
                                robotsinroom=[id1,id2];
                                robotoutroom=setdiff(robotoutroom,[id1,id2]);
                                LandmarkIDSet=id1;
                                OccupiedRobotIDSet=id1;
                                stime=3;
                            end 
                        case 3      % deployment                      
                           DCAinRoom(robotsinroom,g(idg,:)); 
                           s=1;
                           if ~isempty(robotoutroom)
                               dist=[];                       
                               for id=robotoutroom,
                                   dist=[dist,norm(Robot(id).x-GateP(3,:))];
                               end
                               [value,Index]=sort(dist);
                               robotoutroom=robotoutroom(Index);
                               id1=robotoutroom(1);
                               xd1=GateP(3,:);
                               if norm(Robot(id1).x-xd1)>rms
                                   HDC(id1,xd1);
                               else
                                   RobotStop(id1); % Stop robot id1
                               end

                               Nc=intersect(Robot(id1).Nc,robotoutroom);
                               if ~isempty(Nc)  % Lam the nao de dung khi ko con la critical nodes?
                                    for i=Nc,
                                       HDC(i,Robot(id1).x); 
                                    end
                               else
                                   Ni=setdiff(Robot(id1).N,robotsinroom);
                                   if ~isempty(Ni)
                                       for i=Ni,
                                          RobotStop(i); % Stop robot i 
                                       end
                                   end
                               end                       

                               Ni=setdiff(Robot(id1).N,robotsinroom);
                               if ~isempty(Ni)
                                  dist=[];
                                  for i=Ni,
                                      dist=[dist,G.R(id1,i)];                              
                                  end                        
                                  [value,idm]=min(dist);
                                  j=Ni(idm);
                                  if value>0.3
                                      HDC(j,Robot(id1).x);
                                  else
                                      RobotStop(j);; % Stop robot j                           
                                  end
                               end

                               if (norm(Robot(id1).x-xd1)<=0.1 ) & (size(robotsinroom,2) <size(xd,1))
                                   if size(robotoutroom,2)>1
                                       robotsinroom=[robotsinroom,id1]
                                   end
                                   robotoutroom=setdiff(robotoutroom,id1)
                               end
                               ido=id1;                       
                           end
                            if ((size(OccupiedRobotIDSet,2)==size(robotsinroom,2)) & s==1)|(size(robotsinroom,2)==2 & size(xd,1)==1)
                                RobotStop(100);; % Stop all Robots
                                for id=robotsinroom,
                                    Robot(id).status=0;
                                    Robot(id).target=[];
                                    Robot(id).route=[];
                                    Robot(id).xdtemp=[];
                                end
                                if size(xd,1)==1,
                                    ido=robotsinroom(2);
                                    robotsinroom=robotsinroom(1);

                                end
                                LandmarkIDSet=[ido];
                                FreeTargetIDSet=[1:1:size(xd,1)];
                                OccupiedRobotIDSet=[ido,robotsinroom];
                                Robot(ido).status=2;
                                A=[g(idg,1), g(idg,2)];
                                B=[g(idg,3), g(idg,4)];
                                AB=B-A;
                                Po=GateP(3,:);
                                Pi=GateP(2,:);
                                n=Po-Pi;
                                xd=[];
                                k=0;
                                for i=0:5,
                                    if mod(i,2)==0
                                        row=[-2 -1 0 1];
                                    else
                                        row=[1 0 -1 -2];
                                    end                            

                                    for j=row,
                                        P=Po+i*0.3*n/norm(n)+j*0.3*AB/norm(AB);
                                        if k<size(robotsinroom,2)
                                            robotoutroom=setdiff(robotset,robotsinroom);
                                            check=1;
                                            for m=robotoutroom
                                                if norm(P-Robot(m).x)<=0.2,
                                                    check=0;
                                                end
                                            end
                                            if check==1,
                                                xd=[xd;P];
                                                k=k+1;
                                            end
                                        end
                                    end
                                end
                                plot(hf,xd(:,1),xd(:,2),'*b');
                                robotoutroom=[ido,robotsinroom(size(robotsinroom,2))]; %Tap robot thu hien di chuyen ra ngoai room, robot id1 vai tro lam cot moc dan duong.
                                idt=size(xd,1);
                                idn=0;
                                stime=4;
                            end 
                        case 4 % withdraw 
                            idn=idn+1;
                            if idn<5
                               for i=robotsinroom,
                                   LED(i,'R','on');
                                   buzzer(i);
                               end  
                            elseif idn<10
                                for i=robotsinroom,
                                    LED(i,'R','off');
                                    buzzer(i);
                                end  
                            end
                            if idn>10
                                id=robotoutroom(size(robotoutroom,2));
                                MTTC(robotoutroom); 
                                idc=size(robotsinroom,2)-(size(robotoutroom,2)-1);
                                if distP2Seg(Robot(id).x,g(idg,:))<=0.4 & idc>0                       
                                    robotoutroom=[robotoutroom,robotsinroom(idc)];
                                end
                                if size(OccupiedRobotIDSet,2)-size(robotsinroom,2)==size(xd,1)+1
                                    RobotStop(100); % Stop all Robots
                                    for id=robotoutroom,
                                        Robot(id).status=0;
                                        Robot(id).target=[];
                                        Robot(id).route=[];
                                        Robot(id).xdtemp=[];
                                    end
                                    xd=[1,-0.5]; 
                                    stime=0;
                                    status=0;
                                end
                            end
                    end
                case 4  % Home
                    % Swarm movement
                    % Out room
                    switch stime
                        case 0
                            xdk=xd(idt,:);
                            reachtarget=SM(xdk,robotset) 
                            if reachtarget==1,
                               if idt>1
                                   idt=idt-1;
                               else
                                   stime=1;
                                   xd=[];
                                   for i=0:6
                                       for j=[0,1,2,-1], 
                                           if size(xd,1)<size(robotset,2),
                                               x=[0.15+0.25*i,-0.5+0.25*j];
                                               xd=[xd;x];
                                           end
                                       end
                                   end
                                   plot(hf,xd(:,1),xd(:,2),'.r');
                                   DestID=[1:1:size(xd,1)]; %dao nguov DestID va robo
                                   rs=robotset;
                                   for i=DestID,
                                       dist=[];
                                       for id=rs,
                                           dist=[dist,norm(xd(i,:)-Robot(id).x)];
                                       end
                                       [val,idx]=min(dist);    
                                       idm=rs(idx);
                                       Robot(idm).xo=xd(i,:);
                                       xo=Robot(idm).xo
                                       rs=setdiff(rs,idm);
                                   end
                               end                        
                            end
                        case 1
                            c=0;
                            for id=robotset,
                                xdk=Robot(id).xo;
                                if norm (Robot(id).x-xdk)>rms,
                                    HDC(id,xdk);
                                else
                                    c=c+1;
                                    RobotStop(id);
                                end
                            end  
                            
                            if c==size(robotset,2)
                                for i=1:5
                                    RobotStop(100);
                                end
                                xd=[];
                                status=0;
                            end
                    end
                    
%                    xdk=xd(idt,:); 
%                    if norm(Robot(id).x-xdk)>0.2
%                        HDC(id,xdk);   
%                    else
%                        if idt>1
%                            idt=idt-1;
%                        else
%                            RobotStop(id);
%                            xd=[];
%                            status=0;
%                        end
%                    end
            end
        case 8 % cyclic pursuit
            
            %n=numberofRobots;
            if mod(step,100)==0&size(robotset,2)<numberofRobots,
                nr=size(robotset,2);
                robotset=[robotset,nr+1];
            end
            cstr=['r','g','b','m','y','k','r','g','b','m'];
            for aa=1:size(robotset,2),
                if aa==size(robotset,2)
                    bb=1;
                else
                    bb=aa+1;
                end
                i=robotset(aa);
                j=robotset(bb);
                posi=Robot(i).x;
                posj=Robot(j).x;
                anglei=Robot(i).angle;
                anglej=Robot(j).angle;
                Robot(i).alpha=cal_alpha(i,j,posi,posj,anglei,anglej); % Related Orientation
                vs=[0 0];
                sigma=0;
                s=1;
                if ~isempty(Robot(i).Na)                       % Separation vector for collision avoidance
                    for a=Robot(i).Na,
                        rij=Robot(a).x-Robot(i).x;
                        cosphi=dot(rij,Robot(i).ni)/norm(rij);
                        wij=sigma+(1-sigma)*(1+cosphi)/2;
                        vs=vs-wij*rij*exp(-25*(norm(rij)-L/2-ra)); %wij*rij*exp(-25*(norm(rij)-ra)); 
                    end
                    dist=norm(vs);
                    idset=[i];
                    for a=Robot(i).Na,
                        if norm(Robot(a).vs),
                            dist=[dist,norm(Robot(a).vs)];
                            idset=[idset,a];
                        end
                    end
                    if norm(vs)==max(dist)
                        s=0;
                    end
                end                 
                Robot(i).vs=vs;                
                if s==0
                      vel=0;
                end        
                
                switch mode
                    case 1
                        [dleft,dright,vleft,vright]=linear_vel_2_angular_vel_of_ws(i);
                        CommandSend(i,dleft,dright,vleft,vright);
                    case 2
                        [wl,wr]=linear_vel_2_angular_vel_of_ws_simu(i)
                        WMR(i,[wl;wr]);
                end
                plot(hf,Robot(i).x(1,1),Robot(i).x(1,2),'.','color',cstr(i));
            end
        case 9  % Circular formation
            cstr=['r','g','b','m','y','k','r','g','b','m'];
            n=numberofRobots;
            %vmax=0.8;
            k1=vmax/2*ones(n,1);  % k1 in range (0 vmax)
            K=0.4;  % K in range (0, v/r) 1.25
            h=0.05;
            c=zeros(2*n,1);
            c1=zeros(2,1);
            jj=1;            
            x3
            r=1;
            for i=1:n,
                x_alpha=x3'-alpha;
                angle=Ln(i,:)*x_alpha;
                u1(i)=vmax-(k1(i)*tanh(angle));
                for j=1:n
                    c1=[Robot(j).x(1)-r*sin(Robot(j).angle), Robot(j).x(2)+r*cos(Robot(j).angle)];
                    c(jj:jj+1,1)=c1;
                    jj=jj+2;
                end 
                
                jj=1;
                L2=kron(Ln(i,:),I2);
                y(i)=[cos(Robot(i).angle) sin(Robot(i).angle)]*L2*c;
                phi=1/(1+norm(y(i))); %1.5
                u2(i)=(u1(i)/r)+K*phi*y(i);
                x3(i)=x3(i)+h*u2(i);                
% % %                 Robot(i).x(1)=Robot(i).x(1)+(h*u1(i))*cos(Robot(i).angle);
% % %                 Robot(i).x(2)=Robot(i).x(2)+(h*u1(i))*sin(Robot(i).angle);
% % %                 Robot(i).angle=Robot(i).angle+(h*u2(i));  
% %                 x=xo(i,:)+[(h*u1(i))*cos(x3(i)),(h*u1(i))*sin(x3(i))];
% %                 xo(i,:)=x;
% %                 x3(i)=x3(i)+h*u2(i);
% %                 Robot(i).v=(x-Robot(i).x);
% %                 if norm(Robot(i).v)<0.4
% %                     Robot(i).v=0.4*Robot(i).v/norm(Robot(i).v);
% %                 end
% %                 vs=[0 0];
% %                 sigma=0;
% %                 s=1;
% %                 if ~isempty(Robot(i).Na)                       % Separation vector for collision avoidance
% %                     for a=Robot(i).Na,
% %                         rij=Robot(a).x-Robot(i).x;
% %                         cosphi=dot(rij,Robot(i).ni)/norm(rij);
% %                         wij=sigma+(1-sigma)*(1+cosphi)/2;
% %                         vs=vs-wij*rij*exp(-25*(norm(rij)-L/2-ra)); %wij*rij*exp(-25*(norm(rij)-ra)); 
% %                     end
% %                     dist=norm(vs);
% %                     idset=[i];
% %                     for a=Robot(i).Na,
% %                         if norm(Robot(a).vs),
% %                             dist=[dist,norm(Robot(a).vs)];
% %                             idset=[idset,a];
% %                         end
% %                     end
% %                     if norm(vs)==max(dist)
% %                         s=0;
% %                     end
% %                 end                 
% %                 Robot(i).vs=vs;                
% %                 if s==0
% %                       Robot(i).v=[0,0];
% %                 end                
% %                 L2AVel_Circular1(i);
% %                 plot(hf,xo(i,1),xo(i,2),'.','color',cstr(i));                  

                vel=norm([(h*u1(i))*cos(Robot(i).angle),(h*u1(i))*sin(Robot(i).angle)])/Tc
                Robot(i).alpha=h*u2(i);
                vs=[0 0];
                sigma=0;
                s=1;
                if ~isempty(Robot(i).Na)                       % Separation vector for collision avoidance
                    for a=Robot(i).Na,
                        rij=Robot(a).x-Robot(i).x;
                        cosphi=dot(rij,Robot(i).ni)/norm(rij);
                        wij=sigma+(1-sigma)*(1+cosphi)/2;
                        vs=vs-wij*rij*exp(-25*(norm(rij)-L/2-ra)); %wij*rij*exp(-25*(norm(rij)-ra)); 
                    end
                    dist=norm(vs);
                    idset=[i];
                    for a=Robot(i).Na,
                        if norm(Robot(a).vs),
                            dist=[dist,norm(Robot(a).vs)];
                            idset=[idset,a];
                        end
                    end
                    if norm(vs)==max(dist)
                        s=0;
                    end
                end                 
                Robot(i).vs=vs;                
                if s==0
                      vel=0;
                end
                
                L2AVel_Circular(i,vel);
 
                plot(hf,Robot(i).x(1,1),Robot(i).x(1,2),'.','color',cstr(i));
            end
        case 10 % Concentric formation
            
    end
    for id=1:G.NodesN,
        fprintf(fid,'%d %1.4f %1.4f %1.4f %d %d\n',id,Robot(id).x(1),Robot(id).x(2),Robot(id).angle, Robot(id).dLED(1), Robot(id).dLED(2));
    end
    
    if recordvideoStatus==1 
        vid.writeVideo(getframe(hf));
    end  
    
    Tp=etime(clock,ts);
    set(tp_obj,'string',num2str(Tp,'%2.4f'));
    ts=clock;
    if st==0
        fclose(fid);
        fclose(fr);
        RobotStop(100); % Stop all Robots
        t3=timerfind('Name','timer3');
        if ~isempty(t3)
            stop(t3);% stop timer 3 for robotcontrol
        end
    end

function RobotTestPro(mTimer,~)
global editrunstep_obj IDSet comobj
global G Robot rms xd dest_index step   


step=step+1;
set(editrunstep_obj,'string',step);
G = generateGraph();
id=1;
xdi=xd(dest_index,:);
HDC(id,xdi);

if norm(Robot(id).x-xdi)<rms
    IDl=11;
    status=0;
    frame=['L',IDSet(id),IDl,status];
    fwrite(comobj,frame);    % Green LED on
    frame=['S',IDSet(id)];     
    fwrite(comobj,frame);    % Buzzer sounds out
    
    if dest_index<size(xd,1)
        dest_index=dest_index+1;
    else
        dest_index=1;
%         frame=['R',id,'0','0','0','0']; 
%         fwrite(comobj,frame);            % stop command is sent to robot
%         stop(timerfind('Name','timer3'));% stop timer 3 for robotcontrol
    end
    IDl=11;
    status=1;
    frame=['L',IDSet(id),IDl,status];
    fwrite(comobj,frame);    % Green LED off    
end



function SwarmMovementPro(mTimer,~)
global editrunstep_obj comobj
global G Robot rms xd dest_index step
global sigma IDSet

step=step+1;
set(editrunstep_obj,'string',step);
G = generateGraph();

for id=1:G.NodesN,
   Robot(id).va=xd(dest_index,:)-Robot(id).x; % robot is aware of common direction for movement
   synV=[0 0];
   ni=[cos(Robot(id).angle),sin(Robot(id).angle)];
   for i=Robot(id).N,                           % Orientation Synchronization with neighbours
       rij=Robot(i).x-Robot(id).x;
       cosphi=dot(rij,ni)/norm(rij);
       wij=sigma+(1-sigma)*(1-cosphi)/2;
       synV=synV+wij*Robot(i).v;        
   end
   synV=synV/norm(size(Robot(id).N,2));
   xdk=Robot(id).x+Robot(id).va+synV;
   HDC(id,xdk); 
   
   if norm(Robot(id).x-xd(dest_index,:))<rms
        IDl=11;
        status=0;
        frame=['L',IDSet(id),IDl,status];
        fwrite(comobj,frame);    % Green LED on
        frame=['S',IDSet(id)];     
        fwrite(comobj,frame);    % Buzzer sounds out
       
       if dest_index<size(xd,1)
           dest_index=dest_index+1;
       else
           dest_index=1;
           %         id=100;
           %         frame=['R',id,'0','0','0','0']; 
           %         fwrite(comobj,frame);            % stop command is sent to robot
           %         stop(timerfind('Name','timer3'));% stop timer 3 for robotcontrol
       end 
       IDl=11;
       status=1;
       frame=['L',IDSet(id),IDl,status];
       fwrite(comobj,frame);    % Green LED off    
       break;
   end
end

function Multi_target_Tracking(mTimer,~)
%%Multi-Target Tracking
global step editrunstep_obj OccupiedRobotIDSet 
global G Robot LandmarkIDSet FreeTargetIDSet fid 
    step=step+1
    OccupiedRobotIDSet
    FreeTargetIDSet
    LandmarkIDSet
    set(editrunstep_obj,'string',step);
    G = generateGraph();
    MTT();
    for id=1:G.NodesN,
       dist=[];
       for i=1:size(Robot(id).N,2),
           j=Robot(id).N(i);
           dist=[dist,norm(Robot(id).x-Robot(j).x)];
       end
       [valuei,idmin]=min(dist);
       j=Robot(id).N(idmin);
       fprintf(fid,'%d %1.4f %1.4f %d %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f\n',id,Robot(id).x(1),Robot(id).x(2),j,Robot(j).x(1),Robot(j).x(2),Robot(id).va(1), Robot(id).va(2),Robot(id).vs(1),Robot(id).vs(2),Robot(id).v(1),Robot(id).v(2));
    end
    if size(OccupiedRobotIDSet,2)==G.NodesN,
        fclose(fid);
        for id=1:G.NodesN,
            CommandSend(id,0,0,0,0); % Stop Robot
        end
        t5=timerfind('Name','timer5');
        if ~isempty(t5)
            stop(timerfind('Name','timer5'));% stop timer 5 for target tracking
        end   
    end


% --- Executes on button press in stop_button.
function stop_button_Callback(hObject, eventdata, handles)
% hObject    handle to stop_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global numberofRobots fid mode Robot fid_data fr recordvideoStatus vid

if recordvideoStatus==1    
    close(vid);    
end  

if ~isempty(fid)
    fclose(fid);
    fid=[];
end

if ~isempty(fr)
    fclose(fr);
    fr=[];
end

if ~isempty(fid_data),
    fclose(fid_data);
    fid_data=[];
end

if mode==1
    CommandSend(100,0,0,0,0); % stop all robots
end

if numberofRobots~=0
    for id=1:numberofRobots,
        LED(id,'G','off'); 
        LED(id,'R','off'); 
        LED(id,'B','off'); 
        Robot(id).v=[0 0];
        Robot(id).va=[0 0];
        Robot(id).vc=[0 0];
        Robot(id).vs=[0 0];
    end
end

t3=timerfind('Name','timer3');
if ~isempty(t3)
    stop(timerfind('Name','timer3'));% stop timer 3 for main program
end

t4=timerfind('Name','timer4');
if ~isempty(t4)
    stop(timerfind('Name','timer4'));% stop timer 4 for swarm movement
end

t5=timerfind('Name','timer5');
if ~isempty(t5)
    stop(timerfind('Name','timer5'));% stop timer 4 for swarm movement
end




% --- Executes on button press in fw_button.
function fw_button_Callback(hObject, eventdata, handles)
% hObject    handle to fw_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
IDstr=get(handles.idc,'string');
vleftstr=get(handles.vleft,'string');
vrightstr=get(handles.vright,'string');
ID=str2num(IDstr);
vleft=str2num(vleftstr);
vright=str2num(vrightstr);
dleft=0;    % forward
dright=0;   % forward
IvIDSet(ID)
CommandSend(IvIDSet(ID),dleft,dright,vleft,vright);

% id=45;
% cmd=1;
% v1=0;
% v2=0;
% v3=4;
% frame=['m',id,cmd,v1,v2,v3]
% fwrite(comobj,frame);



% id=45;
% cmd=3;
% v1=5;%pick between 0 - 5
% v2=0;
% v3=0;
% frame=['m',id,cmd,v1,v2,v3]
% fwrite(comobj,frame); 




% --- Executes on button press in tl_button.
function tl_button_Callback(hObject, eventdata, handles)
% hObject    handle to tl_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

IDstr=get(handles.idc,'string');
vleftstr=get(handles.vleft,'string');
vrightstr=get(handles.vright,'string');
ID=str2num(IDstr);
vleft=str2num(vleftstr);
vright=str2num(vrightstr);
dleft=1;    % backward
dright=0;   % forward
CommandSend(IvIDSet(ID),dleft,dright,vleft,vright);

% --- Executes on button press in tr_button.
function tr_button_Callback(hObject, eventdata, handles)
% hObject    handle to tr_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

IDstr=get(handles.idc,'string');
vleftstr=get(handles.vleft,'string');
vrightstr=get(handles.vright,'string');
ID=str2num(IDstr);
vleft=str2num(vleftstr);
vright=str2num(vrightstr);
dleft=0;    % forward
dright=1;   % backward
CommandSend(IvIDSet(ID),dleft,dright,vleft,vright);

% --- Executes on button press in bw_button.
function bw_button_Callback(hObject, eventdata, handles)
% hObject    handle to bw_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

IDstr=get(handles.idc,'string');
vleftstr=get(handles.vleft,'string');
vrightstr=get(handles.vright,'string');
ID=str2num(IDstr);
vleft=str2num(vleftstr);
vright=str2num(vrightstr);
dleft=1;    % backward
dright=1;   % backward
CommandSend(IvIDSet(ID),dleft,dright,vleft,vright);

 

% --- Executes on button press in cl_pushbutton.
function cl_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to cl_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% stop timer 1 (update robot's pos)
t1=timerfind('Name','timer1');
if ~isempty(t1)
    stop(t1);
    delete(t1);
end

% stop timer 3 (Robot Test)
t3=timerfind('Name','timer3');
if ~isempty(t3)
    stop(t3);
    delete(t3);
end

% stop timer 4 for robot tracking points
t4=timerfind('Name','timer4');
if ~isempty(t4)
    stop(t4);
    delete(t4);
end

% stop timer 5 for homing
t5=timerfind('Name','timer5');
if ~isempty(t5)
    stop(t5);
    delete(t5);
end


close all


% --- Outputs from this function are returned to the command line.
function varargout = MRS_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function id1_Callback(hObject, eventdata, handles)
% hObject    handle to id1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id1 as text
%        str2double(get(hObject,'String')) returns contents of id1 as a double


% --- Executes during object creation, after setting all properties.
function id1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function id2_Callback(hObject, eventdata, handles)
% hObject    handle to id2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id2 as text
%        str2double(get(hObject,'String')) returns contents of id2 as a double


% --- Executes during object creation, after setting all properties.
function id2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function id3_Callback(hObject, eventdata, handles)
% hObject    handle to id3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id3 as text
%        str2double(get(hObject,'String')) returns contents of id3 as a double


% --- Executes during object creation, after setting all properties.
function id3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function id4_Callback(hObject, eventdata, handles)
% hObject    handle to id4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id4 as text
%        str2double(get(hObject,'String')) returns contents of id4 as a double


% --- Executes during object creation, after setting all properties.
function id4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function id5_Callback(hObject, eventdata, handles)
% hObject    handle to id5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id5 as text
%        str2double(get(hObject,'String')) returns contents of id5 as a double


% --- Executes during object creation, after setting all properties.
function id5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function id6_Callback(hObject, eventdata, handles)
% hObject    handle to id6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id6 as text
%        str2double(get(hObject,'String')) returns contents of id6 as a double


% --- Executes during object creation, after setting all properties.
function id6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px1_Callback(hObject, eventdata, handles)
% hObject    handle to px1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px1 as text
%        str2double(get(hObject,'String')) returns contents of px1 as a double


% --- Executes during object creation, after setting all properties.
function px1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px2_Callback(hObject, eventdata, handles)
% hObject    handle to px2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px2 as text
%        str2double(get(hObject,'String')) returns contents of px2 as a double


% --- Executes during object creation, after setting all properties.
function px2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px3_Callback(hObject, eventdata, handles)
% hObject    handle to px3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px3 as text
%        str2double(get(hObject,'String')) returns contents of px3 as a double


% --- Executes during object creation, after setting all properties.
function px3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px4_Callback(hObject, eventdata, handles)
% hObject    handle to px4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px4 as text
%        str2double(get(hObject,'String')) returns contents of px4 as a double


% --- Executes during object creation, after setting all properties.
function px4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px5_Callback(hObject, eventdata, handles)
% hObject    handle to px5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px5 as text
%        str2double(get(hObject,'String')) returns contents of px5 as a double


% --- Executes during object creation, after setting all properties.
function px5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px6_Callback(hObject, eventdata, handles)
% hObject    handle to px6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px6 as text
%        str2double(get(hObject,'String')) returns contents of px6 as a double


% --- Executes during object creation, after setting all properties.
function px6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py1_Callback(hObject, eventdata, handles)
% hObject    handle to py1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py1 as text
%        str2double(get(hObject,'String')) returns contents of py1 as a double


% --- Executes during object creation, after setting all properties.
function py1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py2_Callback(hObject, eventdata, handles)
% hObject    handle to py2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py2 as text
%        str2double(get(hObject,'String')) returns contents of py2 as a double


% --- Executes during object creation, after setting all properties.
function py2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py3_Callback(hObject, eventdata, handles)
% hObject    handle to py3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py3 as text
%        str2double(get(hObject,'String')) returns contents of py3 as a double


% --- Executes during object creation, after setting all properties.
function py3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py4_Callback(hObject, eventdata, handles)
% hObject    handle to py4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py4 as text
%        str2double(get(hObject,'String')) returns contents of py4 as a double


% --- Executes during object creation, after setting all properties.
function py4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py5_Callback(hObject, eventdata, handles)
% hObject    handle to py5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py5 as text
%        str2double(get(hObject,'String')) returns contents of py5 as a double


% --- Executes during object creation, after setting all properties.
function py5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py6_Callback(hObject, eventdata, handles)
% hObject    handle to py6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py6 as text
%        str2double(get(hObject,'String')) returns contents of py6 as a double


% --- Executes during object creation, after setting all properties.
function py6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or1_Callback(hObject, eventdata, handles)
% hObject    handle to or1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or1 as text
%        str2double(get(hObject,'String')) returns contents of or1 as a double


% --- Executes during object creation, after setting all properties.
function or1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or2_Callback(hObject, eventdata, handles)
% hObject    handle to or2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or2 as text
%        str2double(get(hObject,'String')) returns contents of or2 as a double


% --- Executes during object creation, after setting all properties.
function or2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or3_Callback(hObject, eventdata, handles)
% hObject    handle to or3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or3 as text
%        str2double(get(hObject,'String')) returns contents of or3 as a double


% --- Executes during object creation, after setting all properties.
function or3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or4_Callback(hObject, eventdata, handles)
% hObject    handle to or4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or4 as text
%        str2double(get(hObject,'String')) returns contents of or4 as a double


% --- Executes during object creation, after setting all properties.
function or4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or5_Callback(hObject, eventdata, handles)
% hObject    handle to or5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or5 as text
%        str2double(get(hObject,'String')) returns contents of or5 as a double


% --- Executes during object creation, after setting all properties.
function or5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or6_Callback(hObject, eventdata, handles)
% hObject    handle to or6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or6 as text
%        str2double(get(hObject,'String')) returns contents of or6 as a double


% --- Executes during object creation, after setting all properties.
function or6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl1_Callback(hObject, eventdata, handles)
% hObject    handle to bl1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl1 as text
%        str2double(get(hObject,'String')) returns contents of bl1 as a double


% --- Executes during object creation, after setting all properties.
function bl1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl2_Callback(hObject, eventdata, handles)
% hObject    handle to bl2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl2 as text
%        str2double(get(hObject,'String')) returns contents of bl2 as a double


% --- Executes during object creation, after setting all properties.
function bl2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl3_Callback(hObject, eventdata, handles)
% hObject    handle to bl3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl3 as text
%        str2double(get(hObject,'String')) returns contents of bl3 as a double


% --- Executes during object creation, after setting all properties.
function bl3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl4_Callback(hObject, eventdata, handles)
% hObject    handle to bl4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl4 as text
%        str2double(get(hObject,'String')) returns contents of bl4 as a double


% --- Executes during object creation, after setting all properties.
function bl4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl5_Callback(hObject, eventdata, handles)
% hObject    handle to bl5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl5 as text
%        str2double(get(hObject,'String')) returns contents of bl5 as a double


% --- Executes during object creation, after setting all properties.
function bl5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl6_Callback(hObject, eventdata, handles)
% hObject    handle to bl6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl6 as text
%        str2double(get(hObject,'String')) returns contents of bl6 as a double


% --- Executes during object creation, after setting all properties.
function bl6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function id7_Callback(hObject, eventdata, handles)
% hObject    handle to id7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id7 as text
%        str2double(get(hObject,'String')) returns contents of id7 as a double


% --- Executes during object creation, after setting all properties.
function id7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px7_Callback(hObject, eventdata, handles)
% hObject    handle to px7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px7 as text
%        str2double(get(hObject,'String')) returns contents of px7 as a double


% --- Executes during object creation, after setting all properties.
function px7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py7_Callback(hObject, eventdata, handles)
% hObject    handle to py7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py7 as text
%        str2double(get(hObject,'String')) returns contents of py7 as a double


% --- Executes during object creation, after setting all properties.
function py7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or7_Callback(hObject, eventdata, handles)
% hObject    handle to or7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or7 as text
%        str2double(get(hObject,'String')) returns contents of or7 as a double


% --- Executes during object creation, after setting all properties.
function or7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl7_Callback(hObject, eventdata, handles)
% hObject    handle to bl7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl7 as text
%        str2double(get(hObject,'String')) returns contents of bl7 as a double


% --- Executes during object creation, after setting all properties.
function bl7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function id8_Callback(hObject, eventdata, handles)
% hObject    handle to id8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id8 as text
%        str2double(get(hObject,'String')) returns contents of id8 as a double


% --- Executes during object creation, after setting all properties.
function id8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px8_Callback(hObject, eventdata, handles)
% hObject    handle to px8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px8 as text
%        str2double(get(hObject,'String')) returns contents of px8 as a double


% --- Executes during object creation, after setting all properties.
function px8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py8_Callback(hObject, eventdata, handles)
% hObject    handle to py8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py8 as text
%        str2double(get(hObject,'String')) returns contents of py8 as a double


% --- Executes during object creation, after setting all properties.
function py8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or8_Callback(hObject, eventdata, handles)
% hObject    handle to or8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or8 as text
%        str2double(get(hObject,'String')) returns contents of or8 as a double


% --- Executes during object creation, after setting all properties.
function or8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl8_Callback(hObject, eventdata, handles)
% hObject    handle to bl8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl8 as text
%        str2double(get(hObject,'String')) returns contents of bl8 as a double


% --- Executes during object creation, after setting all properties.
function bl8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function id9_Callback(hObject, eventdata, handles)
% hObject    handle to id9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id9 as text
%        str2double(get(hObject,'String')) returns contents of id9 as a double


% --- Executes during object creation, after setting all properties.
function id9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px9_Callback(hObject, eventdata, handles)
% hObject    handle to px9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px9 as text
%        str2double(get(hObject,'String')) returns contents of px9 as a double


% --- Executes during object creation, after setting all properties.
function px9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py9_Callback(hObject, eventdata, handles)
% hObject    handle to py9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py9 as text
%        str2double(get(hObject,'String')) returns contents of py9 as a double


% --- Executes during object creation, after setting all properties.
function py9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or9_Callback(hObject, eventdata, handles)
% hObject    handle to or9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or9 as text
%        str2double(get(hObject,'String')) returns contents of or9 as a double


% --- Executes during object creation, after setting all properties.
function or9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl9_Callback(hObject, eventdata, handles)
% hObject    handle to bl9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl9 as text
%        str2double(get(hObject,'String')) returns contents of bl9 as a double


% --- Executes during object creation, after setting all properties.
function bl9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function id10_Callback(hObject, eventdata, handles)
% hObject    handle to id10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id10 as text
%        str2double(get(hObject,'String')) returns contents of id10 as a double


% --- Executes during object creation, after setting all properties.
function id10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px10_Callback(hObject, eventdata, handles)
% hObject    handle to px10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px10 as text
%        str2double(get(hObject,'String')) returns contents of px10 as a double


% --- Executes during object creation, after setting all properties.
function px10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py10_Callback(hObject, eventdata, handles)
% hObject    handle to py10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py10 as text
%        str2double(get(hObject,'String')) returns contents of py10 as a double


% --- Executes during object creation, after setting all properties.
function py10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or10_Callback(hObject, eventdata, handles)
% hObject    handle to or10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or10 as text
%        str2double(get(hObject,'String')) returns contents of or10 as a double


% --- Executes during object creation, after setting all properties.
function or10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl10_Callback(hObject, eventdata, handles)
% hObject    handle to bl10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl10 as text
%        str2double(get(hObject,'String')) returns contents of bl10 as a double


% --- Executes during object creation, after setting all properties.
function bl10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editmes_Callback(hObject, eventdata, handles)
% hObject    handle to editmes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editmes as text
%        str2double(get(hObject,'String')) returns contents of editmes as a double


% --- Executes during object creation, after setting all properties.
function editmes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editmes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editserial_Callback(hObject, eventdata, handles)
% hObject    handle to editserial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editserial as text
%        str2double(get(hObject,'String')) returns contents of editserial as a double


% --- Executes during object creation, after setting all properties.
function editserial_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editserial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editother_Callback(hObject, eventdata, handles)
% hObject    handle to editother (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editother as text
%        str2double(get(hObject,'String')) returns contents of editother as a double


% --- Executes during object creation, after setting all properties.
function editother_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editother (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rc_Callback(hObject, eventdata, handles)
% hObject    handle to rc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rc as text
%        str2double(get(hObject,'String')) returns contents of rc as a double


% --- Executes during object creation, after setting all properties.
function rc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ra_Callback(hObject, eventdata, handles)
% hObject    handle to ra (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ra as text
%        str2double(get(hObject,'String')) returns contents of ra as a double


% --- Executes during object creation, after setting all properties.
function ra_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ra (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function epsilon_Callback(hObject, eventdata, handles)
% hObject    handle to epsilon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of epsilon as text
%        str2double(get(hObject,'String')) returns contents of epsilon as a double


% --- Executes during object creation, after setting all properties.
function epsilon_CreateFcn(hObject, eventdata, handles)
% hObject    handle to epsilon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function beta_Callback(hObject, eventdata, handles)
% hObject    handle to beta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of beta as text
%        str2double(get(hObject,'String')) returns contents of beta as a double


% --- Executes during object creation, after setting all properties.
function beta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to beta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function phi_Callback(hObject, eventdata, handles)
% hObject    handle to phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of phi as text
%        str2double(get(hObject,'String')) returns contents of phi as a double


% --- Executes during object creation, after setting all properties.
function phi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tc_Callback(hObject, eventdata, handles)
% hObject    handle to tc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tc as text
%        str2double(get(hObject,'String')) returns contents of tc as a double


% --- Executes during object creation, after setting all properties.
function tc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function idc_Callback(hObject, eventdata, handles)
% hObject    handle to idc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of idc as text
%        str2double(get(hObject,'String')) returns contents of idc as a double


% --- Executes during object creation, after setting all properties.
function idc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to idc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vleft_Callback(hObject, eventdata, handles)
% hObject    handle to vleft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vleft as text
%        str2double(get(hObject,'String')) returns contents of vleft as a double


% --- Executes during object creation, after setting all properties.
function vleft_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vleft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vright_Callback(hObject, eventdata, handles)
% hObject    handle to vright (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vright as text
%        str2double(get(hObject,'String')) returns contents of vright as a double


% --- Executes during object creation, after setting all properties.
function vright_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vright (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editrunstep_Callback(hObject, eventdata, handles)
% hObject    handle to editrunstep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editrunstep as text
%        str2double(get(hObject,'String')) returns contents of editrunstep as a double


% --- Executes during object creation, after setting all properties.
function editrunstep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editrunstep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in protable_popmen.
function protable_popmen_Callback(hObject, eventdata, handles)
% hObject    handle to protable_popmen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns protable_popmen contents as cell array
%        contents{get(hObject,'Value')} returns selected item from protable_popmen


% --- Executes during object creation, after setting all properties.
function protable_popmen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to protable_popmen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in bl.
function bl_Callback(hObject, eventdata, handles)
% hObject    handle to bl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global comobj blobjname
    % Indicate battery level of robot (reading from com port)
    
IDstr=get(handles.idc,'string');
ID=str2num(IDstr);    
frame=['p',ID];
fwrite(comobj,frame);
flushoutput(comobj);

if comobj.BytesAvailable>0
    byte12=fread(comobj,2);
    if (byte12(1)=='A')&(byte12(2)==ID) 
        byte34=fread(comobj,2);
        millivoltvalue=(byte34(1)+byte34(2)*256)*3/2;       
        if millivoltvalue<4000
            set(blobjname(byte12(2)),'string',num2str(millivoltvalue,'%4.1f'),'BackgroundColor','r');
        else
            set(blobjname(byte12(2)),'string',num2str(millivoltvalue,'%4.1f'),'BackgroundColor','w');
        end
    end
end

% --- Executes on button press in led.
function led_Callback(hObject, eventdata, handles)
% hObject    handle to led (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% LEDID=11 (Green); LEDID=12 (Red); LEDID=17 (Blue)
% LEDStatus = 0 (on); LEDStatus = 0 (off) 
global comobj
IDstr=get(handles.idc,'string');
idr=str2num(IDstr);
idled=get(handles.LedID,'string');
status=get(handles.LedStatus,'string');
LED(IvIDSet(idr),idled,status)
flushoutput(comobj);



% --- Executes on button press in buzzer.
function buzzer_Callback(hObject, eventdata, handles)
% hObject    handle to buzzer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global comobj

IDstr=get(handles.idc,'string');
ID=str2num(IDstr);
frame=['s',ID];
fwrite(comobj,frame);
flushoutput(comobj);



function LedID_Callback(hObject, eventdata, handles)
% hObject    handle to LedID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LedID as text
%        str2double(get(hObject,'String')) returns contents of LedID as a double


% --- Executes during object creation, after setting all properties.
function LedID_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LedID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LedStatus_Callback(hObject, eventdata, handles)
% hObject    handle to LedStatus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LedStatus as text
%        str2double(get(hObject,'String')) returns contents of LedStatus as a double


% --- Executes during object creation, after setting all properties.
function LedStatus_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LedStatus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function id11_Callback(hObject, eventdata, handles)
% hObject    handle to id11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id11 as text
%        str2double(get(hObject,'String')) returns contents of id11 as a double


% --- Executes during object creation, after setting all properties.
function id11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px11_Callback(hObject, eventdata, handles)
% hObject    handle to px11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px11 as text
%        str2double(get(hObject,'String')) returns contents of px11 as a double


% --- Executes during object creation, after setting all properties.
function px11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py11_Callback(hObject, eventdata, handles)
% hObject    handle to py11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py11 as text
%        str2double(get(hObject,'String')) returns contents of py11 as a double


% --- Executes during object creation, after setting all properties.
function py11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or11_Callback(hObject, eventdata, handles)
% hObject    handle to or11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or11 as text
%        str2double(get(hObject,'String')) returns contents of or11 as a double


% --- Executes during object creation, after setting all properties.
function or11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl11_Callback(hObject, eventdata, handles)
% hObject    handle to bl11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl11 as text
%        str2double(get(hObject,'String')) returns contents of bl11 as a double


% --- Executes during object creation, after setting all properties.
function bl11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function start_button_CreateFcn(hObject, eventdata, handles)
% hObject    handle to start_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function id12_Callback(hObject, eventdata, handles)
% hObject    handle to id12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id12 as text
%        str2double(get(hObject,'String')) returns contents of id12 as a double


% --- Executes during object creation, after setting all properties.
function id12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px12_Callback(hObject, eventdata, handles)
% hObject    handle to px12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px12 as text
%        str2double(get(hObject,'String')) returns contents of px12 as a double


% --- Executes during object creation, after setting all properties.
function px12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py12_Callback(hObject, eventdata, handles)
% hObject    handle to py12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py12 as text
%        str2double(get(hObject,'String')) returns contents of py12 as a double


% --- Executes during object creation, after setting all properties.
function py12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or12_Callback(hObject, eventdata, handles)
% hObject    handle to or12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or12 as text
%        str2double(get(hObject,'String')) returns contents of or12 as a double


% --- Executes during object creation, after setting all properties.
function or12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl12_Callback(hObject, eventdata, handles)
% hObject    handle to bl12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl12 as text
%        str2double(get(hObject,'String')) returns contents of bl12 as a double


% --- Executes during object creation, after setting all properties.
function bl12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function id13_Callback(hObject, eventdata, handles)
% hObject    handle to id13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id13 as text
%        str2double(get(hObject,'String')) returns contents of id13 as a double


% --- Executes during object creation, after setting all properties.
function id13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px13_Callback(hObject, eventdata, handles)
% hObject    handle to px13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px13 as text
%        str2double(get(hObject,'String')) returns contents of px13 as a double


% --- Executes during object creation, after setting all properties.
function px13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py13_Callback(hObject, eventdata, handles)
% hObject    handle to py13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py13 as text
%        str2double(get(hObject,'String')) returns contents of py13 as a double


% --- Executes during object creation, after setting all properties.
function py13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or13_Callback(hObject, eventdata, handles)
% hObject    handle to or13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or13 as text
%        str2double(get(hObject,'String')) returns contents of or13 as a double


% --- Executes during object creation, after setting all properties.
function or13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl13_Callback(hObject, eventdata, handles)
% hObject    handle to bl13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl13 as text
%        str2double(get(hObject,'String')) returns contents of bl13 as a double


% --- Executes during object creation, after setting all properties.
function bl13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function id14_Callback(hObject, eventdata, handles)
% hObject    handle to id14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id14 as text
%        str2double(get(hObject,'String')) returns contents of id14 as a double


% --- Executes during object creation, after setting all properties.
function id14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px14_Callback(hObject, eventdata, handles)
% hObject    handle to px14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px14 as text
%        str2double(get(hObject,'String')) returns contents of px14 as a double


% --- Executes during object creation, after setting all properties.
function px14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py14_Callback(hObject, eventdata, handles)
% hObject    handle to py14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py14 as text
%        str2double(get(hObject,'String')) returns contents of py14 as a double


% --- Executes during object creation, after setting all properties.
function py14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or14_Callback(hObject, eventdata, handles)
% hObject    handle to or14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or14 as text
%        str2double(get(hObject,'String')) returns contents of or14 as a double


% --- Executes during object creation, after setting all properties.
function or14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl14_Callback(hObject, eventdata, handles)
% hObject    handle to bl14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl14 as text
%        str2double(get(hObject,'String')) returns contents of bl14 as a double


% --- Executes during object creation, after setting all properties.
function bl14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function id15_Callback(hObject, eventdata, handles)
% hObject    handle to id15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id15 as text
%        str2double(get(hObject,'String')) returns contents of id15 as a double


% --- Executes during object creation, after setting all properties.
function id15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px15_Callback(hObject, eventdata, handles)
% hObject    handle to px15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px15 as text
%        str2double(get(hObject,'String')) returns contents of px15 as a double


% --- Executes during object creation, after setting all properties.
function px15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py15_Callback(hObject, eventdata, handles)
% hObject    handle to py15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py15 as text
%        str2double(get(hObject,'String')) returns contents of py15 as a double


% --- Executes during object creation, after setting all properties.
function py15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or15_Callback(hObject, eventdata, handles)
% hObject    handle to or15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or15 as text
%        str2double(get(hObject,'String')) returns contents of or15 as a double


% --- Executes during object creation, after setting all properties.
function or15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl15_Callback(hObject, eventdata, handles)
% hObject    handle to bl15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl15 as text
%        str2double(get(hObject,'String')) returns contents of bl15 as a double


% --- Executes during object creation, after setting all properties.
function bl15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function id16_Callback(hObject, eventdata, handles)
% hObject    handle to id16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id16 as text
%        str2double(get(hObject,'String')) returns contents of id16 as a double


% --- Executes during object creation, after setting all properties.
function id16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px16_Callback(hObject, eventdata, handles)
% hObject    handle to px16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px16 as text
%        str2double(get(hObject,'String')) returns contents of px16 as a double


% --- Executes during object creation, after setting all properties.
function px16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py16_Callback(hObject, eventdata, handles)
% hObject    handle to py16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py16 as text
%        str2double(get(hObject,'String')) returns contents of py16 as a double


% --- Executes during object creation, after setting all properties.
function py16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or16_Callback(hObject, eventdata, handles)
% hObject    handle to or16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or16 as text
%        str2double(get(hObject,'String')) returns contents of or16 as a double


% --- Executes during object creation, after setting all properties.
function or16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl16_Callback(hObject, eventdata, handles)
% hObject    handle to bl16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl16 as text
%        str2double(get(hObject,'String')) returns contents of bl16 as a double


% --- Executes during object creation, after setting all properties.
function bl16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function id17_Callback(hObject, eventdata, handles)
% hObject    handle to id17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id17 as text
%        str2double(get(hObject,'String')) returns contents of id17 as a double


% --- Executes during object creation, after setting all properties.
function id17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px17_Callback(hObject, eventdata, handles)
% hObject    handle to px17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px17 as text
%        str2double(get(hObject,'String')) returns contents of px17 as a double


% --- Executes during object creation, after setting all properties.
function px17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py17_Callback(hObject, eventdata, handles)
% hObject    handle to py17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py17 as text
%        str2double(get(hObject,'String')) returns contents of py17 as a double


% --- Executes during object creation, after setting all properties.
function py17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or17_Callback(hObject, eventdata, handles)
% hObject    handle to or17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or17 as text
%        str2double(get(hObject,'String')) returns contents of or17 as a double


% --- Executes during object creation, after setting all properties.
function or17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl17_Callback(hObject, eventdata, handles)
% hObject    handle to bl17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl17 as text
%        str2double(get(hObject,'String')) returns contents of bl17 as a double


% --- Executes during object creation, after setting all properties.
function bl17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function id18_Callback(hObject, eventdata, handles)
% hObject    handle to id18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id18 as text
%        str2double(get(hObject,'String')) returns contents of id18 as a double


% --- Executes during object creation, after setting all properties.
function id18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px18_Callback(hObject, eventdata, handles)
% hObject    handle to px18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px18 as text
%        str2double(get(hObject,'String')) returns contents of px18 as a double


% --- Executes during object creation, after setting all properties.
function px18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py18_Callback(hObject, eventdata, handles)
% hObject    handle to py18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py18 as text
%        str2double(get(hObject,'String')) returns contents of py18 as a double


% --- Executes during object creation, after setting all properties.
function py18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or18_Callback(hObject, eventdata, handles)
% hObject    handle to or18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or18 as text
%        str2double(get(hObject,'String')) returns contents of or18 as a double


% --- Executes during object creation, after setting all properties.
function or18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl18_Callback(hObject, eventdata, handles)
% hObject    handle to bl18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl18 as text
%        str2double(get(hObject,'String')) returns contents of bl18 as a double


% --- Executes during object creation, after setting all properties.
function bl18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function id19_Callback(hObject, eventdata, handles)
% hObject    handle to id19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id19 as text
%        str2double(get(hObject,'String')) returns contents of id19 as a double


% --- Executes during object creation, after setting all properties.
function id19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px19_Callback(hObject, eventdata, handles)
% hObject    handle to px19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px19 as text
%        str2double(get(hObject,'String')) returns contents of px19 as a double


% --- Executes during object creation, after setting all properties.
function px19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py19_Callback(hObject, eventdata, handles)
% hObject    handle to py19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py19 as text
%        str2double(get(hObject,'String')) returns contents of py19 as a double


% --- Executes during object creation, after setting all properties.
function py19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or19_Callback(hObject, eventdata, handles)
% hObject    handle to or19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or19 as text
%        str2double(get(hObject,'String')) returns contents of or19 as a double


% --- Executes during object creation, after setting all properties.
function or19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl19_Callback(hObject, eventdata, handles)
% hObject    handle to bl19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl19 as text
%        str2double(get(hObject,'String')) returns contents of bl19 as a double


% --- Executes during object creation, after setting all properties.
function bl19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function id20_Callback(hObject, eventdata, handles)
% hObject    handle to id20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id20 as text
%        str2double(get(hObject,'String')) returns contents of id20 as a double


% --- Executes during object creation, after setting all properties.
function id20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px20_Callback(hObject, eventdata, handles)
% hObject    handle to px20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px20 as text
%        str2double(get(hObject,'String')) returns contents of px20 as a double


% --- Executes during object creation, after setting all properties.
function px20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py20_Callback(hObject, eventdata, handles)
% hObject    handle to py20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py20 as text
%        str2double(get(hObject,'String')) returns contents of py20 as a double


% --- Executes during object creation, after setting all properties.
function py20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or20_Callback(hObject, eventdata, handles)
% hObject    handle to or20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or20 as text
%        str2double(get(hObject,'String')) returns contents of or20 as a double


% --- Executes during object creation, after setting all properties.
function or20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl20_Callback(hObject, eventdata, handles)
% hObject    handle to bl20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl20 as text
%        str2double(get(hObject,'String')) returns contents of bl20 as a double


% --- Executes during object creation, after setting all properties.
function bl20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function id21_Callback(hObject, eventdata, handles)
% hObject    handle to id21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id21 as text
%        str2double(get(hObject,'String')) returns contents of id21 as a double


% --- Executes during object creation, after setting all properties.
function id21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px21_Callback(hObject, eventdata, handles)
% hObject    handle to px21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px21 as text
%        str2double(get(hObject,'String')) returns contents of px21 as a double


% --- Executes during object creation, after setting all properties.
function px21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py21_Callback(hObject, eventdata, handles)
% hObject    handle to py21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py21 as text
%        str2double(get(hObject,'String')) returns contents of py21 as a double


% --- Executes during object creation, after setting all properties.
function py21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or21_Callback(hObject, eventdata, handles)
% hObject    handle to or21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or21 as text
%        str2double(get(hObject,'String')) returns contents of or21 as a double


% --- Executes during object creation, after setting all properties.
function or21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl21_Callback(hObject, eventdata, handles)
% hObject    handle to bl21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl21 as text
%        str2double(get(hObject,'String')) returns contents of bl21 as a double


% --- Executes during object creation, after setting all properties.
function bl21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function id22_Callback(hObject, eventdata, handles)
% hObject    handle to id22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id22 as text
%        str2double(get(hObject,'String')) returns contents of id22 as a double


% --- Executes during object creation, after setting all properties.
function id22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px22_Callback(hObject, eventdata, handles)
% hObject    handle to px22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px22 as text
%        str2double(get(hObject,'String')) returns contents of px22 as a double


% --- Executes during object creation, after setting all properties.
function px22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py22_Callback(hObject, eventdata, handles)
% hObject    handle to py22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py22 as text
%        str2double(get(hObject,'String')) returns contents of py22 as a double


% --- Executes during object creation, after setting all properties.
function py22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or22_Callback(hObject, eventdata, handles)
% hObject    handle to or22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or22 as text
%        str2double(get(hObject,'String')) returns contents of or22 as a double


% --- Executes during object creation, after setting all properties.
function or22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl22_Callback(hObject, eventdata, handles)
% hObject    handle to bl22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl22 as text
%        str2double(get(hObject,'String')) returns contents of bl22 as a double


% --- Executes during object creation, after setting all properties.
function bl22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function id23_Callback(hObject, eventdata, handles)
% hObject    handle to id23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id23 as text
%        str2double(get(hObject,'String')) returns contents of id23 as a double


% --- Executes during object creation, after setting all properties.
function id23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px23_Callback(hObject, eventdata, handles)
% hObject    handle to px23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px23 as text
%        str2double(get(hObject,'String')) returns contents of px23 as a double


% --- Executes during object creation, after setting all properties.
function px23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py23_Callback(hObject, eventdata, handles)
% hObject    handle to py23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py23 as text
%        str2double(get(hObject,'String')) returns contents of py23 as a double


% --- Executes during object creation, after setting all properties.
function py23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or23_Callback(hObject, eventdata, handles)
% hObject    handle to or23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or23 as text
%        str2double(get(hObject,'String')) returns contents of or23 as a double


% --- Executes during object creation, after setting all properties.
function or23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl23_Callback(hObject, eventdata, handles)
% hObject    handle to bl23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl23 as text
%        str2double(get(hObject,'String')) returns contents of bl23 as a double


% --- Executes during object creation, after setting all properties.
function bl23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function id24_Callback(hObject, eventdata, handles)
% hObject    handle to id24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id24 as text
%        str2double(get(hObject,'String')) returns contents of id24 as a double


% --- Executes during object creation, after setting all properties.
function id24_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px24_Callback(hObject, eventdata, handles)
% hObject    handle to px24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px24 as text
%        str2double(get(hObject,'String')) returns contents of px24 as a double


% --- Executes during object creation, after setting all properties.
function px24_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py24_Callback(hObject, eventdata, handles)
% hObject    handle to py24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py24 as text
%        str2double(get(hObject,'String')) returns contents of py24 as a double


% --- Executes during object creation, after setting all properties.
function py24_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or24_Callback(hObject, eventdata, handles)
% hObject    handle to or24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or24 as text
%        str2double(get(hObject,'String')) returns contents of or24 as a double


% --- Executes during object creation, after setting all properties.
function or24_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl24_Callback(hObject, eventdata, handles)
% hObject    handle to bl24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl24 as text
%        str2double(get(hObject,'String')) returns contents of bl24 as a double


% --- Executes during object creation, after setting all properties.
function bl24_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function id25_Callback(hObject, eventdata, handles)
% hObject    handle to id25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id25 as text
%        str2double(get(hObject,'String')) returns contents of id25 as a double


% --- Executes during object creation, after setting all properties.
function id25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px25_Callback(hObject, eventdata, handles)
% hObject    handle to px25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px25 as text
%        str2double(get(hObject,'String')) returns contents of px25 as a double


% --- Executes during object creation, after setting all properties.
function px25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py25_Callback(hObject, eventdata, handles)
% hObject    handle to py25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py25 as text
%        str2double(get(hObject,'String')) returns contents of py25 as a double


% --- Executes during object creation, after setting all properties.
function py25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or25_Callback(hObject, eventdata, handles)
% hObject    handle to or25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or25 as text
%        str2double(get(hObject,'String')) returns contents of or25 as a double


% --- Executes during object creation, after setting all properties.
function or25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl25_Callback(hObject, eventdata, handles)
% hObject    handle to bl25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl25 as text
%        str2double(get(hObject,'String')) returns contents of bl25 as a double


% --- Executes during object creation, after setting all properties.
function bl25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function id26_Callback(hObject, eventdata, handles)
% hObject    handle to id26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id26 as text
%        str2double(get(hObject,'String')) returns contents of id26 as a double


% --- Executes during object creation, after setting all properties.
function id26_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px26_Callback(hObject, eventdata, handles)
% hObject    handle to px26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px26 as text
%        str2double(get(hObject,'String')) returns contents of px26 as a double


% --- Executes during object creation, after setting all properties.
function px26_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py26_Callback(hObject, eventdata, handles)
% hObject    handle to py26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py26 as text
%        str2double(get(hObject,'String')) returns contents of py26 as a double


% --- Executes during object creation, after setting all properties.
function py26_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or26_Callback(hObject, eventdata, handles)
% hObject    handle to or26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or26 as text
%        str2double(get(hObject,'String')) returns contents of or26 as a double


% --- Executes during object creation, after setting all properties.
function or26_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl26_Callback(hObject, eventdata, handles)
% hObject    handle to bl26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl26 as text
%        str2double(get(hObject,'String')) returns contents of bl26 as a double


% --- Executes during object creation, after setting all properties.
function bl26_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function id27_Callback(hObject, eventdata, handles)
% hObject    handle to id27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id27 as text
%        str2double(get(hObject,'String')) returns contents of id27 as a double


% --- Executes during object creation, after setting all properties.
function id27_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px27_Callback(hObject, eventdata, handles)
% hObject    handle to px27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px27 as text
%        str2double(get(hObject,'String')) returns contents of px27 as a double


% --- Executes during object creation, after setting all properties.
function px27_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py27_Callback(hObject, eventdata, handles)
% hObject    handle to py27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py27 as text
%        str2double(get(hObject,'String')) returns contents of py27 as a double


% --- Executes during object creation, after setting all properties.
function py27_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or27_Callback(hObject, eventdata, handles)
% hObject    handle to or27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or27 as text
%        str2double(get(hObject,'String')) returns contents of or27 as a double


% --- Executes during object creation, after setting all properties.
function or27_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl27_Callback(hObject, eventdata, handles)
% hObject    handle to bl27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl27 as text
%        str2double(get(hObject,'String')) returns contents of bl27 as a double


% --- Executes during object creation, after setting all properties.
function bl27_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function id28_Callback(hObject, eventdata, handles)
% hObject    handle to id28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id28 as text
%        str2double(get(hObject,'String')) returns contents of id28 as a double


% --- Executes during object creation, after setting all properties.
function id28_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px28_Callback(hObject, eventdata, handles)
% hObject    handle to px28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px28 as text
%        str2double(get(hObject,'String')) returns contents of px28 as a double


% --- Executes during object creation, after setting all properties.
function px28_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py28_Callback(hObject, eventdata, handles)
% hObject    handle to py28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py28 as text
%        str2double(get(hObject,'String')) returns contents of py28 as a double


% --- Executes during object creation, after setting all properties.
function py28_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or28_Callback(hObject, eventdata, handles)
% hObject    handle to or28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or28 as text
%        str2double(get(hObject,'String')) returns contents of or28 as a double


% --- Executes during object creation, after setting all properties.
function or28_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl28_Callback(hObject, eventdata, handles)
% hObject    handle to bl28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl28 as text
%        str2double(get(hObject,'String')) returns contents of bl28 as a double


% --- Executes during object creation, after setting all properties.
function bl28_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function id29_Callback(hObject, eventdata, handles)
% hObject    handle to id29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id29 as text
%        str2double(get(hObject,'String')) returns contents of id29 as a double


% --- Executes during object creation, after setting all properties.
function id29_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px29_Callback(hObject, eventdata, handles)
% hObject    handle to px29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px29 as text
%        str2double(get(hObject,'String')) returns contents of px29 as a double


% --- Executes during object creation, after setting all properties.
function px29_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py29_Callback(hObject, eventdata, handles)
% hObject    handle to py29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py29 as text
%        str2double(get(hObject,'String')) returns contents of py29 as a double


% --- Executes during object creation, after setting all properties.
function py29_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or29_Callback(hObject, eventdata, handles)
% hObject    handle to or29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or29 as text
%        str2double(get(hObject,'String')) returns contents of or29 as a double


% --- Executes during object creation, after setting all properties.
function or29_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl29_Callback(hObject, eventdata, handles)
% hObject    handle to bl29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl29 as text
%        str2double(get(hObject,'String')) returns contents of bl29 as a double


% --- Executes during object creation, after setting all properties.
function bl29_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function id30_Callback(hObject, eventdata, handles)
% hObject    handle to id30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id30 as text
%        str2double(get(hObject,'String')) returns contents of id30 as a double


% --- Executes during object creation, after setting all properties.
function id30_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px30_Callback(hObject, eventdata, handles)
% hObject    handle to px30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px30 as text
%        str2double(get(hObject,'String')) returns contents of px30 as a double


% --- Executes during object creation, after setting all properties.
function px30_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py30_Callback(hObject, eventdata, handles)
% hObject    handle to py30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py30 as text
%        str2double(get(hObject,'String')) returns contents of py30 as a double


% --- Executes during object creation, after setting all properties.
function py30_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or30_Callback(hObject, eventdata, handles)
% hObject    handle to or30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or30 as text
%        str2double(get(hObject,'String')) returns contents of or30 as a double


% --- Executes during object creation, after setting all properties.
function or30_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl30_Callback(hObject, eventdata, handles)
% hObject    handle to bl30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl30 as text
%        str2double(get(hObject,'String')) returns contents of bl30 as a double


% --- Executes during object creation, after setting all properties.
function bl30_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function id31_Callback(hObject, eventdata, handles)
% hObject    handle to id31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id31 as text
%        str2double(get(hObject,'String')) returns contents of id31 as a double


% --- Executes during object creation, after setting all properties.
function id31_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px31_Callback(hObject, eventdata, handles)
% hObject    handle to px31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px31 as text
%        str2double(get(hObject,'String')) returns contents of px31 as a double


% --- Executes during object creation, after setting all properties.
function px31_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py31_Callback(hObject, eventdata, handles)
% hObject    handle to py31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py31 as text
%        str2double(get(hObject,'String')) returns contents of py31 as a double


% --- Executes during object creation, after setting all properties.
function py31_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or31_Callback(hObject, eventdata, handles)
% hObject    handle to or31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or31 as text
%        str2double(get(hObject,'String')) returns contents of or31 as a double


% --- Executes during object creation, after setting all properties.
function or31_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl31_Callback(hObject, eventdata, handles)
% hObject    handle to bl31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl31 as text
%        str2double(get(hObject,'String')) returns contents of bl31 as a double


% --- Executes during object creation, after setting all properties.
function bl31_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function id32_Callback(hObject, eventdata, handles)
% hObject    handle to id32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id32 as text
%        str2double(get(hObject,'String')) returns contents of id32 as a double


% --- Executes during object creation, after setting all properties.
function id32_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px32_Callback(hObject, eventdata, handles)
% hObject    handle to px32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px32 as text
%        str2double(get(hObject,'String')) returns contents of px32 as a double


% --- Executes during object creation, after setting all properties.
function px32_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py32_Callback(hObject, eventdata, handles)
% hObject    handle to py32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py32 as text
%        str2double(get(hObject,'String')) returns contents of py32 as a double


% --- Executes during object creation, after setting all properties.
function py32_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or32_Callback(hObject, eventdata, handles)
% hObject    handle to or32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or32 as text
%        str2double(get(hObject,'String')) returns contents of or32 as a double


% --- Executes during object creation, after setting all properties.
function or32_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl32_Callback(hObject, eventdata, handles)
% hObject    handle to bl32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl32 as text
%        str2double(get(hObject,'String')) returns contents of bl32 as a double


% --- Executes during object creation, after setting all properties.
function bl32_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function id33_Callback(hObject, eventdata, handles)
% hObject    handle to id33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id33 as text
%        str2double(get(hObject,'String')) returns contents of id33 as a double


% --- Executes during object creation, after setting all properties.
function id33_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px33_Callback(hObject, eventdata, handles)
% hObject    handle to px33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px33 as text
%        str2double(get(hObject,'String')) returns contents of px33 as a double


% --- Executes during object creation, after setting all properties.
function px33_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py33_Callback(hObject, eventdata, handles)
% hObject    handle to py33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py33 as text
%        str2double(get(hObject,'String')) returns contents of py33 as a double


% --- Executes during object creation, after setting all properties.
function py33_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or33_Callback(hObject, eventdata, handles)
% hObject    handle to or33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or33 as text
%        str2double(get(hObject,'String')) returns contents of or33 as a double


% --- Executes during object creation, after setting all properties.
function or33_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl33_Callback(hObject, eventdata, handles)
% hObject    handle to bl33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl33 as text
%        str2double(get(hObject,'String')) returns contents of bl33 as a double


% --- Executes during object creation, after setting all properties.
function bl33_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function id34_Callback(hObject, eventdata, handles)
% hObject    handle to id34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id34 as text
%        str2double(get(hObject,'String')) returns contents of id34 as a double


% --- Executes during object creation, after setting all properties.
function id34_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px34_Callback(hObject, eventdata, handles)
% hObject    handle to px34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px34 as text
%        str2double(get(hObject,'String')) returns contents of px34 as a double


% --- Executes during object creation, after setting all properties.
function px34_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py34_Callback(hObject, eventdata, handles)
% hObject    handle to py34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py34 as text
%        str2double(get(hObject,'String')) returns contents of py34 as a double


% --- Executes during object creation, after setting all properties.
function py34_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or34_Callback(hObject, eventdata, handles)
% hObject    handle to or34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or34 as text
%        str2double(get(hObject,'String')) returns contents of or34 as a double


% --- Executes during object creation, after setting all properties.
function or34_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl34_Callback(hObject, eventdata, handles)
% hObject    handle to bl34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl34 as text
%        str2double(get(hObject,'String')) returns contents of bl34 as a double


% --- Executes during object creation, after setting all properties.
function bl34_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function id35_Callback(hObject, eventdata, handles)
% hObject    handle to id35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id35 as text
%        str2double(get(hObject,'String')) returns contents of id35 as a double


% --- Executes during object creation, after setting all properties.
function id35_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px35_Callback(hObject, eventdata, handles)
% hObject    handle to px35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px35 as text
%        str2double(get(hObject,'String')) returns contents of px35 as a double


% --- Executes during object creation, after setting all properties.
function px35_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py35_Callback(hObject, eventdata, handles)
% hObject    handle to py35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py35 as text
%        str2double(get(hObject,'String')) returns contents of py35 as a double


% --- Executes during object creation, after setting all properties.
function py35_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or35_Callback(hObject, eventdata, handles)
% hObject    handle to or35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or35 as text
%        str2double(get(hObject,'String')) returns contents of or35 as a double


% --- Executes during object creation, after setting all properties.
function or35_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl35_Callback(hObject, eventdata, handles)
% hObject    handle to bl35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl35 as text
%        str2double(get(hObject,'String')) returns contents of bl35 as a double


% --- Executes during object creation, after setting all properties.
function bl35_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function id36_Callback(hObject, eventdata, handles)
% hObject    handle to id36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id36 as text
%        str2double(get(hObject,'String')) returns contents of id36 as a double


% --- Executes during object creation, after setting all properties.
function id36_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px36_Callback(hObject, eventdata, handles)
% hObject    handle to px36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px36 as text
%        str2double(get(hObject,'String')) returns contents of px36 as a double


% --- Executes during object creation, after setting all properties.
function px36_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py36_Callback(hObject, eventdata, handles)
% hObject    handle to py36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py36 as text
%        str2double(get(hObject,'String')) returns contents of py36 as a double


% --- Executes during object creation, after setting all properties.
function py36_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or36_Callback(hObject, eventdata, handles)
% hObject    handle to or36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or36 as text
%        str2double(get(hObject,'String')) returns contents of or36 as a double


% --- Executes during object creation, after setting all properties.
function or36_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl36_Callback(hObject, eventdata, handles)
% hObject    handle to bl36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl36 as text
%        str2double(get(hObject,'String')) returns contents of bl36 as a double


% --- Executes during object creation, after setting all properties.
function bl36_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function id37_Callback(hObject, eventdata, handles)
% hObject    handle to id37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id37 as text
%        str2double(get(hObject,'String')) returns contents of id37 as a double


% --- Executes during object creation, after setting all properties.
function id37_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px37_Callback(hObject, eventdata, handles)
% hObject    handle to px37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px37 as text
%        str2double(get(hObject,'String')) returns contents of px37 as a double


% --- Executes during object creation, after setting all properties.
function px37_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py37_Callback(hObject, eventdata, handles)
% hObject    handle to py37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py37 as text
%        str2double(get(hObject,'String')) returns contents of py37 as a double


% --- Executes during object creation, after setting all properties.
function py37_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or37_Callback(hObject, eventdata, handles)
% hObject    handle to or37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or37 as text
%        str2double(get(hObject,'String')) returns contents of or37 as a double


% --- Executes during object creation, after setting all properties.
function or37_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl37_Callback(hObject, eventdata, handles)
% hObject    handle to bl37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl37 as text
%        str2double(get(hObject,'String')) returns contents of bl37 as a double


% --- Executes during object creation, after setting all properties.
function bl37_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function id38_Callback(hObject, eventdata, handles)
% hObject    handle to id38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id38 as text
%        str2double(get(hObject,'String')) returns contents of id38 as a double


% --- Executes during object creation, after setting all properties.
function id38_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px38_Callback(hObject, eventdata, handles)
% hObject    handle to px38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px38 as text
%        str2double(get(hObject,'String')) returns contents of px38 as a double


% --- Executes during object creation, after setting all properties.
function px38_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py38_Callback(hObject, eventdata, handles)
% hObject    handle to py38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py38 as text
%        str2double(get(hObject,'String')) returns contents of py38 as a double


% --- Executes during object creation, after setting all properties.
function py38_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or38_Callback(hObject, eventdata, handles)
% hObject    handle to or38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or38 as text
%        str2double(get(hObject,'String')) returns contents of or38 as a double


% --- Executes during object creation, after setting all properties.
function or38_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl38_Callback(hObject, eventdata, handles)
% hObject    handle to bl38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl38 as text
%        str2double(get(hObject,'String')) returns contents of bl38 as a double


% --- Executes during object creation, after setting all properties.
function bl38_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function id39_Callback(hObject, eventdata, handles)
% hObject    handle to id39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id39 as text
%        str2double(get(hObject,'String')) returns contents of id39 as a double


% --- Executes during object creation, after setting all properties.
function id39_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px39_Callback(hObject, eventdata, handles)
% hObject    handle to px39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px39 as text
%        str2double(get(hObject,'String')) returns contents of px39 as a double


% --- Executes during object creation, after setting all properties.
function px39_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py39_Callback(hObject, eventdata, handles)
% hObject    handle to py39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py39 as text
%        str2double(get(hObject,'String')) returns contents of py39 as a double


% --- Executes during object creation, after setting all properties.
function py39_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or39_Callback(hObject, eventdata, handles)
% hObject    handle to or39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or39 as text
%        str2double(get(hObject,'String')) returns contents of or39 as a double


% --- Executes during object creation, after setting all properties.
function or39_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl39_Callback(hObject, eventdata, handles)
% hObject    handle to bl39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl39 as text
%        str2double(get(hObject,'String')) returns contents of bl39 as a double


% --- Executes during object creation, after setting all properties.
function bl39_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function id40_Callback(hObject, eventdata, handles)
% hObject    handle to id40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of id40 as text
%        str2double(get(hObject,'String')) returns contents of id40 as a double


% --- Executes during object creation, after setting all properties.
function id40_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px40_Callback(hObject, eventdata, handles)
% hObject    handle to px40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px40 as text
%        str2double(get(hObject,'String')) returns contents of px40 as a double


% --- Executes during object creation, after setting all properties.
function px40_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py40_Callback(hObject, eventdata, handles)
% hObject    handle to py40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py40 as text
%        str2double(get(hObject,'String')) returns contents of py40 as a double


% --- Executes during object creation, after setting all properties.
function py40_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function or40_Callback(hObject, eventdata, handles)
% hObject    handle to or40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of or40 as text
%        str2double(get(hObject,'String')) returns contents of or40 as a double


% --- Executes during object creation, after setting all properties.
function or40_CreateFcn(hObject, eventdata, handles)
% hObject    handle to or40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bl40_Callback(hObject, eventdata, handles)
% hObject    handle to bl40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bl40 as text
%        str2double(get(hObject,'String')) returns contents of bl40 as a double


% --- Executes during object creation, after setting all properties.
function bl40_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bl40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global xp Robot numberofRobots

[FileName, PathName] = uiputfile('*.txt', 'Save As');
Name = fullfile(PathName,FileName);
if PathName==0,
    return; 
end
fid=fopen(Name,'w');


if ~isempty(xp)
    for i=1:size(xp,1),
        fprintf(fid,'%1.4f %1.4f\n',xp(i,1),xp(i,2));
    end
    fclose(fid); 
else
    if numberofRobots>0
        for i=1:numberofRobots,
            fprintf(fid,'%1.4f %1.4f %1.4f\n',Robot(i).x(1),Robot(i).x(2), Robot(i).angle);
        end
        fclose(fid); 
    end
end

% --- Executes on button press in analysis.
function analysis_Callback(hObject, eventdata, handles)
% hObject    handle to analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global numberofRobots IDSet ra
[file path]=uigetfile('*.txt','Open');
if path==0,
    return; 
end
filename = fullfile(path,file);
data=load(filename);
for id=1:numberofRobots
    posm{id}=[];
    for i=1:size(data,1),
        if data(i,1)==id
            Pos1=[data(i,2),data(i,3)];
            Pos2=[data(i,5),data(i,6)];
            posm{id}=[posm{id},norm(Pos1-Pos2)];            
        end
    end
end

for id=1:numberofRobots,
    figure;
    hold on;
    x=[1:1:size(posm{id},2)];
    y=ra*ones(1,size(posm{id},2));
    plot(x,posm{id},'.r');
    plot(x,y,'-b');
    %xlabel('Collision Avoidance Parameter of Robot %d',IDSet(id));
end



function rh_Callback(hObject, eventdata, handles)
% hObject    handle to rh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rh as text
%        str2double(get(hObject,'String')) returns contents of rh as a double


% --- Executes during object creation, after setting all properties.
function rh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Trackingpoints.
function Trackingpoints_Callback(hObject, eventdata, handles)
% hObject    handle to Trackingpoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global idc 

idc=1;
Tc=0.1;
t4=timerfind('Name','timer4');
if isempty(t4)
    t4=timer('Name','timer4','TimerFcn',@trackingpoints,'Period',Tc,'ExecutionMode','fixedSpacing');
else
    stop(t4);
end
start(t4);


function trackingpoints(mTimer,~)
global G Robot idc xp idc_obj IDSet numberofRobots
global mode vmax


G = generateGraph();
id=str2num(get(idc_obj,'string'));

for i=1:numberofRobots
   if id==IDSet(i),
       k=i;
   end
end
xdk=xp(idc,:);


vmax=0.5;
Robot(k).v=vmax*(xdk-Robot(k).x)/norm(xdk-Robot(k).x);

v=norm(Robot(k).v);
if norm(norm(Robot(k).v))>0
    u=[1,0];    
    angle=acos(dot(Robot(k).v,u)/norm(Robot(k).v));
else
    angle=Robot(k).angle;
end

if Robot(k).v(2)<0     
    angle=-angle;
end
theta=diffangle(angle,Robot(k).angle);
Robot(k).alpha=theta;

%L2AVel_Circular(k,v)

switch mode
    case 1
        [dleft,dright,vleft,vright]=L2AVel(k);
        CommandSend(k,dleft,dright,vleft,vright);
    case 2
        [wl,wr]=L2AVel_simu(k);
        CommandSend_simu(k,[wl;wr]);        
end
if norm(Robot(k).x-xdk)<=0.1,
    if idc<size(xp,1),
       idc=idc+1; 
    else
        RobotStop(k)
        t4=timerfind('Name','timer4');
        if ~isempty(t4)
            stop(t4);
        end        
    end
end

function victim_Callback(hObject, eventdata, handles)
% hObject    handle to victim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of victim as text
%        str2double(get(hObject,'String')) returns contents of victim as a double


% --- Executes during object creation, after setting all properties.
function victim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to victim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in getpoints.
function getpoints_Callback(hObject, eventdata, handles)
% hObject    handle to getpoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global xp hf
xp=[]
set(hf,'ButtonDownFcn',@getpointsforHome);

function getpointsforHome(obj,evt)
global xp hf

pos=get(obj,'CurrentPoint');
xp=[xp;[pos(1,1),pos(1,2)]];
plot(hf,pos(1,1),pos(1,2),'or'); 


% --- Executes on button press in homeall.
function homeall_Callback(hObject, eventdata, handles)
% hObject    handle to homeall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global numberofRobots xd Robot homeset count runset c
xd=[];
count=0;
c=2;
for i=0:5
    for j=2:-1:-2,
        if size(xd,1)<numberofRobots,
            x=[0.15+0.3*i,-0.5+0.3*j];
            %x=[0.15+0.25*i,-0.3+0.25*j];
            %x=[6.3-0.25*i,-0.3+0.25*j];
            xd=[xd;x];
        end
    end
end
plot(xd(:,1),xd(:,2),'.r');
dist=[];
for id=1:numberofRobots,
    dist=[dist,norm(xd(1,:)-Robot(id).x)];
end

[value,homeset]=sort(dist);
runset=homeset(1);


Tc=0.04;
t5=timerfind('Name','timer5');
if isempty(t5)
    t5=timer('Name','timer5','TimerFcn',@homeallrobots,'Period',Tc,'ExecutionMode','fixedSpacing');
else
    stop(t5);
end
start(t5);


function homeallrobots(mTimer,~)
global G Robot homeset xd runset count c
global Pr wmax vmax

Pr
wmax
vmax

G = generateGraph();

count=count+1;
if mod(count,10)==0
   if c<=size(homeset,2)
       runset=[runset,homeset(c)];
       c=c+1;
   end
end

for i=1:size(runset,2),
    xdk=xd(i,:);
    id=runset(i);
     if norm(xdk-Robot(id).x)>0.1
         HDCforHome(id,xdk);
     else
        CommandSend(id,0,0,0,0);
%         LED(id,'G','off'); 
%         LED(id,'R','off'); 
%         LED(id,'B','off'); 
    end
end


% --- Executes on selection change in modeselection.
function modeselection_Callback(hObject, eventdata, handles)
% hObject    handle to modeselection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns modeselection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from modeselection


% --- Executes during object creation, after setting all properties.
function modeselection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to modeselection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in optionforRobotinit.
function optionforRobotinit_Callback(hObject, eventdata, handles)
% hObject    handle to optionforRobotinit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of optionforRobotinit



function NumofRobots_Callback(hObject, eventdata, handles)
% hObject    handle to NumofRobots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NumofRobots as text
%        str2double(get(hObject,'String')) returns contents of NumofRobots as a double


% --- Executes during object creation, after setting all properties.
function NumofRobots_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NumofRobots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in dtime.
function dtime_Callback(hObject, eventdata, handles)
% hObject    handle to dtime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fig=figure;
hold on;
xlim([0 300]);
ylim([0 0.5]);
datafolder=pwd;
path=[datafolder,'\ReData'];
k='1';
legendname=[];
for i=1:40,
    filename=fullfile(path,['datatime',num2str(i),'_',k,'.txt']);
    y=load(filename)';
    x=[1:1:size(y,2)];
    plot(x,y,'-','color',[0,0,0+0.025*i]);
end

xlabel('steps');
ylabel('Sampling time (s)');
if k=='0'
    title('Sampling time of Motion Tracking System');
else
    title('Sampling time of Motion Tracking System (including DT)');
end
figname=['Samplingtime',k];
figname=fullfile(path,figname);
set(fig,'PaperPosition',[-0.2 0 6.5 4]);
set(fig,'PaperSize',[6.0 4]);
saveas(fig,figname,'pdf'); 

% --- Executes on button press in pushbutton29.
function pushbutton29_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in dup.
function dup_Callback(hObject, eventdata, handles)
% hObject    handle to dup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global comobj
id=45;
cmd=2;
v1=1;
v2=0;
v3=0;
frame=['m',id,cmd,v1,v2,v3]
fwrite(comobj,frame); 


% --- Executes on button press in ddown.
function ddown_Callback(hObject, eventdata, handles)
% hObject    handle to ddown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global comobj
id=45;
cmd=2;
v1=0;
v2=0;
v3=0;
frame=['m',id,cmd,v1,v2,v3]
fwrite(comobj,frame);



function Hop_Callback(hObject, eventdata, handles)
% hObject    handle to Hop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Hop as text
%        str2double(get(hObject,'String')) returns contents of Hop as a double


% --- Executes during object creation, after setting all properties.
function Hop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Hop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in swarm.
function swarm_Callback(hObject, eventdata, handles)
% hObject    handle to swarm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns swarm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from swarm


% --- Executes during object creation, after setting all properties.
function swarm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to swarm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in video.
function video_Callback(hObject, eventdata, handles)
% hObject    handle to video (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of video



function vm_Callback(hObject, eventdata, handles)
% hObject    handle to vm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vm as text
%        str2double(get(hObject,'String')) returns contents of vm as a double


% --- Executes during object creation, after setting all properties.
function vm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Td_Callback(hObject, eventdata, handles)
% hObject    handle to Td (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Td as text
%        str2double(get(hObject,'String')) returns contents of Td as a double


% --- Executes during object creation, after setting all properties.
function Td_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Td (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Tp_Callback(hObject, eventdata, handles)
% hObject    handle to Tp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tp as text
%        str2double(get(hObject,'String')) returns contents of Tp as a double


% --- Executes during object creation, after setting all properties.
function Tp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
