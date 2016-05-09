
function [ro] = armgen(params)
% 
% Sensor-Joint Kinematic Equation Generator
% Emanuele Ruffaldi and Lorenzo Peppoloni @ PERCRO, Scuola Superiore Sant'Anna Pisa
% 2012-2014
%
% Generates the direct transformation and optionally the
% Jacobian of the kinematic equation measured by sensors. 
%
% The returned functions accept two vectorial arguments. The first is the
% state of the joints (q,dq,ddq) ordered by joint index, the second are all
% the parameters needed.
%
% The intermediate functions are instead:
%   h(x_(p+1),params,z_p)
%
% The output for joint ones are: (omega,ddx,m,domega,
%
% For the intermediate 
%
% The output function contains, per sensor: omega, acc, mag.
%
% Syntax:
%   params = armgen()
%        Returns all the default options 
%   ro = armgen(params)
%        Builds the structure
% Parameters:
%   params is a structure describing the structure of the kinematic chain
%   and the terms
%   .joints       joint tree expressed by the parent index (root is 0)
%   .sensors      index of the joint to which is attached solidally
%   .gaxis        axis of the gravity (set gaxis=0 to disable gravity)
%   .onlysym      generate only the symbolic version, not the numeric
%   .genJacob=1   build Jacobian
%   .DH           DH parameters as N x 4 (see below)
%   .filename     prefix of the filename of the generated functions
%   .multim0      marks multi m0
%   .genManual    manual generator (default)
%   .genInt       generate intermediate
%   .genJoint     disable generation of joints (in genInt=0)
%
% Default
%   "J##filename and h##filename"
%
% Intermediate Joints versions
%   "J%s%d##filename and h%s%d##filename"
%
%
% Each row of DH contains a, alpha, d and q0. If a coefficient of the
% matrix is specified by NaN it becomes a symbolic parameter. Remember also
% that it is better to use sym(pi) for the pi, because this helps the
% symbolic generator.
%
% This holds for the joints, while the sensors are connected to their
% parents by means of a transformation matrix:
%
%   S = rotz(param + q) roty(param) rotx(param) translate(param)
%
%
%
% Return: returns the structure in ro
%   .hsym  symbolic expression
%   .Jsym  symbolic expression  of jacobian
%   .h     same expression as function (only if filename is not specified)
%   .J     function handle (only if filename is not specified)
%   .allparams list of all the parameters that are passed as second
%         argument
%   .x     variables used as first parameter
%   .outputs list of outputs
%   .hjsym joint h functions for (omega,acc,mag)
%   .hj    joint h function explicit
%   .int   intermediate (structure,see below)
%
% Intermediate:
%   .sensors cell array (hsym,hnumeric) per joint
%   .joints cell array (hsym,hnumeric) per sensor
%
% TODO:
% - Export the full transformation matrices of joints and sensors (cell
% element #6)
% - Make measurable
%
% SIMPLEMODE: frames == joints
%
% joints specifies the joints and their parents
% DH are related to the jointlist
% sensors are related
%
% FRAMEMODE:
% joints NOT USED
% frames IS USED
% njoints IS NEEDED
% 
%
% Dependencies: none

if nargin == 0
    ro = [];
    ro.joints = [];
    ro.frames = [];
    ro.sensors = [];
    ro.njoints = 0;
    ro.gaxis = 2;
    ro.onlysym = 0;
    ro.genJacob = 1;
    ro.genManual = 1;
    ro.DH = [];
    ro.filename = '';
    ro.testsym = 0;
    ro.genJoints = 0;
    ro.multim0 = 0;
    ro.genInt = 0;
    ro.genJoint = 0;
    ro.T0 = eye(4);
    return
end


% ADJUST PARAMETERS
% ----------------------------------------
w = params;
defs = [];
defs.gaxis = 2;
defs.onlysym = 0;
defs.genJacob = 1;
defs.joints = [];
defs.multim0 = 0;
defs.genManual = 1;
defs.useexplicit = 0;
defs.genJoints = 0;
defs.njoints = 0;
defs.selectedsensoroutputs = [1,3,4]; % omega acc mag
defs.sensorstx = ones(length(w.sensors),1);
defs.sensortype = ones(1,length(w.sensors));
defs.genInt = 0;
defs.genJoint = 0;
defs.T0 = eye(4);
defs.testsym = 0;


fdefs = fieldnames(defs);
for I=1:length(fdefs)
    if isfield(w,fdefs{I}) == 0
        w.(fdefs{I}) = defs.(fdefs{I});
    end
end

assert((w.genInt && w.multim0) == 0,'MultiM0 is NOT compatible with Intermediate Z');

if isempty(w.frames) == 0
    assert(isempty(w.joints),'Frame Mode requires joints empty');
    assert(isfield(w,'njoints'),'Frame Mode requires explicit number of joints');
    assert(length(w.frames) == size(w.DH,1),'Number of frames same as DH params');
    assert(min(w.sensors) >= 1 && max(w.sensors) <= length(w.frames),'Sensors attached should reference valid frames'); 
    assert(size(w.DH,2) >= 6,'Frame Mode requires DH to be 6 columns: a alpha d theta rotoidal jointnumber');
    assert(max(double(w.DH(:,6))) <= w.njoints,'Frame Mode requires DH joint to be max to njoints');
else    
    assert(length(w.joints) == size(w.DH,1),'Number of joints same as DH params');    
    
    % now build the frames as in the frame mode
    if size(w.DH,2) < 5
        w.DH = [w.DH,ones(size(w.DH,1),1)];
    end
    if size(w.DH,2) < 6
        w.DH = [w.DH,(1:length(w.joints))'];
    end
    w.frames = w.joints;
    w.njoints = length(w.joints);
end
F = length(w.frames);
N = w.njoints;
M = length(w.sensors);

w.N = N;
w.M = M;
w.F = F;

for I=1:length(w.frames)
    if w.frames(I) == -1
        w.frames(I) = I-1;
    end
end

w.showmsg = M > 1 | N > 3;


allparams = [];

if w.multim0 == 0
    syms m0x m0y m0z real
    m00 = [m0x m0y m0z]';
    allparams = [allparams; m00];
else
    m00 = [0 0 0]';
end
w.m0 = m00;

w.g0 = [0 0 0]';
if w.gaxis > 0
    w.g0(w.gaxis) = -9.81;
end

[allparams,qs,qf,qp,x] = setupParams(allparams,w.DH,w);
if w.onlysym < 2
[hp,hpi] = setupJoints(qp,qf,qs,w);
else
    hp = [];
    hpi = [];
end
[allparams,ors] = setupSensorsRpi(allparams,w);
if w.onlysym < 2
[sp,spi] = setupSensors(qs,qf,qp,ors,hp,w);
else
    sp = [];
    spi=[];
end


% Aggregate all the output sensor variables: o ddx m
% TODO: make selectable
sensornames = {'o','do','ddx','m','g','Ti2p','Ti2w','dx'};
jointnames = {'o','do','ddx','m','g','Ti2p','Ti2w','dx'};
selectedjointoutputs = [1,3,4]; % omega acc mag

h = [];
if w.onlysym < 2
for J=1:M
    for K=w.selectedsensoroutputs
        h = [h; sp{J,K}];
    end
end
end

% for every row of h
if w.genJacob
    % Jv 
    Jv = cell(N*3);
    for J=1:N
        Jv{(J-1)*3+1} = qs{J,1};
        Jv{(J-1)*3+2} = qs{J,2};
        Jv{(J-1)*3+3} = qs{J,3};
    end
    
    Ja = [];
    for J=1:length(h)
        Jr = [];
        % and every variable interested
        for Q=1:length(Jv)
            yy = (diff(h(J),Jv{Q}));
            Jr = [Jr , yy];
        end
        Ja = [Ja ; Jr];
    end
else
    Ja = [];
    Jv = [];
end

% flatten
x = x(:);
allparams = allparams(:);
if w.onlysym > 0
    hf = [];
    Jf = [];
else
    if w.genManual
        [hf,Jf] = genfunctionsManual(w,qp,'sensors','h','J');
    else
        [hf,Jf] = genfunctions(w,h,Ja,{x, allparams},'sensors','h','J');
    end
end

w.symparams = allparams;

% Prepare output
ro = [];
ro.hsym = h; % all simbolic measures
ro.h = hf;   % measures as function
ro.Jsym = Ja; % Jacobian symbolic
ro.J = Jf;  % Jacobian as function
ro.symparams = allparams; % all parameters
ro.x = x; % input types
ro.outputs = sensornames(w.selectedsensoroutputs);
ro.x_desc = '(q_i,dq_i,ddq_i)';
ro.xsize = 3*N;
ro.h_desc = 'output of all sensors (x state, params)';
ro.params = w;
ro.qp = qp;
ro.qs = qs;
ro.frames = w.frames;
ro.sensors = w.sensors;
ro.sensorstx = w.sensorstx;
ro.hp = hp;
ro.sp = sp;

if w.genJoint
    % Aggregate all the output sensor variables: o ddx m
    %for J=1:N
    %    jh = [jh; hp(J,selectedjointoutputs)'];
    %end
    jh = [];
    for J=1:N
        for K=selectedjointoutputs
            jh = [jh; hp{J,K}];
        end
    end
    jhf = genfunctions(w,jh,[],{x, allparams},'joints','h','');
    ro.jhsym = jh;
    ro.jh = jhf;
    ro.joutputs = jointnames(selectedjointoutputs);
end

% Intermediate Joints generator
% h_i(fullstate_x,params,allintermediates)
if w.genInt
    ro.ii = genIntermediate(hpi,spi,qs,qf,qp,w);
end

if isempty(w.filename) == 0
    save([w.filename '.mat'],'ro');
end

end

function R = rotx(t)
ct = cos(t);
st = sin(t);
R = [
    1   0    0 0 ; 
    0   ct  -st 0; 
    0   st   ct 0;
    0 0 0 1;
];
end

function R = rotz(t)
ct = cos(t);
st = sin(t);
R = [
    ct  -st  0 0;
    st   ct  0 0;
    0    0   1 0;
      0 0 0 1;
    ];
end    
    
function R = roty(t)
ct = cos(t);
st = sin(t);
R = [
    ct  0   st 0;
    0   1   0  0;
   -st  0   ct 0;
   0 0 0 1;
   ];
end
   
% Propertries of Skew http://www8.tfe.umu.se/courses/elektro/RobotControl/Lecture06_5EL158.pdf
% S + S' = 0
% S(alpha a + beta b) = alpha S(a) + beta S(a)
% S(a) p = a cross p = - p cross a = - S(p) a
% R S(a) R' = S(R a)  for any SO(3) and any vector a
% R S(a) R' p = R (a cross R'p) = R (a cross p) = Ra cross Rb = Ra cross p = S(Ra) p
% x' S x = 0 for any S and x
%
% Derivative:
% d/dtheta R(theta) with R in SO(3)
% d/dtheta R(theta) = S(omega(t)) R(t)
%      of rotating frame wrt to the ?xed frame at time

function S = skew(v)
S = [  0   -v(3)  v(2)
      v(3)  0    -v(1)
     -v(2) v(1)   0];
end

function S = skew4(v)
S = [  0   -v(3)  v(2) 0;
      v(3)  0    -v(1) 0;
     -v(2) v(1)   0 0; 0 0 0 1];
end

function T = trans(v)
T = [1 0 0 v(1);0 1 0 v(2); 0 0 1 v(3); 0 0 0 1];
end


function [ho,Jo] = genfunctions(w,h,Ja,vars,text,hpre,Jpre)
%
% Generalized Generator of Code from h and the optional Jacobian
%
% w is the struct of options (onlysym,filename)
% h,Ja are the function and Jacobian symbolic
%
% vars are the input variable
% text is the description
% hpre Jpre are the prefix of the generated filename
% showmsg is for printing details


    showmsg = w.showmsg;

    if isfield(w,'filename') && isempty(w.filename) == 0
        if showmsg
            disp(['Generating File Version of h for ' text]);
            tic;
        end
        matlabFunction(h,'vars',vars,'file',[w.filename hpre '.m']);
        if showmsg
            toc
        end
        if w.genJacob & isempty(Ja) == 0
            if showmsg
                disp(['Generating File Version of J for ' text]);
                tic;
            end
            matlabFunction(Ja,'vars',vars,'file',[w.filename Jpre '.m']);
            if showmsg
                toc
            end
        end
        ho = [];
        Jo = [];
    else
        if showmsg
            disp(['Generating Numeric Version of h for ' text]);
            tic;
        end
        ho = matlabFunction(h,'vars',vars);
        if showmsg
            toc;
        end

        if w.genJacob & isempty(Ja) == 0
            if showmsg
                disp(['Generating Numeric Version of J for ' text]);
                tic;
            end
            Jo = matlabFunction(Ja,'vars',vars);
            if showmsg
                toc
            end
        else
            Jo = [];
        end
    end
end


function [ho,Jo] = genfunctionsManual(w,qp,text,hpre,Jpre)
%
% Generalized Generator of Code from h and the optional Jacobian
%
% w is the struct of options (onlysym,filename)
% h,Ja are the function and Jacobian symbolic
%
% vars are the input variable
% text is the description
% hpre Jpre are the prefix of the generated filename
% showmsg is for printing details


    showmsg = w.showmsg;

    if isfield(w,'filename') && isempty(w.filename) == 0
        if showmsg
            disp(['Generating File Version of h for ' text]);
            tic;
        end
        % TODO matlabFunction(h,'vars',vars,'file',[w.filename hpre '.m']);
        armmatlabFunction(w,qp,text,hpre,Jpre);
        if showmsg
            toc
        end
        if w.genJacob 
            error('not implemented');
            if showmsg
                disp(['Generating File Version of J for ' text]);
                tic;
            end
            matlabFunction(Ja,'vars',vars,'file',[w.filename Jpre '.m']);
            if showmsg
                toc
            end
        end
        ho = [];
        Jo = [];
    else
        error('Not supported in manual mode');
    end
end



function S = genvar(prefix,xsize,index)
% Generates a new variable with given vectorial size and prefix
% the index is appended to the name
    if nargin > 2
      prefix = sprintf([prefix '%d'],index);
    end              
    if xsize == 1
        S = sym(prefix,'real');
    else
        if length(xsize) == 1
          xsize = [xsize,1];
        end
        S = sym(sym(prefix,xsize),'real');        
    end
end

function [sp,spi] = setupSensors(qs,qf,qp,ors,hp,w)
%OLD convention: sensorparent was the parent of the attachment
% but this not works well with bifurcation because we need ca=pa+1
%NEW convention: sensorparent=ca and pa=parent
        
M = w.M;

z0 = [0 0 1]';

% o do ddx m g Tpi T0i dx
sp = (cell(M,8));

% o do ddx m g Tpi T0i dx 
spi = (cell(M,8));
        
% sensors
for I=1:M
    ca = w.sensors(I); % associated frame
    pa = w.frames(ca); % parent frame
    
    rt = qp{ca,5};
    cjoint = qp{ca,6};
    
    q = qs{cjoint,1};
    dq = qs{cjoint,2};
    ddq = qs{cjoint,3};    
           
    % Transformation is: Z Y X T using coefficients XYZ
    if w.sensorstx(I) == 1
        rs = ors{I,1};
        ts = ors{I,2};
        rx = rs(1);
        ry = rs(2);        
        % rotation comprises the attachment segment over parent       
        if rt==1
            rz = rs(3)+q;  
        else
            rz = rs(3);
        end
        S11 = [cos(rz) -sin(rz) 0 0;sin(rz) cos(rz) 0 0;0 0 1 0;0 0 0 1];
        S12 = [cos(ry) 0 sin(ry) 0; 0 1 0 0;-sin(ry) 0 cos(ry) 0; 0 0 0 1];
        S13 = [1 0 0 0;0 cos(rx) -sin(rx) 0;0 sin(rx) cos(rx) 0; 0 0 0 1];
        S14 = [1 0 0 ts(1); 0 1 0 ts(2); 0 0 1 ts(3); 0 0 0 1];
        S = S11*S12*S13*S14; 
        if rt==1
            S = trans([ 0 0 q])*S;
        end
        rpi = [ts(1),ts(2),ts(3)]'; %= Rp2i*S(1:3,4);
    else                
        if rt==1
            rz = q;
            S = [cos(rz) -sin(rz) 0 0;sin(rz) cos(rz) 0 0;0 0 1 0;0 0 0 1];
        else
            S = eye(4);
        end
        rpi = [0,0,0]';
    end
    
    Ti2p = S;
    Rp2i = S(1:3,1:3)';
    
   if pa == 0
       pahp = w.N+1;
   else
       pahp = pa;
   end
    op   = hp{pahp,1};
    dop  = hp{pahp,2};
    ddxp = hp{pahp,3};
    mp   = hp{pahp,4};
    gp   = hp{pahp,5};
    %Tp2pp = hp{pahp,6);
    Tp2w = hp{pahp,7};
    dxp  = hp{pahp,8};
    
    g = Rp2i*gp;
   
   if rt == 1   
        o = Rp2i*(z0*dq + op);
        do = Rp2i*(ddq*z0-dq*skew(z0)*op+dop);        
        ddx = Rp2i * ddxp - skew(rpi) * do + skew(o)*skew(o)*rpi;        
        dx = Rp2i*dxp - skew(rpi) * o;
   else
        o = Rp2i*op;
        do = Rp2i*dop;        
        ddx = Rp2i * (ddxp +ddq*z0) + 2*dq*skew(o)*Rp2i*z0 - skew(rpi)*do + skew(o)*skew(o)*rpi;
        dx = Rp2i*(dxp + dq*z0) - skew(rpi) * o;
   end
   ddx = ddx + g;
   Ti2w = Tp2w*Ti2p; % i to world sensors
   
    if w.multim0 
        m0 = ors{I,3}; % m0 global
        m = Ti2w(1:3,1:3)'*m0;      
    else
        m = Rp2i*mp;
    end

    % First 7 as same of hp
    sp{I,1} = o;
    sp{I,2} = do;
    sp{I,3} = ddx;
    sp{I,4} = m;
    sp{I,5} = g;
    sp{I,6} = Ti2p;
    sp{I,7} = Ti2w;
    sp{I,8} = dx;
    
    % Special Version for Intermediate Mode
    if w.genInt
        if pa > 0
            % using explicit parent from joint
            op = qf{pa,1};
            dop = qf{pa,2};
            ddxp = qf{pa,3};
            mp = qf{pa,4};
            gp = qf{pa,5};
            dxp = qf{pa,8};
            
            if rt == 1
                o = Rp2i*(z0*dq + op);
                do = Rp2i*(ddq*z0-dq*skew(z0)*op+dop);
                ddx = Rp2i * ddxp - skew(rpi) * do + skew(o)*skew(o)*rpi;                    
                dx = Rp2i*dxp - skew(rpi) * o;
            else
                o = Rp2i*op;
                do = Rp2i*dop;
                ddx = Rp2i * (ddxp +ddq*z0) + 2*dq*skew(o)*Rp2i*z0 - skew(rpi)*do + skew(o)*skew(o)*rpi;
                dx = Rp2i*(dxp + dq*z0) - skew(rpi) * o;
            end
            if w.multim0
                m0 = ors{I,3}; % m0 global
                m = Ti2w(1:3,1:3)'*m0;      
            else
                m = Rp2i*mp;                
            end
            g = Rp2i*gp;
            ddx = ddx + g;

            spi(I,1) = o;
            spi(I,2) = do;
            spi(I,3) = ddx;
            spi(I,4) = m;
            spi(I,5) = g;
            % TODO Ti2p
            % TODO Ti2w
            spi(I,8) = dx;
        else
            % otherwise just use itself
            spi(I,:) = sp(I,:); 
        end
    end
end
end


function [allparams,ors] = setupSensorsRpi(allparams,w)

    M = w.M;
    ors = cell(M,3);
    for I=1:M
        if w.sensorstx(I) == 1
            ors{I,1} = genvar('rs_pi',3,I);
            ors{I,2} = genvar('ts_pi',3,I);

            tmp = ors{I,1};
            allparams = [allparams; tmp(:)];
            tmp = ors{I,2};
            allparams = [allparams; tmp(:)];
        end

        if w.multim0
            ors{I,3} = genvar('m0_',3,I);
            allparams = [allparams; ors{I,3}];
        end

    end
end


function r = isinfx(x)
    if isinf(double(x))
        r = 1;
    else
        r = 0;
    end
end


function [allparams, qs,qf,qp,x] = setupParams(allparams,DH,w)

N = w.N;
F = w.F;

prefixes = {'pa%d','palpha%d','pd%d','ptheta%d'};

qs = (cell(N,3));
qf = (cell(F,(6*(w.genInt == 1)))); %oi doi ddxi mi gi dxi
qp = (cell(F,6)); % params of frames: a alpha d qoff joint

if w.useexplicit == 0
    x = sym(sym('x',[N*3,1]),'real'); 
else
    x = sym(zeros(N*3,1));
end

for I=1:N
    if w.useexplicit == 1      
        % generate the qi
        qs{I,1} = genvar('q',1,I);
        qs{I,2}= genvar('dq',1,I);
        qs{I,3} = genvar('ddq',1,I);
        x((I-1)*3+1) = qs{I,1};
        x((I-1)*3+2) = qs{I,2};
        x((I-1)*3+3) = qs{I,3};
    else
        % use the full vector x
        qs{I,1} = x((I-1)*3+1);
        qs{I,2} = x((I-1)*3+2);
        qs{I,3} = x((I-1)*3+3);
    end
end

for I=1:F    
    if w.genInt    
        % generate the intermediates
        qf{I,1} = genvar('oi',3,I);
        qf{I,2} = genvar('doi',3,I);
        qf{I,3} = genvar('ddxi',3,I);
        qf{I,4} = genvar('mi',3,I);
        qf{I,5} = genvar('gi',3,I);
        qf{I,8} = genvar('dxi',3,I);
    end
    
    % explore the DH assigning properly. Eventually the params are symbolic
    for K=1:4    
        s = DH(I,K);
        try
            sd = double(s);
            sdu = 1;
            ss = {};
        catch me
            ss = symvar(s);
            sd = 0;
            sdu = 0;
        end
        if sdu && isnan(sd)            
            si = sym(sprintf(prefixes{K},I),'real');
            s = si;
        elseif sdu && isinf(sd)
            si = sym(sprintf(prefixes{K},I),'real');
            s = sign(sd)*si;
        elseif isempty(ss) == 0
            % just add what is there
            for K=1:length(ss)
                allparams = [allparams; sym(ss{K})];
            end
            si=[];
        else
            % pure number
            si=[];
        end
        if isempty(si)==0
            allparams = [allparams; si];
        end
        qp{I,K} = s;
    end
    qp{I,5} = double(DH(I,5)); % type: 1=rotoidal 0=prismatic
    qp{I,6} = double(DH(I,6));
end


end


function [hp,hpi] = setupJoints(qp,qf,qs,w)

F = w.F;

% o do ddx m g 
hpi = (cell(F,8));
hp = (cell(F+1,8));
z0 = [0 0 1]';

% Extra
ie = F+1;
hp{ie,1} = sym([0,0,0]');
hp{ie,2}= sym([0,0,0]');
hp{ie,3} = sym([0,0,0]');
hp{ie,4} = sym(w.m0); % unique and used except for multim0
hp{ie,5} = sym(w.g0);
hp{ie,6} = sym(w.T0);
hp{ie,7} = sym(w.T0);
hp{ie,8} = sym([0,0,0]');

for I=1:F
    pa = w.frames(I); % parent frasym
    ca = I;           % current frame
    
    % extract the DH properties
    a = qp{ca,1};
    alpha = qp{ca,2};
    d = qp{ca,3};
    theta = qp{ca,4};
    rt = qp{ca,5}; % type: 1 is prismatic    
    cjoint = qp{ca,6};
    
    % extract the joint variable associated to the given frame
    q = qs{cjoint,1};
    dq = qs{cjoint,2};
    ddq = qs{cjoint,3};
    
    
    % build effective elements of the DH
    if rt == 1
        qreal = theta+q;
        dreal = d;
    else
        qreal = theta;
        dreal = d+q;
    end
    
    % build the DH as 4x4
    Ri2p44 = rotz(qreal)*rotx(alpha); 
    Rp2i = Ri2p44(1:3,1:3)';
    Ti2p = trans([ 0 0 dreal]) * rotz(qreal) * trans([a 0 0]) * rotx(alpha);
    rpi = [a, dreal*sin(alpha),dreal*cos(alpha)]';
    
    % adjust the parent index to the dummy
    if pa == 0
        pahp = ie;
    else
        pahp = pa;
    end
    
    % extract all that is needed
    op   = hp{pahp,1};
    dop  = hp{pahp,2};
    ddxp = hp{pahp,3};
    mp   = hp{pahp,4};
    gp   = hp{pahp,5};
    %Tp2pp = hp{pahp,6};
    Tp2w = hp{pahp,7};
    dxp  = hp{pahp,8};
    
    % Rotoidal
    if rt == 1
        o = Rp2i*(z0*dq + op);
        do = Rp2i*(ddq*z0-dq*skew(z0)*op+dop);
        dx = Rp2i * dxp - skew(rpi)*o;
        ddx = Rp2i * ddxp    - skew(rpi) * do + skew(o)*skew(o)*rpi;        
    else
        o = Rp2i*op;
        do = Rp2i*dop; 
        dx = Rp2i * (dxp + dq*z0) - skew(rpi)*o;
        ddx = Rp2i * (ddxp + ddq*z0) + 2*dq*skew(o)*Rp2i*z0   - skew(rpi)*do + skew(o)*skew(o) *rpi;
    end
    
    m = Rp2i*mp;
    g = Rp2i*gp;
    Ti2w = Tp2w*Ti2p;
       
    % Assign
    hp{I,1} = o;
    hp{I,2} = do;
    hp{I,3} = ddx;
    hp{I,4} = m;
    hp{I,5} = g;
    hp{I,6} = Ti2p;
    hp{I,7} = Ti2w;
    hp{I,8} = dx;
    
    if w.genInt 
        if pa ~= 0
            % using explicit parent instead of full expression
            op = qf{pa,1};
            dop = qf{pa,2};
            ddxp = qf{pa,3};
            mp = qf{pa,4};            
            gp = qf{pa,5};
            dxp = qf{pa,8};
            
            if rt == 1
                o = Rp2i*(z0*dq + op);
                do = Rp2i*(ddq*z0-dq*skew(z0)*op+dop);
                dx = Rp2i * dxp  - skew(rpi)*o;
                ddx = Rp2i * ddxp - skew(rpi) * do + skew(o)*skew(o)*rpi;        
            else
                o = Rp2i*op;
                do = Rp2i*dop; 
                dx = Rp2i * (dxp + dq*z0) - skew(rpi)*o;
                ddx = Rp2i * (ddxp + ddq*z0) + 2*dq*skew(o)*Rp2i*z0   - skew(rpi)*do + skew(o)*skew(o) *rpi;
            end
            m = Rp2i*mp;
            g = Rp2i*gp;
            hpi{I,1} = o;
            hpi{I,2} = do;
            hpi{I,3} = ddx;
            hpi{I,4} = m;
            hpi{I,5} = g;
            % Ti2pi not implemented
            % Ti2w not implemented
            hpi{I,8} = dx;
        else
            % same
            hpi(I,:) = hp(I,:); % same
        end        
    end
   
end




end


function ii = genIntermediate(spi,hpi,qs,qf,qp,w)
    M = w.M;
    N = w.N;
    F = w.F;
    showmsg  = w.showmsg;
    
    % the first input is the same, that is: x
    % the third input corresponds to the selected parent vars
    
    parentvars = 1:5; % all: o do ddx m g (optional TODO dx)
    joutputvars = 1:5; % all outputs: o do ddx m g (optional TODO dx)
    soutputvars = [1,3,4]; % just measured as usual: o ddx m
    
    % Extract the relevant variables from the qf
    zall = [];
    zallidx = cell(F,1);
    for I=1:F
        pre = length(zall);
        warning('fix');
        zall = [zall; qf(I,parentvars)'];
        K = length(zall)-pre+1;
        zallidx{I} = (pre+1):(pre+K-1);
    end
       

    % one per joint and one per sensors
    % note: no Jacobian, although it could be possible
    ii = [];
    ii.frames = cell(N,5); % function and then subset of 
    ii.sensors = cell(M,5);
    ii.zall = zall; % all outputs of joints that are needed

    % first iterate over frames
    for I=1:F
        % Build h using frame functions
        warning('fix');
        jh = hpi(I,joutputvars)';
        
        cajoint = qp{I,6};
        
        xlocal = qs{cajoint,1:3}';
        
        % build z from parentvars
        pa = w.frames(I);
        if pa ~= 0
            warning('fix');
            zlocal = qf(pa,parentvars)';
            zi = zallidx{pa};
        else
            zlocal = dummy;
            zi = [];
        end
        
        % this is for global function style => make it option
       % vars = {x, allparams, zall};
        vars = {xlocal, allparams, zlocal};
        jhf = genfunctions(w,jh,[],vars,'jointsint',sprintf('hi%d',I),'');
        ii.joints{I,1} = jh;
        ii.joints{I,2} = jhf;
        ii.joints{I,3} = []; % qinputs
        ii.joints{I,4} = zi; % input parent indices
        ii.joints{I,5} = zallidx{I}; % output z indices (in zall)
        
    end
       
    % then iterate over sensors
    for I=1:M
        % Build h using joint functions
        jh = spi(I,soutputvars)';
        
        
        pa = w.sensors(I); % parent joint
        if pa ~= 0
            warning('fix')
            zlocal = qf(pa,parentvars)';
            zi = zallidx{pa};
        else
            zi = [];
            zlocal= dummy;
        end
         cajoint = qp{ca,6};
         
        if cajoint > 0
            xlocal = qs(cajoint,1:3)';
        else
            xlocal = dummy;
        end
        
        % this is for global function style => make it option
        %vars = {x, allparams, zall};
        vars = {xlocal, allparams, zlocal};
            
        assert(isempty(jh) == 0,sprintf('Sensor jh is empty! %d',I));
        jhf = genfunctions(w,jh,[],vars,'sensorsint',sprintf('his%d',I),'',showmsg);
        ii.sensors{I,1} = jh;
        ii.sensors{I,2} = jhf;
        ii.sensors{I,3} = []; %qinputs
        ii.sensors{I,4} = zi; %zi inputs
        ii.sensors{I,5} = []; %zoutput indices
        
        persensor = length(jh);
        first = persensor*(I-1)+1;
        last = first+persensor-1;
        ii.sensors{I,4} = first:last;
    end
    ii.sensors_desc = 'per sensor: symbolic, functionhandle, zall input indices, zall output indices';
    ii.joints_desc = 'per joint: symbolic, functionhandle, zall input indices, h all output indices';
    ii.joutputs = '(o_i ddx_i m_i do_i g_i)';
    ii.soutputs = '(o_i ddx_i m_i)';
end

function r = cell2sym(c)
c = c(:);
r = sym(zeros(numel(c),1));
for I=1:length(c)
    r(I) = c{I};
end

end
