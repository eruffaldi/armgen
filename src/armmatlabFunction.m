function armmatlabFunction(w,qp,text,hpre,Jpre)
%
% Function Generator - To overcome matlab limitation
% 

sensortypes = armgenSensorTypes();
frametype = [];
frametype.rotoidal = 1;
frametype.prismatic = 0;


prefix = {
'function [ Z ] = hRow( x, params )'
'%codegen'
'Z0 = [0 0 1]'';'
'q = x(1:3:end);'
'dq = x(2:3:end);'
'ddq = x(3:3:end);'
};
if w.testsym 
    prefix{end+1} = sprintf('pi = sym(''pi'');');
end
prefix{end+1} = sprintf('T0 = %s;',mat2str(w.T0));
prefix{end+1} = sprintf('g0 = %s;',mat2str(w.g0));

prefix{end+1} = '% TODO: expand the parameters present in params related to the DH';

allparams = w.symparams;

for I=1:length(allparams)
    prefix{end+1} = sprintf('%s = params(%d);',char(allparams(I)),I);
end

for I=1:w.F
    si = I;
    
    DHa = char(qp{si,1});
    DHalpha = char(qp{si,2});
    DHd = char(qp{si,3});
    DHtheta = char(qp{si,4});
    rt = qp{si,5};
    qi = qp{si,6};
    pai = w.frames(si);

    prefix{end+1} = sprintf('\n%% Frame T%d as children of T%d using q%d ',I,pai,qi);
    if rt == frametype.rotoidal 
        part = {
            sprintf('ct = cos(q(%d)+%s); st = sin(q(%d)+%s);',qi,DHtheta,qi,DHtheta);
            sprintf('ca = cos(%s); sa = sin(%s);',DHalpha,DHalpha);
            sprintf('A%d = [ct -st*ca st*sa %s*ct; st ct*ca -ct*sa %s*st; 0 sa ca %s; 0 0 0 1];',I,DHa,DHa,DHd);
            sprintf('T%d = T%d*A%d;',I,pai,I);
       };
    else
        part = {
            sprintf('ct = cos(%s); st = sin(%s);',DHtheta,DHtheta);
            sprintf('ca = cos(%s); sa = sin(%s);',DHalpha,DHalpha);
            sprintf('A%d = [ct -st*ca st*sa %s*ct; st ct*ca -ct*sa %s*st; 0 sa ca q(%d)+%s; 0 0 0 1];',I,DHa,DHa,qi,DHd);
            sprintf('T%d = T%d*A%d;',I,pai,I);
        };
    end
    
    prefix = [prefix; part(:)];
    
    
        
        
    prefix{end+1} = sprintf('p = A%d(1:3,1:3)''*A%d(1:3,4);',I,I);
    prefix{end+1} = sprintf('Rpi = A%d(1:3,1:3)'';',I);
    parts = prepareSub(I,pai,qi,rt);
    prefix = [prefix; parts(:)];

end

for I=1:w.M
    parts = cell(9,1);
    si = w.sensors(I); % associated friend frame
    
    DHa = char(qp{si,1});
    DHalpha = char(qp{si,2});
    DHd = char(qp{si,3});
    DHtheta = char(qp{si,4});
    rt = qp{si,5};
    qi = qp{si,6};
    pai = w.frames(si);
    
    st = w.sensortype(I);
    
        II = I+1000;
    useoffset = 0;
    
    if useoffset == 0
        DHtheta = '0';
        DHd = '0';
    end
    
    if st ~= sensortypes.joint
            prefix{end+1} = sprintf('\n%% SH Sensor S%d as T%d children of T%d using q%d (sibling of T%d)',I,II,pai,qi,si);
            %parts{1} = sprintf('m%d = params(%d,1:3);',I,I+1);
            parts{1} = '';
            if rt == frametype.rotoidal 
                parts{2} = sprintf('zs = q(%d)+rs_pi%d3+%s;',qi,I,DHtheta);
            else
                parts{2} = sprintf('zs = rs_pi%d3;',I);
            end
            parts{3} = sprintf('ys = rs_pi%d2;',I);
            parts{4} = sprintf('xs = rs_pi%d1;',I);
            parts{5} = sprintf('Sa = [cos(zs) -sin(zs) 0 0;sin(zs) cos(zs) 0 0;0 0 1 0;0 0 0 1];');
            parts{6} = sprintf('Sb = [cos(ys) 0 sin(ys) 0; 0 1 0 0;-sin(ys) 0 cos(ys) 0; 0 0 0 1];');
            parts{7} = sprintf('Sc = [1 0 0 0;0 cos(xs) -sin(xs) 0;0 sin(xs) cos(xs) 0; 0 0 0 1];');
            parts{8} = sprintf('Sd = [1 0 0 ts_pi%d1; 0 1 0 ts_pi%d2; 0 0 1 ts_pi%d3; 0 0 0 1];',I,I,I);
            
            if rt == frametype.rotoidal 
                parts{9} = sprintf('A%d = Sa*Sb*Sc*Sd;',II);
            else
                parts{9} = sprintf('A%d = [1 0 0 0; 0 1 0 0; 0 0 1 q(%d)+%s; 0 0 0 1]*Sa*Sb*Sc*Sd;',II,qi,DHd);
            end
            prefix = [prefix; parts(:)];
        
            parts1 = {
                sprintf('p = A%d(1:3,1:3)''*A%d(1:3,4);',II,II);
                sprintf('% USE poses');
                sprintf('Rpi = A%d(1:3,1:3)'';',II);
                sprintf('T%d = T%d*A%d;',II,pai,II);
            };
            parts2 = prepareSub(II,pai,qi,rt);
            prefix = [prefix; parts1(:); parts2(:)];        
    end    
   
end

prefix{end+1} = '%sensoroutputs not used';
prefix{end+1} = 'gg = 9.81;';
prefix{end+1} = 'Z = [...';
for I=1:w.M
    II = I+1000;    
    si = w.sensors(I); % associated friend frame
    qi = qp{si,6};
    if w.sensortype(I) == sensortypes.inertial
        if w.multim0 == 1
            prefix{end+1} = sprintf('omega%d; (ddx%d+T%d(1:3,1:3)''*g0)/gg; T%d(1:3,1:3)''*[m0_%d1,m0_%d2,m0_%d3]'';',II,II,II,II,I,I,I);
        else
            prefix{end+1} = sprintf('omega%d; (ddx%d+T%d(1:3,1:3)''*g0)/gg; T%d(1:3,1:3)''*[m0x,m0y,0z]'';',II,II,II,II);
        end    
    elseif w.sensortype(I) == sensortypes.joint
        prefix{end+1} = sprintf('q(%d);',qi);
    elseif w.sensortype(I) == sensortypes.globalpos
        prefix{end+1} = sprintf('T%d(1:3,4)',II);
    elseif w.sensortype(I) == sensortypes.globalz
        prefix{end+1} = sprintf('T%d(1:3,3)',II);
    elseif w.sensortype(I) == sensortypes.globalposAndz
        % special combination of position and z
        prefix{end+1} = sprintf('T%d(1:3,3);T%d(1:3,4)',II,II);
    end
end
prefix{end+1} = '];';

f = fopen(sprintf('%s%s.m',hpre,text),'w');
for I=1:length(prefix)
    fprintf(f,'%s\n',prefix{I});
end
fclose(f);

% then
if isempty(Jpre)
    % TODO: jacobian of skew matrix
    % Given:
    % Jdqj Jddqj Jqj
    
    % Rotoidal
    % oi = Rpi(op + dqj*z0)
    % Joi = JRpi (op + dqj*z0) + Rpi(Jop + Jdqj*z0) = JRpi Rpi' oi + Rpi
    % (Jop + Jdqj*z0);
    %
    % doi = Rpi(dop - dqj*S(z0)*oi + ddqj*z0)
    % Jdoi = JRpi (dop - dqj*S(z0)*oi + ddqj*z0)  + Rpi(Jdop - Jdqj*S(z0)*oi - dqj*S(z0)*Joi + Jddqj*z0)
    %  = JRpi Rpi' doi + ...
    %
    % ddxi = Rpi*ddxp - S(p) doi - S(p) S(oi) oi 
    % Jddxi = JRpi*Rpi'ddxi - S(p) Jdoi - S(p) JS(oi) oi - S(p) S(oi) Joi
    
    % Prismatic
    % oi = Rpi op
    % Joi = JRpi Rpi' oi + Rpi Jop;
    %
    % doi = Rpi(dop - dqj*S(z0)*oi + ddqj*z0)
    % Jdoi = JRpi Rpi' doi - Rpi Jdop;
    %
    % ddxi = Rpi*(ddxp+ddq*z0) - S(p) doi - S(p) S(oi) oi + 2 dqj S(o) Rpi z0
    % Jddxi =  ... complex really!
    
    return;
end

% generate the Jacobians

function out = join(glue, strs)
strs = strs(:)';
strs(2,:) = {glue};
strs = strs(:)';
strs(end) = [];
out = char(1, strs{:});

function prefix = prepareSub(I,pai,qi,rt)
    prefix = {};
    if rt == 1  % rotoidal
        if pai > 0
            prefix{end+1} = sprintf('omega%d = Rpi*(omega%d + dq(%d)*Z0);',I,pai,qi);
            prefix{end+1} = sprintf('domega%d = Rpi*(domega%d + dq(%d)*cross2(omega%d,Z0)+ddq(%d)*Z0);',I,pai,qi,pai,qi);
            prefix{end+1} = sprintf('dx%d = Rpi*dx%d + cross2(omega%d,p);',I,pai,I);
            prefix{end+1} = sprintf('ddx%d = Rpi*ddx%d + cross2(domega%d,p)+cross2(omega%d,cross(omega%d,p));',I,pai,I,I,I);
        else
            prefix{end+1} = sprintf('omega%d = Rpi*(dq(%d)*Z0);',I,qi);
            prefix{end+1} = sprintf('domega%d = ddq(%d)*Rpi*Z0;',I,qi);
            prefix{end+1} = sprintf('dx%d = cross2(omega%d,p);',I,I);
            prefix{end+1} = sprintf('ddx%d = cross2(domega%d,p)+cross2(omega%d,cross2(omega%d,p));',I,I,I,I);
        end
    else
        if pai > 0
            prefix{end+1} = sprintf('omega%d = Rpi*omega%d;',I,pai);
            prefix{end+1} = sprintf('domega%d = Rpi*domega%d;',I,pai);
            prefix{end+1} = sprintf('dx%d = Rpi*(dx%d +dq(%d)*Z0) + cross2(omega%d,p);',I,pai,qi,I);
 %       ddx = Rp2i * (ddxp +ddq*z0) + 2*dq*skew(o)*Rp2i*z0 - skew(rpi)*do + skew(o)*skew(o)*rpi;
            prefix{end+1} = sprintf('ddx%d = Rpi*(ddx%d+ddq(%d)*Z0)+2*dq(%d)*cross(omega%d,Rpi*Z0) + cross2(domega%d,p)+cross2(omega%d,cross2(omega%d,p));',I,pai,qi,qi,I,I,I,I);
        else
            prefix{end+1} = sprintf('omega%d = [0,0,0]'';',I);
            prefix{end+1} = sprintf('domega%d = [0,0,0]'';',I);            
            prefix{end+1} = sprintf('dx%d = Rpi*dq(%d)*Z0;',I,qi);
            prefix{end+1} = sprintf('ddx%d = ddq(%d)*Rpi*Z0;',I,qi);
          
        end        
    end