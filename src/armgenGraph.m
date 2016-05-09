function armgenGraph(w,qp,usedot)


if nargin < 3
usedot = 0;
end

% build the graph for graph output
%
% use 

if usedot
    c1 = '_';
    c2 = ' ';
else
    c1 = '';
    c2 = ' ';
end



A = zeros(1+w.F+w.M);
off1 = 1;
off2 = w.F+off1;

nl = cell(length(A),1);
al = cell(size(A));

if usedot
nl{1} = '$root$';
else
nl{1} = 'root';
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
    
    %TODO shape for prismatic
    
    if pai > 0
        A(I+off1,pai+off1) = 1;
        al{I+off1,pai+off1} = '';
    else
        A(I+off1,off1) = 1;
        al{I+off1,off1} = '';
    end
    if rt == 0
        shape = '",shape="rect';
    else
        shape = '';
    end
    if usedot
        nl{I+off1} = sprintf('T_{%d} q_{%d}%s',I,qi,shape);
    else
        nl{I+off1} = sprintf('T%d q%d%s',I,qi,shape);
    end
end



for I=1:w.M
    si = w.sensors(I); % associated friend frame    
    DHa = char(qp{si,1});
    DHd = char(qp{si,3});
    qi = qp{si,6};
    pai = w.frames(si);
    
    rt = qp{si,5};
    st = w.sensortype(I);
    
    %TODO shape for sensor type
    
    
    if pai > 0
        A(I+off2,pai+off1) = 1;
        al{I+off2,pai+off1} = '';
    else
        A(I+off2,off1) = 1;
        al{I+off2,off1} = '';
    end
%    A(I+off2,si+off1) = 1;
%    A(si+off1,I+off2) = 1;
    if rt == 0
        shape = '",shape="rect';
    else
        shape = '';
    end
    style ='",style="filled';
    if usedot
        nl{I+off2} = sprintf('S_{%d} q_{%d}%s%s',I,qi,shape,style);
    else
        nl{I+off2} = sprintf('S%d q%d%s%s',I,qi,shape,style);
    end
end

A
nl


graph_to_dot(A','filename','out.dot','node_label',nl); %'arc_label',al);
