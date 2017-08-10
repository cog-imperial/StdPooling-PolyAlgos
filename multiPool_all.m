function [] = multiPool_all( )

clear;
tic;
%options = optimoptions('linprog','Algorithm','interior-point','Display','off'); %'iter'
%{ 
inputVectConcs = [3   1    4    2   1.5];
inputVectCosts = [10   8   10    7   11];
inputPoolConns = {
                    [1 2 3];      % inputs to pool 1
                    [1 2 5];      % inputs to pool 2
                    [1 2 3 4];    % inputs to pool 3
                 };
outputMaxConcVect = [3 2.5];
outputProfitVect  = [14 15];
poolOutputConns = {
                    [1 2];      % pools to output 1
                    [2 3];      % pools to output 2
                  };
outputDemandsU    = [100, 200];

directVectConcs  = [2 6 2];
directVectCosts  = [10 3 8];
directOutputConns ={   
                    [1 2];    % directs to output 1
                    [1 2];    % directs to output 2
                   };
pStep = 0.05;
%}

%%{ 
% Original parameters
%{
inputVectConcs = [3   1    4    2   1.5];
inputVectCosts = [6   13   10   9   11];
inputPoolConns = {
                    [1 2 3];    % inputs to pool 1
                    [1 2 3 4];  % inputs to pool 2
                    [1 2 5];    % inputs to pool 3
                 };
outputMaxConcVect = [2.5];
outputProfitVect  = [15];
poolOutputConns = {
                    [1 2];      % pools to output 1
                  };
outputDemandsU    = [200];

directVectConcs  = [2 6 4];
directVectCosts  = [10 3 12];
directOutputConns ={   
                    [1 2 3];      % directs to output 1
                   };
pStep = 0.05;
%}

%{
Ci = [3   1];
gi = [6   13];
A_L = [0 0 ]; % feed availab inputs*
A_U = 100*[100 100];
inputPoolConns = {
                    [1 2];      % inputs to pool 1
                    %[1 2 3];    % inputs to pool 2
                 };
P_L = [1, 1]; %P_L = [0 0]; *
P_U = [2, 2];
outputProfitVect  = [15, 15];
poolOutputConns = {
                    [1];      % pools to output 1
                    [1];      % pools to output 2
                  };
D_L    = [0, 0]; %*
D_U    = [200, 150];

Cj  = [1.5  2.5];
gj  = [12   3];
%Cj  = [1.5  5 5.5];
%gj  = [10  5 4.6];
A_Lj = [0 0]; % feed availab directs*
A_Uj = 100*[100 100];
directOutputConns ={   
                    [1];    % directs to output 1
                    [2];    % directs to output 2
                   };
S_L = [0]; %lower pool capacities *
S_U = [350*100]; %upper pool capacities *

%S_L = [0, 0]; %lower pool capacities
%S_U = [350,350]; %upper pool capacities
pStep = 0.01;
%}

%{
Ci = [3   1];
gi = [6   16];
A_L = [0 0 ]; % feed availab inputs*
A_U = 100*[100 100];
inputPoolConns = {
                    [1 2];      % inputs to pool 1
                    %[1 2 3];    % inputs to pool 2
                 };
P_L = [0, 0]; %P_L = [0 0]; *
P_U = [2.5, 1.5];
outputProfitVect  = [9, 15];
poolOutputConns = {
                    [1];      % pools to output 1
                    [1];      % pools to output 2
                  };
D_L    = [0, 0]; %*
D_U    = [100, 200];

Cj  = [2];
gj  = [10];
%Cj  = [1.5  5 5.5];
%gj  = [10  5 4.6];
A_Lj = [0]; % feed availab directs*
A_Uj = 100*[100];
directOutputConns ={   
                    [1];    % directs to output 1
                    [1];    % directs to output 2
                   };
S_L = [0]; %lower pool capacities *
S_U = [350*100]; %upper pool capacities *

%S_L = [0, 0]; %lower pool capacities
%S_U = [350,350]; %upper pool capacities
pStep = 0.01;
%}

%{
Ci = [0.1    0.25  0.3  0.35 0.45   0.5];
gi = [22     10      4    4    4    4];
A_L = [0 0 0 0 0 0]; % feed availab inputs*
A_U = 100*[100 100 100 100 100 100];
inputPoolConns = {
                    [3 4];      % inputs to pool 1
                    [1 2];    % inputs to pool 2
                    [5 6];    % inputs to pool 3
                 };
P_L = [0.15, 0.3]; %P_L = [0 0]; *
P_U = [0.20, 0.4];
outputProfitVect  = [20, 20];
poolOutputConns = {
                    [1,2];      % pools to output 1
                    [2,3];      % pools to output 2
                  };
D_L    = [100, 100]; %*
D_U    = [100, 100];

Cj  = [1000];
gj  = [100000];
%Cj  = [1.5  5 5.5];
%gj  = [10  5 4.6];
A_Lj = [0]; % feed availab directs*
A_Uj = [0];
directOutputConns ={   
                    [1];    % directs to output 1
                    [1];    % directs to output 2
                   };
S_L = [0, 0, 0]; %lower pool capacities *
S_U = [350*100, 350*100, 350*100]; %upper pool capacities *
pStep = 0.01;
%}

%
Ci = [3   1];
gi = [6   13];
A_L = [0 0]; % feed availab inputs*
A_U = 100*[100 100];
inputPoolConns = {
                    [1 2];      % inputs to pool 1
                 };
P_L = [1, 1]; %P_L = [0 0]; *
P_U = [2, 2];
outputProfitVect  = [15, 15];
poolOutputConns = {
                    [1];      % pools to output 1
                    [1];      % pools to output 2
                  };
D_L    = [0, 0]; %*
D_U    = [200, 150];

Cj  = [1.5 2.5];
gj  = [12 3];
%Cj  = [1.5  5 5.5];
%gj  = [10  5 4.6];
A_Lj = [0 0]; % feed availab directs*
A_Uj = 100*[100 100];
directOutputConns ={   
                    [1];    % directs to output 1
                    [2];    % directs to output 2
                   };
S_L = [0]; %lower pool capacities *
S_U = [350*100]; %upper pool capacities *
pStep = 0.01;
%}

%{
Ci = [3   1    4    2   1.5];
gi = [6   13   10   9   11];
A_L = [0 0 0 0 0]; % feed availab inputs*
A_U = 100*[100 100 100 100 100];
inputPoolConns = {
                    [1 2 3 4 5];      % inputs to pool 1
                    %[1 2 3];    % inputs to pool 2
                 };
P_L = [1.5]; %P_L = [0 0]; *
P_U = [3];
outputProfitVect  = [14];
poolOutputConns = {
                    [1];      % pools to output 1
                    %[1];      % pools to output 2
                  };
D_L    = [0]; %*
D_U    = [250];

Cj  = [1.5  3 5];
gj  = [10 11 9.5];
%Cj  = [1.5  5 5.5];
%gj  = [10  5 4.6];
A_Lj = [0 0 0]; % feed availab directs*
A_Uj = 100*[100 100 100];
directOutputConns ={   
                    [2 3];    % directs to output 1
                    %[1 2];    % directs to output 2
                   };
S_L = [0]; %lower pool capacities *
S_U = [350*100]; %upper pool capacities *

%S_L = [0, 0]; %lower pool capacities
%S_U = [350,350]; %upper pool capacities
pStep = 0.01;
%}

%{
Ci = [3   1    4    2   1.5];
gi = [6   13   10   9   11];
A_L = [0 0 0 0 0]; % feed availab inputs*
A_U = 100*[100 100 100 100 100];
inputPoolConns = {
                    [1 2 3 4 5];      % inputs to pool 1
                    %[1 2 3];    % inputs to pool 2
                 };
P_L = [2, 1.5]; %P_L = [0 0]; *
P_U = [3.5, 3];
outputProfitVect  = [15, 14];
poolOutputConns = {
                    [1];      % pools to output 1
                    [1];      % pools to output 2
                  };
D_L    = [0, 0]; %*
D_U    = [200, 250];

Cj  = [1.5  3 5];
gj  = [10 11 9.5];
%Cj  = [1.5  5 5.5];
%gj  = [10  5 4.6];
A_Lj = [0 0 0]; % feed availab directs*
A_Uj = 100*[100 100 100];
directOutputConns ={   
                    [1 2 3];    % directs to output 1
                    [2 3];    % directs to output 2
                   };
S_L = [0]; %lower pool capacities *
S_U = [350*100]; %upper pool capacities *

%S_L = [0, 0]; %lower pool capacities
%S_U = [350,350]; %upper pool capacities
pStep = 0.01;
%}

%{ 
Ci = [3   1    4    2   1.5];
gi = [6   13   10   9   11];
A_L = [0 0 0 0 0]; % feed availab inputs*
A_U = 0.5*[100 100 100 100 100];
inputPoolConns = {
                    [1 2 3 4 5];      % inputs to pool 1
                    %[1 2 3];    % inputs to pool 2
                 };
P_L = [1.5 2.5]; %P_L = [0 0]; *
P_U = [2 4];
outputProfitVect  = [15 15];
poolOutputConns = {
                    [1];      % pools to output 1
                    [1];      % pools to output 2
                  };
D_L    = [150, 100]; %*
D_U    = [200, 150];

Cj  = [1.5  3];
gj  = [10 11];
A_Lj = [30 30]; % feed availab directs*
A_Uj = 1*[100 100];
directOutputConns ={   
                    [1 2];    % directs to output 1
                    [1 2];    % directs to output 2
                   };
S_L = [20]; %lower pool capacities *
S_U = [350]; %upper pool capacities *

%S_L = [0, 0]; %lower pool capacities
%S_U = [350,350]; %upper pool capacities
pStep = 0.01;
%}

excelExport = 0;
%%%%%%%%%%%%%%%%%%
nbOutputs = length(P_U);
allResults = [];
all2SamePools = all(cellfun('length',poolOutputConns)==2)&& nnz(diff(cell2mat(poolOutputConns),1,1))==0;
all1Pool = all(cellfun('length',poolOutputConns)==1);
% indices of how many x, y, z vars are present in obj func for an output
%fNbVars = [];
%%%%%%%%%% For each output
%%{
for o=1:nbOutputs
    %%%%%% Subproblems (1 subproblem/output node, with all corresponding
    %%%%%% x,y,z in this order as variables)
    
    %Determine the pools going into the output and
    %relevant costs and conc% for inputs, directs for this output
    poolsConnected = poolOutputConns{o};
    nbPoolsAtOutput = length(poolsConnected);
    inputVectConcPerPool = arrayfun(@(pool) Ci(inputPoolConns{pool}), poolsConnected','UniformOutput', false);
    inputVectCostPerPool = arrayfun(@(pool) gi(inputPoolConns{pool}), poolsConnected','UniformOutput', false);
    directVectConc = Cj(directOutputConns{o});
    directVectCost = gj(directOutputConns{o});
    
    inputCostsPerOutput = horzcat(inputVectCostPerPool{:});
    inputConcsPerOutput = horzcat(inputVectConcPerPool{:});
    nbInputsAtOutput = length(inputConcsPerOutput);
    nbDirectsAtOutput = length(directVectConc);
    inputOnesVect = ones(1,nbInputsAtOutput);
    poolsOnesVect = ones(1,nbPoolsAtOutput);
    directOnesVect = ones(1,nbDirectsAtOutput);
    totalNbNodes = nbInputsAtOutput + nbPoolsAtOutput + nbDirectsAtOutput;
    
    %objective function vector
    f = [inputCostsPerOutput, -outputProfitVect(o)*poolsOnesVect, -outputProfitVect(o)+directVectCost];
    %fNbVars(o,:) = [nbInputsAtOutput, nbPoolsAtOutput, nbDirectsAtOutput];
    %%% inequalities [ upper+lower product demand constraint; upper+lower product quality contraint; upper+lower pool capacities]
    Ao = [inputOnesVect*0, poolsOnesVect, directOnesVect;
        inputOnesVect*0, -poolsOnesVect, -directOnesVect;
        inputConcsPerOutput, -P_U(o)*poolsOnesVect, directVectConc-P_U(o);
        -inputConcsPerOutput, P_L(o)*poolsOnesVect, -directVectConc+P_L(o);
        zeros(nbPoolsAtOutput*2,nbInputsAtOutput), [diag(ones(1,nbPoolsAtOutput)); -diag(ones(1,nbPoolsAtOutput))], zeros(nbPoolsAtOutput*2,nbDirectsAtOutput)*0;];
        %-diag(ones(1,totalNbNodes))];
    bo = [D_U(o) -D_L(o) zeros(1,2) S_U(poolsConnected) -S_L(poolsConnected)]; %zeros(1,2+totalNbNodes)];
    if o==1
       bo = [D_U(o) -D_L(o) zeros(1,2) S_U(poolsConnected) -S_L(poolsConnected)*0]; 
    end
    %%% equalities [ material balances/ pool, quality balances/ pool ]
    inputOnesVectPerPool = cell2mat(arrayfun(@(p)...
        [zeros(1,(p>1)*length(horzcat(inputVectConcPerPool{1:p-1}))),...
        ones(1,length(inputVectConcPerPool{p})),...
        zeros(1,(p<nbPoolsAtOutput)*length(horzcat(inputVectConcPerPool{p+1:end})))]',...
        1:nbPoolsAtOutput,'UniformOutput', false))';
    inputConcVectPerPool = cell2mat(arrayfun(@(p)...
        [zeros(1,(p>1)*length(horzcat(inputVectConcPerPool{1:p-1}))),...
        inputVectConcPerPool{p},...
        zeros(1,(p<nbPoolsAtOutput)*length(horzcat(inputVectConcPerPool{p+1:end})))]',...
        1:nbPoolsAtOutput,'UniformOutput', false))';
    Aeq = [inputOnesVectPerPool, -diag(ones(1,nbPoolsAtOutput)), zeros(nbPoolsAtOutput, nbDirectsAtOutput);
        inputConcVectPerPool, zeros(nbPoolsAtOutput, nbPoolsAtOutput + nbDirectsAtOutput)];
    pAeq = [zeros(nbPoolsAtOutput,totalNbNodes);
        inputOnesVectPerPool, zeros(nbPoolsAtOutput, nbPoolsAtOutput + nbDirectsAtOutput)];
    beq = zeros(1,2*nbPoolsAtOutput);
    
    % building iteration mesh over p's
    origPoolsIters= arrayfun(@(L,B) L:pStep:B,...
        cellfun(@min,inputVectConcPerPool),...
        cellfun(@max,inputVectConcPerPool),...
        'UniformOutput',false);
    
    origPoolsIters = combvec(origPoolsIters{:})';
    iters = size(origPoolsIters,1);
    poolsIters = [origPoolsIters*0 origPoolsIters];
    poolsIters = mat2cell(poolsIters,ones(1,size(poolsIters,1)),size(poolsIters,2));
    %%%% calling linear solver on poolsIters mesh of p values for all pools
    %%%% p values make a difference in the Aeq equality matrix
    res = zeros(size(origPoolsIters,1),1+totalNbNodes);
    b = [bo,beq]';
    e = [ones(1,size(Ao,1)), 3*ones(1,size(Aeq,1))]';
    lp = mxlpsolve('make_lp', length(e), totalNbNodes);
    mxlpsolve('set_verbose', lp, 3);    
    mxlpsolve('set_rh_vec', lp, b);
    mxlpsolve('set_obj_fn', lp, -f);
    mxlpsolve('set_maxim', lp);
    for i = 1:totalNbNodes
        mxlpsolve('set_lowbo', lp, i, 0);
    end
    for i = 1:length(e)
        mxlpsolve('set_constr_type', lp, i, e(i));
    end  
    for i=1:iters
        mxlpsolve('set_mat', lp, [Ao;Aeq - bsxfun(@times,poolsIters{i}',pAeq)]);
        mxlpsolve('solve', lp);
        [res(i,1), res(i,2:end), ~] = mxlpsolve('get_solution', lp);
    end
    mxlpsolve('delete_lp', lp);
    
    if nbPoolsAtOutput >= 2
        %plotDataTable(excelExport,[origPoolsIters res],o,sprintf('Numeric values for output %d',o),sprintf('fval(%d)',o),poolsConnected,inputPoolConns,directOutputConns,poolOutputConns);
    end
    if nbPoolsAtOutput == 2
        %plot3D(sprintf('Mesh for output %d', o),[origPoolsIters, res(:,1)],pStep,poolsConnected,sprintf('Obj function for output %d', o));
    elseif (nbPoolsAtOutput == 1)   
        %plotDataTable(excelExport,[origPoolsIters res],o,sprintf('Numeric values for output %d',o),sprintf('fval(%d)',o),poolsConnected,inputPoolConns,directOutputConns,poolOutputConns);
        if (~all1Pool)
            %plot 2D if 1 pool
            figure2d = figure('name','Results','units','normalized','outerposition',[0 0.15 1 0.85]);
            axes1 = axes('Parent',figure2d,'YGrid','on','XGrid','on',...
                'Position',[0.02 0.079 0.70 0.84]);
            box(axes1,'on');
            hold(axes1,'on');
            plot(origPoolsIters,res(:,1),'LineWidth',2,'Parent',axes1);
            xlabel(sprintf('Pool p%d', poolsConnected(1)),'FontSize',11);
            zlabel(sprintf('Obj function for output %d', o),'FontSize',11);
        else
            allResults(:,o)=res(:,1);
        end
    end
end
%}

%%%%%% Original Full Problem (all x's y's and z's in this order as variables)
%%%%%%%%%%%%%
directsConnected = unique(cat(2,directOutputConns{:}));
nbDirectsUnique = length(directsConnected);
poolsConnected2 = unique(reshape(cell2mat(poolOutputConns),1,[]));
nbPoolsTotal = cellfun(@length, poolOutputConns);
inputVectConcPerPool2 = arrayfun(@(pool) Ci(inputPoolConns{pool}), poolsConnected2','UniformOutput', false);
inputVectCostPerPool2 = arrayfun(@(pool) gi(inputPoolConns{pool}), poolsConnected2','UniformOutput', false);
directVectConc = arrayfun(@(o) Cj(directOutputConns{o}),(1:nbOutputs)','UniformOutput',false);
directVectCost = arrayfun(@(o) gj(directOutputConns{o}),(1:nbOutputs)','UniformOutput',false);

inputCostsPerOutput = horzcat(inputVectCostPerPool2{:});
inputConcsPerOutput = horzcat(inputVectConcPerPool2{:});
nbInputsTotal = length(inputConcsPerOutput);
nbDirectsTotal = sum(cellfun(@length, directVectConc));
poolsOnesVect = arrayfun(@(p) ones(1,p), nbPoolsTotal,'UniformOutput',false);
totalNbNodes = nbInputsTotal + sum(nbPoolsTotal) + nbDirectsTotal;
objYpart = cell2mat(arrayfun(@(o) outputProfitVect(o)*poolsOnesVect{o},1:nbOutputs,'UniformOutput',false));
objZpart = cell2mat(arrayfun(@(o) outputProfitVect(o)*ones(1,length(directOutputConns{o}))-directVectCost{o}, 1:nbOutputs, 'UniformOutput',false));
%objective function vector
fOrig = [inputCostsPerOutput, -objYpart, -objZpart];

allpoolsOutputsConns = cell2mat(poolOutputConns');
YonesEq = -cell2mat(arrayfun(@(pool) allpoolsOutputsConns==pool, poolsConnected2', 'UniformOutput',false));

%%% inequalities [ upper+lower product demand constraint; upper+lower product quality contraint;  
%                  upper+lower pool capacities; upper+lower feed availab for inputs+directs]
Yones = blkdiag(poolOutputConns{:});   Yones(Yones>0) = 1;
Zones = blkdiag(directOutputConns{:}); Zones(Zones>0) = 1;

inputRow = cat(2,inputPoolConns{poolsConnected2}); 
inputRow = cell2mat(arrayfun(@(i) (inputRow==i)*1, unique(inputRow)','UniformOutput', false));
directRow = cat(2,directOutputConns{:}); 
directRow = cell2mat(arrayfun(@(i) (directRow==i)*1, unique(directRow)','UniformOutput', false));

AOrig = [zeros(nbOutputs, nbInputsTotal), Yones, Zones;
         zeros(nbOutputs, nbInputsTotal), -Yones, -Zones;
         zeros(nbOutputs, nbInputsTotal), -diag(P_U)*Yones, blkdiag(directVectConc{:})-diag(P_U)*Zones;
         zeros(nbOutputs, nbInputsTotal), diag(P_L)*Yones, -blkdiag(directVectConc{:})+diag(P_L)*Zones;
         zeros(2*length(poolsConnected2), nbInputsTotal), [-YonesEq; YonesEq], zeros(2*length(poolsConnected2),nbDirectsTotal);
         inputRow zeros(nbInputsTotal, nbDirectsTotal+sum(nbPoolsTotal));
         -inputRow zeros(nbInputsTotal, nbDirectsTotal+sum(nbPoolsTotal));
         zeros(nbDirectsUnique, nbInputsTotal+sum(nbPoolsTotal)) directRow;
         zeros(nbDirectsUnique, nbInputsTotal+sum(nbPoolsTotal)) -directRow;
         ];
bOrig = [D_U, -D_L, zeros(1,2*nbOutputs), S_U, -S_L, A_U, -A_L, A_Uj(directsConnected), -A_Lj(directsConnected)];

%%% equalities [ material balances/ pool, quality balances/ pool for all pools]
nbUniquePools = length(poolsConnected2);
XonesEq = blkdiag(inputPoolConns{poolsConnected2});   XonesEq(XonesEq>0) = 1;
XEq     = blkdiag(inputVectConcPerPool2{:});
ZzerosEq = zeros(nbUniquePools, nbDirectsTotal);
AOrigEq = [XonesEq, YonesEq,   ZzerosEq;
    XEq,     YonesEq*0, ZzerosEq ];
pAOrigEq = [zeros(nbUniquePools, totalNbNodes);
    XonesEq, YonesEq*0, ZzerosEq ];
bOrigEq = zeros(1,2*nbUniquePools);

% building iteration mesh over pools
%lowerBounds = zeros(1,totalNbNodes)';
origPoolsIters= arrayfun(@(L,B) L:pStep:B,...
    cellfun(@min,inputVectConcPerPool2),...
    cellfun(@max,inputVectConcPerPool2),...
    'UniformOutput',false);
origPoolsIters = combvec(origPoolsIters{:})';
iters = size(origPoolsIters,1);
%poolConns = cellfun(@(v) {ones(1,length(v))},poolOutputConns);
poolsIters = [origPoolsIters*0 origPoolsIters];
%%%% calling linear solver on poolsIters mesh of p values for all pools
%%%% p values make a difference in the AOrigEq equality matrix
%%%% pIneq values make a difference in the AOrig equality matrix
res = zeros(size(origPoolsIters,1),1+totalNbNodes);
%%{
var = mat2cell(origPoolsIters(:,horzcat(poolOutputConns{:})),ones(iters,1),nbPoolsTotal');
[p2,m2] = cellfun(@size,var(1,:));
p1 = [0, cumsum(p2)];
m1 = [0, cumsum(m2)];
indices = cell(2,nbOutputs);
for k=1:nbOutputs
    indices(1:2,k)= {p1(k)+1:p1(k+1);m1(k)+1:m1(k+1)};
end
y = zeros(p1(end),m1(end)); %Preallocate

b = [bOrig,bOrigEq]';
e = [ones(1,size(AOrig,1)), 3*ones(1,size(AOrigEq,1))]';
lp = mxlpsolve('make_lp', length(e), totalNbNodes);
mxlpsolve('set_verbose', lp, 3);
mxlpsolve('set_rh_vec',  lp, b);
mxlpsolve('set_obj_fn',  lp, -fOrig);
mxlpsolve('set_maxim',   lp);
for i = 1:totalNbNodes
    mxlpsolve('set_lowbo', lp, i, 0);
end
for i = 1:length(e)
    mxlpsolve('set_constr_type', lp, i, e(i));
end
xInd  = 2*nbOutputs+1 : 4*nbOutputs;
xInd2 = 4*nbOutputs+2*length(poolsConnected2)+2*nbInputsTotal+2*nbDirectsUnique+1;
yInd  = nbInputsTotal+1 : nbInputsTotal+sum(nbPoolsTotal);
AOrigCopy = AOrig(xInd,yInd);
Ap = [AOrig;AOrigEq];
for i=1:iters
    for k=1:nbOutputs
        y(indices{1,k},indices{2,k}) = var{i,k};
    end
    Ap(xInd,yInd)   = AOrigCopy + [y; -y];
    Ap(xInd2:end,:) = AOrigEq - bsxfun(@times,poolsIters(i,:)',  pAOrigEq);
    mxlpsolve('set_mat', lp, Ap);
    mxlpsolve('solve',   lp);
    [res(i,1), res(i,2:end), ~] = mxlpsolve('get_solution', lp);
end
mxlpsolve('delete_lp', lp);
%}
    
if all2SamePools
    % plot 3D sum of 1 output subproblems
    %plot3D('Mesh summed over all outputs',resultsSum,pStep,poolsConnected,'Total obj function');
    % plot 3D total problem for 2 pools
    plot3D('Mesh for original problem',[origPoolsIters, res(:,1)],pStep,poolsConnected2,'Original Obj function');
    %plotDataTable(excelExport,[origPoolsIters res],0,'Numeric values for original problem','fval',...
    %    poolsConnected2,inputPoolConns,directOutputConns,poolOutputConns);  
elseif all1Pool
    figure2d = figure('name','Results','units','normalized','outerposition',[0 0.15 1 0.85]);
    % Create axes
    axes1 = axes('Parent',figure2d,'YGrid','on','XGrid','on',...
        'Position',[0.02 0.079 0.80 0.84]);
    box(axes1,'on');
    hold(axes1,'on');
    % plot subproblems
    plot1 = plot(origPoolsIters,allResults,'LineWidth',1,'Parent',axes1);
    % plot total (original problem)
    plot2 = plot(origPoolsIters,res(:,1),'LineWidth',2,'Parent',axes1,...
        'Color', [0.46 0.67 0.188]);
    xlabel('p','FontSize',11);
    ylabel('Opt value','FontSize',11);
    %leg = legend(axes1,'show');
    %set(leg, 'string',[legSubs 'Total Obj'],'FontSize', 10);
    title(sprintf('Orig Objective (in green bold) vs. 1 output subprobs for %d inputs, %d outputs, %d direct inputs',...
        nbInputsTotal,nbOutputs,nbDirectsTotal),'FontSize',11);
    
    % write down configuration of pooling network for the plot
    configText = uicontrol('style','text');
    set(configText,'units','normalized','Position', [0.85 0 0.25 0.84]);
    set(configText,'HorizontalAlignment', 'left', 'FontSize', 10);
    set(configText,'String',sprintf(['Pooling network config:\n\n',...
        'inputVectConc     = \n',sprintf('%5.1f  ', inputVectConcPerPool2{1}),...
        '\ninputVectCost     = \n',sprintf('%5.1f  ', inputVectCostPerPool2{1}),...
        '\nP_U = \n',sprintf('%5.1f  ', P_U),...
        '\noutputProfitVect  = \n',sprintf('%5.1f  ', outputProfitVect),...
        '\nD_U    = \n',sprintf('%5.1f  ', D_U),...
        '\nCj   = \n',sprintf('%5.1f  ', Cj),...
        '\ngj   = \n',sprintf('%5.1f  ', gj),...
        '\ndirectOutputConns = ',...
        strjoin(arrayfun(@(o) strcat('\n[',sprintf('%5.1f  ', directOutputConns{o}),']'),...
        (1:nbOutputs)','UniformOutput',false))]...
        ));
    set(configText,'FontName', 'FixedWidth');
    
    plotDataTable(excelExport, [origPoolsIters res],0,'Numeric values for original problem','fval',...
        poolsConnected2,inputPoolConns,directOutputConns,poolOutputConns);
else
%     allResults = []; xInd = 0; yInd = 0; zInd = 0;
%     fIndices = [0,0,0; cumsum(fNbVars,1)];
%     %fIndices = cumsum(fIndices,2);
%     fInd = cumsum([nbInputsTotal,sum(nbPoolsTotal),nbDirectsTotal]);
%     %nbInputsTotal + sum(nbPoolsTotal) + nbDirectsTotal;
%     for o=1:nbOutputs
%         fOrigAtO = [        fIndices(o,1)+1 :         fIndices(o+1,1),...
%                     fInd(1)+fIndices(o,2)+1 : fInd(1)+fIndices(o+1,2),...
%                     fInd(2)+fIndices(o,3)+1 : fInd(2)+fIndices(o+1,3) ];
%         allResults(:,o) = -res(:,1+fOrigAtO)*fOrig(fOrigAtO)';
% 
%     end
    plotDataTable(excelExport,[origPoolsIters, res],0,'Numeric values for original problem','fval',...
        poolsConnected2,inputPoolConns,directOutputConns,poolOutputConns);
end
toc;
end

function plot3D(figName,data,pStep,poolsConnected,labelZ)
figure('name',figName,'units','normalized','outerposition',[0 0.15 1 0.85]);
x= data(:,1); y= data(:,2); z= data(:,3);
x_edge=floor(min(x)):pStep:ceil(max(x));
y_edge=floor(min(y)):pStep:ceil(max(y));
[X,Y]=meshgrid(x_edge,y_edge);
Z=griddata(x,y,z,X,Y);
surf(X,Y,Z);
xlabel(sprintf('Pool p%d', poolsConnected(1)),'FontSize',11);
ylabel(sprintf('Pool p%d', poolsConnected(2)),'FontSize',11);
zlabel(labelZ,'FontSize',11);
rotate3d on;
end

function plotDataTable(excelExport, res,o,tableTitle,fvals,poolsConnected,inputPoolConns,directOutputConns,poolOutputConns)
nbOutputs = length(poolOutputConns);
ps = arrayfun(@(pool) sprintf('p%d',pool) ,poolsConnected,'UniformOutput',false);
xs = arrayfun(@(pool) arrayfun(@(input) sprintf('x%d%d',input,pool),...
    inputPoolConns{pool},'UniformOutput', false),...
    poolsConnected','UniformOutput', false);
xs = {cat(2, xs{:})};
if o==0 %orig problem
    zs = arrayfun(@(o) arrayfun(@(dirInput) sprintf('z%d%d',dirInput,o),...
        directOutputConns{o},'UniformOutput', false),...
        (1:nbOutputs),'UniformOutput',false);
    zs = {cat(2, zs{:})};
else
    zs = arrayfun(@(dirInput) sprintf('z%d%d',dirInput,o),directOutputConns{o},...
        'UniformOutput', false);
end
if o==0 %orig problem
    ys = arrayfun(@(o) arrayfun(@(pool) sprintf('y%d%d',pool,o),...
        poolOutputConns{o}, 'UniformOutput', false),...
        (1:nbOutputs), 'UniformOutput', false);
    ys = {cat(2, ys{:})};
else
    ys = arrayfun(@(pool) sprintf('y%d%d',pool,o),poolOutputConns{o},...
            'UniformOutput', false);
end
if o==0 %orig problem
    %arrayfun(@(o) sprintf('fval(%d)',o),(1:nbOutputs),'UniformOutput',false),
    header = [ps, 'fval', xs{:}, ys{:}, zs{:}];
else
    header = [ps, sprintf('fval(%d)',o), xs{:}, ys, zs];
end
psVals = arrayfun(@(row) strcat('(',strjoin(arrayfun(@(p) sprintf('%0.2f',res(row,p)), 1:length(poolsConnected),...
                    'UniformOutput',false),','),')'), 1:size(res,1),...
                    'UniformOutput',false);
if excelExport
    fileName = 'Res1PoolAdditiv1.xlsx'; sheet = sprintf('Sheet%d',o+1);
    xlswrite(fileName,res,   sheet,'B2');     %Write data
    xlswrite(fileName,header,sheet,'B1');     %Write column header
    xlswrite(fileName,psVals',sheet,'A2');    %Write row header
else
    figureForOutput = figure('name',tableTitle,'units','normalized','outerposition',[0 0.15 1 0.85]);
    uitable(figureForOutput,'Units','normalized','Position',[0 0 1 1],...
    'FontName','monospace','FontSize',8,'Data',res,'ColumnName',header,...
    'ColumnFormat',repmat({'bank'},1,size(res,2)),...
    'RowStriping','off','ColumnWidth',repmat({50},1,size(res,2)),...
    'RowName',psVals);
end
end
