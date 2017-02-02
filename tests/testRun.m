%{
model = pdeSetup(1);
initGeometry(model);
setBoundary(model, [-3.5, 3.5], [1 2]);
setInterface(model, 0);
setRegion(model,{'B1','I1'});
setRegion(model,{'I1','B2'});
%plotModel(model);
%
setCoefficients(model,0,1,0);
setForcing(model,@(x,t) 0);
setBC(model,@(t) 2*exp(-3200*(-3.5 - t + 1/2).^2),'B1',-1,'Dirichlet');
setIC(model,@(x) 2*exp(-3200*(x + 1/2).^2));
setIF(model,'I1',[1/2,-1/2]);
setSolution(model,@(x,t)2*exp(-3200*(x - t + 1/2).^2));
%
initDiscretisation(model,1);
%setOperator(model,'DRP(4,2,2,8,pi/2)','DRPoperators.mat','../SBP_operators')
setOperator(model,'SBP(8,4)','operators.mat','../SBP_operators')
for k=1:5
    N(k) = 210*2^k;
    disp(['Run ',num2str(k),' of 5. N = ',num2str(2*N(k))])
    setResolution(model,N(k));
    getStep(model);
    %
    pdePrepare(model);
    E(k) = pdeSolve(model,0,1);
end
N = 2*N;
c = polyfit(log(N),log(E),1);
loglog(N,E,'bs')
hold on
loglog(N,exp(c(2))*N.^c(1),'k--')
disp(['Convergence rate: ',num2str(-c(1))])
%}

% {
%% Geometry
model = pdeSetup(1);
initGeometry(model);
setBoundary(model, [-3.5 3.5], [1 2]);
setInterface(model, -1);
setInterface(model, 1);
setRegion(model,{'B1','I1'});
setRegion(model,{'I1','I2'});
setRegion(model,{'I2','B2'});
plotModel(model);
%% Discretisation
initDiscretisation(model,1);
N = 1000;
setRegRes(model,'omega1',N);
setRegRes(model,'omega2',N);
setRegRes(model,'omega3',N);
getStep(model);
setRegOp(model,'omega1','DRP(6,3,1,8,pi/5)','DRPoperators.mat','../SBP_operators');
setRegOp(model,'omega2','DRP(6,3,1,8,pi/5)','DRPoperators.mat','../SBP_operators');
setRegOp(model,'omega3','DRP(6,3,1,8,pi/5)','DRPoperators.mat','../SBP_operators');
%setOperator(model,'SBP(8,4)','operators.mat','../SBP_operators');
%% Equation
setCoefficients(model,0,1,0);
setForcing(model,@(x,t) 0);
setIC(model,@(x) 2*exp(-3200*(x + 3).^2));
setSolution(model,@(x,t)2*exp(-3200*(x - t + 3).^2));
setBC(model,@(t) 2*exp(-3200*(-3.5 - t + 3).^2),'B1',-1,'Dirichlet');
setIF(model,'I1',[1/2,-1/2]);
setIF(model,'I2',[1/2,-1/2]);
pdePrepare(model);
%% Solution
tic
%pdeSolve(model,0,6,'plot',[-0.5 2.5])
res = pdeSolve(model,0,6)
toc

%for k=1:5
%    N(k) = 125*2^(k-1);
%    disp(['Run ',num2str(k),' of 5. N = ',num2str(3*N(k))])
%    setRegRes(model,'omega1',N(k));
%    setRegRes(model,'omega2',N(k));
%    setRegRes(model,'omega3',N(k));
%    getStep(model);
%    %
%    pdePrepare(model);
%    tic
%    res = pdeSolve(model,0,6);
%    toc
%    ES84(k) = res.error;
%end
%}

%%
% {
v = [0 0; 0 0.5; 0.5 0.5; 0.5 1; 1 1; 1 0];
f = 1:6;
model = pdeSetup(2);
initGeometry(model);
setBoundary(model, v, f);
%plotGeometry(model.geometry)
setInterface(model,[0.5 0.5; 1 0.5; 0.5 0],[1 2; 1 3]);
%plotGeometry(model.geometry)
setRegion(model,{'B1','B2','I1','I3'});
setRegion(model,{'I3','I1','I2','B6'});
setRegion(model,{'I1','B4','B5','I2'});
plotModel(model)
%}