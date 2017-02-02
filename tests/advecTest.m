%%
advec = pdeSetup(1);
initGeometry(advec);
setBoundary(advec, [0, 1], [1 2]);
%plotModel(advec);
initDiscretisation(advec);
%setResolution(advec,20);
%getStep(advec);
setOperator(advec,'DRP(2,1,3,8,pi/2)','DRPoperators.mat','../SBP_operators')
%setOperator(advec,'SBP(8,4)','operators.mat','../SBP_operators')
%
setCoefficients(advec,0,1,0);
setForcing(advec,@(x,t) 0);
setBC(advec,@(t) -sin(2*pi*t),'B1',-1,'Dirichlet');
setIC(advec,@(x) sin(2*pi*x));
setSolution(advec,@(x,t) sin(2*pi*(x-t)));
%
for k=1:5
    N(k) = 10*2^k;
    setResolution(advec,N(k));
    getStep(advec);
    pdePrepare(advec);
    E(k)=pdeSolve(advec,0,2.3);
end
%pdeSolve(advec,0,1,'plot',[-1.5 1.5])
c = polyfit(log(N),log(E),1);
loglog(N,E,'bs')
hold on
loglog(N,exp(c(2))*N.^c(1),'k--')
disp(['Convergence rate: ',num2str(-c(1))])