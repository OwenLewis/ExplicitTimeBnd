%This function construts the ``gell averaged'' operator for the final constrained 
%step of a backward Euler iteration of 
%diffusive species in the two-phase gel model. These operators are then used
%to advance the reaction-diffusion-advection system. For the purposes of
%this construction, we assume that the R-D-A system is subject to Robin
%boundary condtions at the left and homogeneous Neumann conditions at the
%right. We also assume that the solvent volume fraction is subject to a
%given advective flux at the left and outflow boundary conditions at the
%right.
%
% function syntax:
%
%     [Lfut,bndterms] = ConstrainedBackEulOperatorConstruct(DiffCoeff,dt,val)
%
%
%     inputs:
%         DiffCoeff is an array which defines the diffusion coefficients for
%           the species which these operators will update the constrained
%           heat equation
%         dt is the time update step size. 
%     output:
%         Lfut is the "(Ncell+2)*size(DiffCoeff)+Nedge" square sparce matrix 
%           which applies  the operation of the SBD scheme to the future time
%           step. This one will be inverted to solve for the new values.
%           (The first and last entries of the Ncell + 2 blocks are related 
%           to ghost cells and bounary condtions). 
%         bndterms is the entry vector which is an approximation of the ion
%           exchanger model terms "evaluate" at the future time step. This 
%           will need to be added to the right hand side vector to take a
%           time step.


function [Lfut,bndterms] = ConstrainedBackEulOperatorConstruct(DiffCoeff,dt,val)

%Lets 'import' the two big global structs
global GelState GelSimParams
% 'time to form gel velocity operator'

%Here are some parameters we need to define the sizes of things
hx = GelSimParams.hx;
Ncell = GelSimParams.Ncell;
Nedge = GelSimParams.Nedges;

%And the diffusion coefficient, whos name i don't want to type over & over
D = DiffCoeff;


%First, lets get the SBD2 Time Step operators for each species
Lfut1 = BackEulOperatorConstruct(D(1),dt,-GelSimParams.HydExchangeRate*GelSimParams.HydExchangerParam);
%  -GelSimParams.HydExchangeRate*GelSimParams.HydExchangerParam
Lfut2 = BackEulOperatorConstruct(D(2),dt,-GelSimParams.BicExchangeRate*GelSimParams.BicExchangerParam);
% -GelSimParams.BicExchangeRate*GelSimParams.BicExchangerParam
Lfut3 = BackEulOperatorConstruct(D(3),dt,-GelSimParams.HydExchangeRate);
% -GelSimParams.HydExchangeRate
Lfut4 = BackEulOperatorConstruct(D(4),dt,-GelSimParams.BicExchangeRate);
% -GelSimParams.BicExchangeRate

%Now we'll construct the first diff. operators which act on the electric 
%potential

%First, the action of electric potential on hydrogen
%Operator requires an approx of future concentration
FutConcExtrap = 2*GelState.Hconc - GelState.Hold; %Heres an extrap to next time step
HydBndFut = mean(GelState.Hconc(1:2)); %This will be necessary for ion exchanger evaluation

%Lets set asside the diagonal of the new operator that takes care of the
%buffering reaction implicitly. 
BicAdjDiag = GelSimParams.dt*GelSimParams.Kbind*FutConcExtrap;
BicAdjDiag(1) = 0;
BicAdjDiag(end) = 0;
BicAdj = spdiags(BicAdjDiag,0,Ncell+2,Ncell+2);

%Now we interpolate the concentration times valence times volume fraction
%to cell edges
weightedges = D(1)*val(1)*interp1(GelState.XcellExtend,FutConcExtrap.*GelState.ThetaS,GelState.XedgeExtend,'linear');

%Make the diagonals from these values
DImpLDiag = dt*weightedges(1:end-1)./(GelState.ThetaS(2:end-1)*hx);
DImpUDiag = -dt*weightedges(2:end)./(GelState.ThetaS(2:end-1)*hx);

DImpUDiag = [-D(1)*val(1)*mean(FutConcExtrap(1:2));DImpUDiag];
DImpLDiag = [DImpLDiag;0];

%And assemble the operator
Dfut1 = spdiags([DImpLDiag';DImpUDiag']',-1:0,Ncell+2,Nedge+2);
Dfut1(1,:) = Dfut1(1,:)*GelSimParams.SolValL;
Dfut1(end,:) = Dfut1(end,:)*mean(GelState.ThetaS(end-1:end));


%Now on species two
%Operator requires an approx of future concentration
FutConcExtrap = 2*GelState.Bconc - GelState.Bold; %Heres an extrap to next time step
BicBndFut = mean(GelState.Bconc(1:2)); %This will be necessary for ion exchanger evaluation

%Lets set asside the diagonal of the new operator that takes care of the
%buffering reaction implicitly. 
HydAdjDiag = GelSimParams.dt*GelSimParams.Kbind*FutConcExtrap;
HydAdjDiag(1) = 0;
HydAdjDiag(end) = 0;
HydAdj = spdiags(HydAdjDiag,0,Ncell+2,Ncell+2);

%Now we interpolate the concentration times valence times volume fraction
%to cell edges
weightedges = D(2)*val(2)*interp1(GelState.XcellExtend,FutConcExtrap.*GelState.ThetaS,GelState.XedgeExtend,'linear');

%Make the diagonals from these values
DImpLDiag = dt*weightedges(1:end-1)./(GelState.ThetaS(2:end-1)*hx);
DImpUDiag = -dt*weightedges(2:end)./(GelState.ThetaS(2:end-1)*hx);

DImpUDiag = [-D(2)*val(2)*mean(FutConcExtrap(1:2));DImpUDiag];
DImpLDiag = [DImpLDiag;0];

%And assemble the operator
Dfut2 = spdiags([DImpLDiag';DImpUDiag']',-1:0,Ncell+2,Nedge+2);
Dfut2(1,:) = Dfut2(1,:)*GelSimParams.SolValL;
Dfut2(end,:) = Dfut2(end,:)*mean(GelState.ThetaS(end-1:end));

%Now on species three
%Operator requires an approx of future concentration
FutConcExtrap = 2*GelState.Iconc - GelState.Iold; %Heres an extrap to next time step
IonBndFut = mean(GelState.Iconc(1:2)); %This will be necessary for ion exchanger evaluation

%Now we interpolate the concentration times valence times volume fraction
%to cell edges
weightedges = D(3)*val(3)*interp1(GelState.XcellExtend,FutConcExtrap.*GelState.ThetaS,GelState.XedgeExtend,'linear');

%Make the diagonals from these values
DImpLDiag = dt*weightedges(1:end-1)./(GelState.ThetaS(2:end-1)*hx);
DImpUDiag = -dt*weightedges(2:end)./(GelState.ThetaS(2:end-1)*hx);

DImpUDiag = [-D(3)*val(3)*mean(FutConcExtrap(1:2));DImpUDiag];
DImpLDiag = [DImpLDiag;0];

%And assemble the operator
Dfut3 = spdiags([DImpLDiag';DImpUDiag']',-1:0,Ncell+2,Nedge+2);
Dfut3(1,:) = Dfut3(1,:)*GelSimParams.SolValL;
Dfut3(end,:) = Dfut3(end,:)*mean(GelState.ThetaS(end-1:end));

%Now species four
%Operator requires an approx of future concentration
FutConcExtrap = 2*GelState.Aconc - GelState.Aold; %Heres an extrap to next time step
AniBndFut = mean(GelState.Aconc(1:2)); %This will be necessary for ion exchanger evaluation

%Now we interpolate the concentration times valence times volume fraction
%to cell edges
weightedges = D(4)*val(4)*interp1(GelState.XcellExtend,FutConcExtrap.*GelState.ThetaS,GelState.XedgeExtend,'linear');

%Make the diagonals from these values
DImpLDiag = dt*weightedges(1:end-1)./(GelState.ThetaS(2:end-1)*hx);
DImpUDiag = -dt*weightedges(2:end)./(GelState.ThetaS(2:end-1)*hx);

DImpUDiag = [-D(4)*val(4)*mean(FutConcExtrap(1:2));DImpUDiag];
DImpLDiag = [DImpLDiag;0];

%And assemble the operator
Dfut4 = spdiags([DImpLDiag';DImpUDiag']',-1:0,Ncell+2,Nedge+2);
Dfut4(1,:) = Dfut4(1,:)*GelSimParams.SolValL;
Dfut4(end,:) = Dfut4(end,:)*mean(GelState.ThetaS(end-1:end));

%Finally, we will need to construct the scaled identity matrices that go on
%the diagonal.
EyeDiag = ones(Nedge+2,1);

ModEye = spdiags(EyeDiag,0,Nedge+2,Ncell+2);

Eye1 = val(1)*ModEye;
Eye2 = val(2)*ModEye;
Eye3 = val(3)*ModEye;
Eye4 = val(4)*ModEye;



%Lets put the whole damn thing together
A = sparse(Ncell+2,Ncell+2);
B = sparse(Nedge+2,Nedge+2);
HtoI = A;
ItoH = A;
BtoA = A;
AtoB = A;

%These are the terms that couple counter-ions through boundary equations
% HtoI(1,1:2) = -GelSimParams.SolValL*GelSimParams.HydExchangeRate/2; 
% ItoH(1,1:2) = -GelSimParams.SolValL*GelSimParams.HydExchangeRate*GelSimParams.HydExchangerParam/2;
% BtoA(1,1:2) = -GelSimParams.SolValL*GelSimParams.BicExchangeRate/2;
% AtoB(1,1:2) = -GelSimParams.SolValL*GelSimParams.BicExchangeRate*GelSimParams.BicExchangerParam/2;
% keyboard


Lfut = [Lfut1+HydAdj,A,A,A,Dfut1;A,Lfut2+BicAdj,A,A,Dfut2;A,A,Lfut3,A,Dfut3;A,A,A,Lfut4,Dfut4;Eye1,Eye2,Eye3,Eye4,B];
bndterms = [HydBndFut;BicBndFut;IonBndFut;AniBndFut];


end