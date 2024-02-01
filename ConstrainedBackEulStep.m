%This function evolves a fully (electrically coupled) 
%reaction-diffusion-advection system one step using my new constrained
%Backward Euler time integrator. 
%
% function syntax:
%
%     ConstrainedBackEulStep
%
%
%     inputs:
%         none 
%     output:
%         none 


function ConstrainedBackEulStep

% Lets 'import' the two big global structs
global GelState GelSimParams

%We'll need the time step
dt = GelSimParams.dt;

%The current solvent velocity field
SolVeloc = GelState.USol;

%First, we will get ready to update the Hydrogen concentration

%The Diffusion Coefficient
D(1) = GelSimParams.Dh;
%And finally the current concentration
conccur = GelState.Hconc;

%Explicitly evaluate the advection term
advcur = AdvectionEvaluate(conccur,SolVeloc);

%Add the reaction term
explcur = GelState.HRHScur - advcur;

%Begin construction of RHS for implicit solve
RHSH = 0*conccur;

%Populate the entries which correspond to cells within the computational
%domain
RHSH(2:end-1) = conccur(2:end-1) + dt*explcur;
    
%Now we need to populate the entries which correspond to ghost cells with
%the appropriate fluxes for Boundary Conditions.
RHSH(1) = GelSimParams.HydFluxL;
RHSH(end) = GelSimParams.HydValR;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Now we'll get ready to update the Bicarbonate Concentration

%The Diffusion Coefficient
D(2) = GelSimParams.Db;
%And finally the current concentration
conccur = GelState.Bconc;

%Explicitly evaluate the advection term
advcur = AdvectionEvaluate(conccur,SolVeloc);

%Add the reaction term
explcur = GelState.BRHScur - advcur;

%Begin construction of RHS for implicit solve
RHSB = 0*conccur;

%Populate the entries which correspond to cells within the computational
%domain
RHSB(2:end-1) = conccur(2:end-1) + dt*explcur;
    
%Now we need to populate the entries which correspond to ghost cells with
%the appropriate fluxes for Boundary Conditions.
RHSB(1) = GelSimParams.BicFluxL;
RHSB(end) = GelSimParams.BicValR;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Now we'll update the (negatively charged) ionic Concentration

%The Diffusion Coefficient
D(3) = GelSimParams.Di;
%And finally the current concentration
conccur = GelState.Iconc;

%Explicitly evaluate the advection term
advcur = AdvectionEvaluate(conccur,SolVeloc);

%Add the reaction term
explcur = GelState.IRHScur - advcur;

%Begin construction of RHS for implicit solve
RHSI = 0*conccur;

%Populate the entries which correspond to cells within the computational
%domain
RHSI(2:end-1) = conccur(2:end-1) + dt*explcur;

%Now we need to populate the entries which correspond to ghost cells with
%the appropriate fluxes for Boundary Conditions.
RHSI(1) = GelSimParams.IonFluxL;
RHSI(end) = GelSimParams.IonValR;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Now we'll update the (positively charged) anion Concentration

%The Diffusion Coefficient
D(4) = GelSimParams.Da;
%And finally the current and 'old' concentrations
conccur = GelState.Aconc;

%Explicitly evaluate the advection terms
advcur = AdvectionEvaluate(conccur,SolVeloc);

%Add the reaction terms
explcur = GelState.ARHScur - advcur;

%Begin construction of RHS for implicit solve
RHSA = 0*conccur;

%Populate the entries which correspond to cells within the computational
%domain
RHSA(2:end-1) = conccur(2:end-1) + dt*explcur;

%Now we need to populate the entries which correspond to ghost cells with
%the appropriate fluxes for Boundary Conditions.
RHSA(1) = GelSimParams.AniFluxL;
RHSA(end) = GelSimParams.AniValR;


val = [1,-1,1,-1];

L = ConstrainedBackEulOperatorConstruct(D,dt,val);

RHS = [RHSH;RHSB;RHSI;RHSA;zeros(GelSimParams.Ncell+1,1)];

newconcs = L\RHS;

%Lets pluck out the entries corresponding to hydrogen
concnew = newconcs(1:GelSimParams.Ncell+2);
%Current becomes old, and new becomes current
GelState.Hold = GelState.Hconc;
GelState.Hconc = concnew;

%Now we'll pluck out the entries corresponding to bicarbonate
concnew = newconcs(GelSimParams.Ncell+3:2*GelSimParams.Ncell+4);
%Current becomes old, and new becomes current
GelState.Bold = GelState.Bconc;
GelState.Bconc = concnew;

%Now we'll pluck out the entries corresponding to negative Ions
concnew = newconcs(2*GelSimParams.Ncell+5:3*GelSimParams.Ncell+6);
%Current becomes old, and new becomes current
GelState.Iold = GelState.Iconc;
GelState.Iconc = concnew;

%Now we'll pluck out the entries corresponding to positive Anions
concnew = newconcs(3*GelSimParams.Ncell+7:4*GelSimParams.Ncell+8);
%Current becomes old, and new becomes current
GelState.Aold = GelState.Aconc;
GelState.Aconc = concnew;

%And finally, lets pluck out the electric potential gradient
GelState.DPsi = newconcs(4*GelSimParams.Ncell+9:end);


 
end