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


%Lets make the time-step operator
val = [1,-1,1,-1];
%The Diffusion Coefficients
D(1) = GelSimParams.Dh;
D(2) = GelSimParams.Db;
D(3) = GelSimParams.Di;
D(4) = GelSimParams.Da;


[L,bndterms] = ConstrainedBackEulOperatorConstruct(D,dt,val);

%First, we will get ready to update the Hydrogen concentration


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
RHSH(1) = GelSimParams.HydFluxL + GelSimParams.SolValL*GelSimParams.HydExchangeRate*bndterms(3);
RHSH(end) = GelSimParams.HydValR;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Now we'll get ready to update the Bicarbonate Concentration

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
RHSB(1) = GelSimParams.BicFluxL + GelSimParams.SolValL*GelSimParams.BicExchangeRate*bndterms(4);
RHSB(end) = GelSimParams.BicValR;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Now we'll update the (negatively charged) ionic Concentration

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
RHSI(1) = GelSimParams.IonFluxL + GelSimParams.SolValL*GelSimParams.HydExchangeRate*GelSimParams.HydExchangerParam*bndterms(1);
RHSI(end) = GelSimParams.IonValR;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Now we'll update the (positively charged) anion Concentration


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
RHSA(1) = GelSimParams.AniFluxL + GelSimParams.SolValL*GelSimParams.BicExchangeRate*GelSimParams.BicExchangerParam*bndterms(2);
RHSA(end) = GelSimParams.AniValR;


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