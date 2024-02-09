%This script begins a simulation. Generally this is passed as an argument
%to a background matlab process via the command line. 
clear all

%These three variables (first 2 are structs) must be declared global for
%the simulation. Comments at the bottom of this file outline the format
%required from each
global GelState GelSimParams rescaled

%Load the three above variables from InitialData.mat (must be in pwd)
load('InitialData.mat')

%Add the directory where all source files reside to the path
addpath('~/Matlab/ExplicitTimeBnd/')

%Call the main time integration loop
MainTimeLoop


%GelSimParam is a struct containing the following variables:
%
% Ncell: Number of cell center points in domain (tied to spatial resolution)
% dt: size of simulation timestep
% SimName: name of simulation. Dictates the filename of output data
% L: Length of spatial domain
% Nedges: Number of cell edges. Must be Ncell-1
% hx: size of spatial step. Must be L/Ncell
% Tmax: Maximum time to simulate to.
% writesteps: how many timesteps between each save to disk
% ThetaTol: Minimum allowed volume fraction
% Dh: Diffusion coefficient of hydrogen species
% Db: Diffusion coefficient of bicarbonate species
% Di: Diffusion coefficient of cation species-05
% Da: Diffusion coefficient of anion species
% De: Diffusion coefficient of enzyme (pepsin) species (unused at this time)
% Dp: Diffusion coefficient of precursor (pepsinogen) species (unused at this time)
% muSol: Viscosity of solvent phase (unused at this time)
% muNet: Viscosity of network phase (unused at this time)
% KbTPerV: Characteristic energy scale for chemical potentials (unused at this time)
% xi: Drag between solvent and network phase (unused at this time)
% SolValL: Volume fraction at which solvent enters the domain from wall
% NetValL: Volume fraction at which network enters the domain from wall (must be 1-SolValL) 
% SolVelValL: Velocity at which solvent enters the domain from wall
% NetVelValL: Velocity at which network enters the domain from wall
% SolVelFluxR: Boundary condition (neumann) at right for solvent velocity solve (unsed at this time)
% NetVelFluxR: Boundary condition (neumann) at right for network velocity
% solve (unsed at this time)
% HydValR: Concentration of hyd at lumen (for right boundary condition)
% BicValR: Concentration of bic at lumen (for right boundary condition)
% IonValR: Concentration of cation at lumen (for right boundary condition)
% AniValR: Concentration of anion at lumen (for right boundary condition)
% ValH: Valence of hydrogen species
% ValB: Valence of bicarbonate species
% ValA: Valence of anion species
% ValI: Valence of cation species
% Kbind: Buffering reaction rate
% StopTol: Tolerance for determining ``steady state''
% BicExchangeRate: Rate constant of bic/anion exchanger at the wall
% BicExchangerParam: Offset constant of bic/anion exchanger at the wall
% HydExchangeRate: Rate constanat of hyd/cation exchanger
% HydExchangerParam: Offset constatne of hyd/cation exchanger


% GelState is a struct containing the following variables
% 
% Xcell: [Ncellx1 array] location of cell center points (interior only)
% Xedge: [Nedgex1 array] location of cell edge points (interior only)
% XcellExtend: [Ncell+2x1 array] cell center points with ghost points
% XedgeExtend: [Nedge+2x1 array] cell edge points with ghost points
% Time: Current time of the simulation
% Hconc: [Ncell+2x1 array] Hyd concentration at cell centers @ current time
% Bconc: [Ncell+2x1 array] Bic concentration at cell centers @ current time
% Iconc: [Ncell+2x1 array] Cation concentration at cell centers @ current time
% Aconc: [Ncell+2x1 array] Anion concentration at cell centers @ current time
% Hold: [Ncell+2x1 array] Hyd concentration at cell centers @ previous time
% Bold: [Ncell+2x1 array] Bic concentration at cell centers @ previous time
% Iold: [Ncell+2x1 array] Cation concentration at cell centers @ previous time
% Aold: [Ncell+2x1 array] Anion concentration at cell centers @ previous time
% DPsi: [Nedge+2x1 array] Electric potential gradient at cell edges
% ThetaS: [Ncell+2x1 array] Solvent volume fraction at cell centers
% ThetaN: [Ncell+2x1 array] Network volume fraction at cell centers (must be 1-ThetaS)
% ThetaSCorr: [Ncellx1 array] Corrected Solvent Volume fraction (to avoid dividing by zero)
% ThetaNCorr: [Ncellx1 array] Corrected Network Volume fraction (to avoid dividing by zero)
% USol: [Nedge+2x1 array] Velocity of solvent at cell edges (current time)
% USolOld: [Nedge+2x1 array] Velocity of solvent at cell edges (previous time)
% UNet: [Nedge+2x1 array] Velocity of network at cell edges (current time)
% UNetOld: [Nedge+2x1 array] Velocity of network at cell edges (previous time)
% Pres: [Ncellx1 array] Hydrodynamic pressure at interior cell edges
% HRHScur: [Ncellx1 array] Holding array for right hand side of hyd evolution equation at current time
% HRHSold: [Ncellx1 array] Holding array for right hand side of hyd evolution equation at previous time
% BRHScur: [Ncellx1 array] Holding array for right hand side of bic evolution equation at current time
% BRHSold: [Ncellx1 array] Holding array for right hand side of bic evolution equation at previous time
% IRHScur: [Ncellx1 array] Holding array for right hand side of cation evolution equation at current time
% IRHSold: [Ncellx1 array] Holding array for right hand side of cation evolution equation at previous time
% ARHScur: [Ncellx1 array] Holding array for right hand side of anion evolution equation at current time
% ARHSold: [Ncellx1 array] Holding array for right hand side of anion evolution equation at previous time

%The variable rescaled is an [Ncellx1 array] which contains the source term
%for hydrogen and anion concentration variables. The simulation will
%automatically include it in HRHScur and ARHScur as necessary. 