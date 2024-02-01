%This function evaluates the (explicit) advection terms of a
%reaction-advection-diffusion system. The scheme used is a very simple
%centered differencing scheme. Obviously, this should not be used as part
%of a Forward-Euler time stepping scheme. It will be unconditionally
%unstable. However, as part of an IMEX time integration scheme, it will be
%quite useful. 
%
% function syntax:
%
%     advectionTerms = AdvectionEvaluate(conc,veloc)
%
%
%     inputs:
%         conc is a length Ncell+2 array which contains the concentrations
%           (at cell centeres) of the species which is being advected. The
%           first and last entries (ghost cells) should have already been 
%           populated with whatever values respect the given boundary
%           conditions
%
%         veloc is a length Nedges+2 array which contains the velocites (at
%           cell edges) with which the species is being advected. Again,
%           the first and last entries (ghost cells) should have aleady
%           been populated with whatever values respect the B.C's
%     output:
%         advectionTerms is a length Ncell array defined at cell centers
%           (excluding ghost cells) which contains the explicit evaluation
%           of the advection of the species 'conc' at speed 'veloc'

function advectionTerms = AdvectionEvaluate(conc,veloc)


%Lets 'import' the two big global structs
global GelState GelSimParams

flag = 2;

%Here are some parameters we need to define the sizes of things
hx = GelSimParams.hx;

%Now we need some data at points outside the computational domain
%Now we put together the expanded spacial grid.
expandedXedge = [0;GelState.Xedge;GelSimParams.L];
expandedXcell = [-hx/2;GelState.Xcell;GelSimParams.L+hx/2];


switch flag
    %Case 1 is a simple centered finite difference scheme. It is known to
    %be unstable for purely advective problems, but works fine for
    %advection diffusion-type problems with small Peclet number. We are not
    %currently using this code at all, but it was helpful for testing
    %purposes
    case 1 
        %We will need to interpolate the concentration function to the cell edges
        %where the velocities live. Lets get it out of the way.
        %We have a choice of how to do that
        % concEdge = interp1(expandedXcell,conc,expandedXedge,'pchip');
        concEdge = interp1(expandedXcell,conc,expandedXedge,'linear');


        %Now we multiply by the velocities
        foo = concEdge.*veloc;

        %And take a standard centered finite difference to evaluate the advection
        advectionTerms = diff(foo)/hx;
        
        
        
        %Case 2 is a very simple upwinding scheme to calculate the
        %advection term. As written, it assumes that all velocities are
        %positive (i.e. the upwind direction is to the left). However there
        %is a check to alert the user if this is not true.
    case 2
        
        %We should check to make sure no downwinding is occuring.
        down = sum(veloc < 0);
        if down > 0
            disp('Woops, you are downwinding')
        end
        %Flux through this edge is velocity times the concentration in the
        %cell center TO THE LEFT
        flux = conc(1:end-1).*veloc;
        %Do a centered difference of flux to get the term you want
        advectionTerms = diff(flux)/hx;
        
        
        %Case 3 is currently unused. It was for testing purposes and gives
        %the exactly correct answer for the case of zero Peclet number
    case 3
        advectionTerms = 0*conc(2:end-1);
        
end

end
