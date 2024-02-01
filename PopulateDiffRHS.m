%This function populates the right hand side of the reaction-diffusion-advection 
%systems for the diffusive volumeless species.
%
% function syntax:
%
%     PopulateDiffRHS
%
%
%     inputs:
%         none 
%     output:
%         none 


function PopulateDiffRHS
    global GelState rescaled


    %Lets move the RHS terms from current to old
    %This is for future versions of the code where we use 2nd order SBD
    %time integration
    GelState.HRHSold = GelState.HRHScur;
    GelState.BRHSold = GelState.BRHScur;
    GelState.IRHSold = GelState.IRHScur;
    GelState.ARHSold = GelState.ARHScur;
    GelState.PRHSold = GelState.PRHScur;
    GelState.ERHSold = GelState.ERHScur;


    %The source term explicitly appears in the RHS of the equation. 
    %Buffering reactions no longer appear here because we are now handling
    %them semi-implicitly (see construction of 'HydAdj' and 'BicAdj' in
    %ConstrainedBackEulOperatorConstruct)
    GelState.HRHScur = rescaled;
    GelState.ARHScur = rescaled;

end
