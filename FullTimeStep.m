%This function will take a single time step of the simulation.  
%
% function syntax:
%
%     FullTimeStep
%
%
%     inputs:
%         none 
%     output:
%         None 


function FullTimeStep
%Variable used to see if this is the first time step. Will only be used in
%future code when we move to 2nd order time integration
persistent FirstStep
%Get those global structs
global GelState GelSimParams

%Finally, we need to populate the reaction and potential terms for the
%Reac-Adv-Diff equations for volumeless species.
PopulateDiffRHS


%And then we'll go ahead and time step those

%Much of this is code written with an eye to the future. We will be using
%2nd order SBD time integration for implicit terms in the future. Right
%now, the commented code ensures that we use Backward Euler at all time
%steps.
switch isempty(FirstStep)
    case 1  %If this is the very first step of the simulation, we cannot use SBD2
            %We will have to evolve the R-D-A equations using a single step
            %method.
        
%         disp('Taking the first time step using SBD1')
        %We will use a 'first order' IMEX method ()    
        ConstrainedBackEulStep
        
        %And now we update 'FirstStep' to make sure that for the rest of
        %the simulation, we'll use SBD2
%         FirstStep = 0;
        
%         disp('Ok, thats done. We shouldnt need to again')

        
    case 0 %If this is not the very first step, we will use SBD2

        %Update the reaction-diffusion-advection equation for the volumeless
        %species
        ConstrainedSBD2Step
        
end

%%%WE'RE DONE!
end
