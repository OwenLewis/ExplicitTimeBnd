%This is the big one! Integrates the system forward in time
%
% function syntax:
%
%     MainTimeLoop
%
%
%     inputs:
%         none 
%     output:
%         None 


function MainTimeLoop

%Global structs are requires
global GelState GelSimParams


%Localize a few standard parameters
dt = GelSimParams.dt;
Tmax = GelSimParams.Tmax - GelState.Time;

%calculate how many time steps we'll need to take
StepsToTake = ceil(Tmax/dt);
time = GelState.Time;

%Write out the system variables before we get started
writeout

%Keep a record of current ionic concentrations
savedh = GelState.Hconc;
savedb = GelState.Bconc;
saveda = GelState.Aconc;
savedi = GelState.Iconc;

%Variables related to how often we will check to see if steady state has
%been reached
checktimes = 10;
bottomcorrect = 1e-10;
changeplot = 0;

for M = 1:StepsToTake
   
   %Take a time step, and incriment our internal clock
   FullTimeStep
   GelState.Time = GelState.Time + dt;
   
   %Lets check and see if we need to write out before evolving the system
   if mod(M,GelSimParams.writesteps) <= dt/2

      sprintf('We have taken %0.7d steps and will write out',M)
      
       %If we do, write everything to file
       writeout
   end
   
   %If it is time to check for steady state..
   if mod(dt*M,checktimes) <= dt/2
       %Calculate absolute change in concentrations
       changeH = abs(GelState.Hconc - savedh);
       changeB = abs(GelState.Bconc - savedb);
       changeI = abs(GelState.Iconc - savedi);
       changeA = abs(GelState.Aconc - saveda);
       
       %Calculate relative change in concentrations
       relH = 2*changeH./(abs(GelState.Hconc + savedh)+bottomcorrect);
       relB = 2*changeB./(abs(GelState.Bconc + savedb)+bottomcorrect);
       relI = 2*changeI./(abs(GelState.Iconc + savedi)+bottomcorrect);
       relA = 2*changeA./(abs(GelState.Aconc + saveda)+bottomcorrect);

       %Calculate maximum relative change in each
       maxH = max(relH(isfinite(relH)));
       maxB = max(relB(isfinite(relB)));
       maxI = max(relI(isfinite(relI)));
       maxA = max(relA(isfinite(relA)));
       
       %Write that information out for post-processing
       changefile = sprintf('%s.%s.%08.1f.mat',GelSimParams.SimName,'change',GelState.Time);
       save(changefile,'relH','relB','relI','relA')
       
       
       %Old legacy code from development stage related to plotting relative
       %changes of ion concentrations
       if changeplot
           figure(1)
           semilogy(GelState.XcellExtend,relH,'r-',GelState.XcellExtend,relB,'b-',GelState.XcellExtend,relI,'c--',GelState.XcellExtend,relA,'m--','LineWidth',2)
           title(sprintf('Time = %4.4f',GelState.Time),'FontSize',16);
           legend('Hydrogen','Bicarbonate','Cations','Anions','Location','Best')
           ylim([1e-10 1e1])
           drawnow
       end
   
       %If we've met stoping criteria, break the time loop
       if maxH < GelSimParams.StopTol && maxB < GelSimParams.StopTol && maxI < GelSimParams.StopTol && maxA < GelSimParams.StopTol
           disp('We are stopping early because movement is less than tolerance')
           break
       end
       
       savedh = GelState.Hconc;
       savedb = GelState.Bconc;
       saveda = GelState.Aconc;
       savedi = GelState.Iconc;
   end

end


%Finally, lets write out at the very end
writeout

end

%Helper function which handles writing the entire system state to disk
function writeout
global GelSimParams GelState rescaled
persistent writes

if isempty(writes)
    writes = 0;
end

%Lets construct the file name based on how many times we've already written
filename = sprintf('%s.%0.7d.mat',GelSimParams.SimName,writes);
       
%Write out and incriment the number of writes we've completed
save(filename,'GelSimParams','GelState','rescaled')


writes = writes+1;
end
