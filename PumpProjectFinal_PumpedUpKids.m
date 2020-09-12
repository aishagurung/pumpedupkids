% -----------------------------------------------
% Program Number: Fluids Pump Project
% Submitted By: Aishwarya Gurung
%
% Credit to:
%
% Submitted On: 4/28/2016
%
% By submitting this program with my name,
% I affirm that the creation and modification
% of this program is primarily my own work.
% -------------------------------------------

clear;

clc;

%VARIABLES USED TO CHECK THE PROGRAM:

%P_1 = 0; %psig
%P_2 = 120; %psig
%d_1 = (1.5)/12; %ft
%d_2 = ((1.5)*0.4)/12; %ft
%d_3 = (d_2)*0.5; %ft
%z_1 = 0; %ft
%z_2 = 10; %ft
%Vol_Flow_Rate = 0.1114; %ft^3/s
%L_1 = 1; %ft
%L_2 = 15; %ft
%L_3 = 75; %ft
%L_4 = 2; %ft
%L_5 = 25; %ft
%L_6 = 3.5; %ft


%OPENING AN EXCEL FILE

filename = 'Pump_project.xlsx'; % Changing this allows the user to work with different files and can use extensions such as C:\

inputsheet = 1;

outputsheet = 2;


%HARDCODED VARIABLES

E = 0.000005; %ft roughness of copper
K_L_90_flanged = 0.3;
K_L_90_threaded = 0.9;
K_L_globe = 10;
K_L_contraction = 0.04; %hardcoded for 45 degrees
K_L_nozzle = 0.04; %hardcoded for 45 degrees
g = 32.174; %ft/s^2
viscosity = 0.00059808; %lb/ft*s
density = 62.30; %lb/ft^3
P_Vapour = (0.3632*144); %psfa


%OPENING MESSAGE

waitfor(msgbox('Welcome to Pump Calculator 2016! Before beginning, please ensure that the MATLAB code has been saved to the desktop to ensure proper printing of the results to an Excel file. In a short while, you will be asked for values pertaining to your diagram. Please input values as they are requested.', 'Pump Calculator 2016'));


%ASK FOR DIAMETERS

%ask for number of diameters

n = 0;

while n <= 0;
    
    x = inputdlg('Input the number of varying diameters: ','Pumped',[1,20]);
    
    n = str2num(x{:}); %#ok<*ST2NM>
    
    clc
    
end

clear x

%ask for values of diameters

for k = 1:n;
    
    diameter(k) = 0;
    
    while diameter(k) <= 0;
        
        prompt = ['Enter diameter',num2str(k),' in feet:'];
        
        dia = inputdlg(prompt,'Pumped',[1,20]);
        
        diameter(k) = str2num(dia{:});
        
    end
    
    k = k+1;
    
    clc
    
end

clear dia


%ASK FOR VOLUMETRIC FLOW RATE

Vol_Flow_Rate = 0;

while Vol_Flow_Rate <= 0;
    
    Vol_Flow = inputdlg('Input volumetric flow rate in gallons per minute: ','Pumped',[1,20]);
    
    VFL = (str2num(Vol_Flow{:}));
    
    Vol_Flow_Rate = (VFL/448.8);
    
    clc
    
end

clear Vol_Flow


%CALCULATE AREAS

Area = pi*(((diameter).^2)/4); %ft^2


%CALCULATE VELOCITIES

for m = 1:n;
    
    Velocity(m) = Vol_Flow_Rate/Area(m); %ft/s
    
    m = m+1;
    
end


%FIND REYNOLDS NUMBER

m = 1;

for m = 1:n;
    
    R_N(m) = (density*Velocity(m)*diameter(m))/viscosity;
    
    m = m+1;
    
end


%CHECK FOR LAMINAR/TURBULENT FLOW & FIND FRICTION FACTOR
%assume no transitioning flow

%Check for convergence and find friction factor

m = 1;

for m = 1:n
    
    if R_N(m) < 4000
        
        f_f(m) = 64/R_N(m);
        
        category_laminar(m,1) = {'YES'};
        
        category_turbulent(m,1) = {'NO'};
        
    elseif R_N(m) >= 4000
        
        category_laminar(m,1) = {'NO'};
        
        category_turbulent(m,1) = {'YES'};
        
        f_min = 0;
        
        f_max = 0.1;
        
        leftover = 2;
        
        while abs(leftover) > 0.00001 %Colebrook's equatiom
            
            f_guess = (f_min + f_max)/2;
            
            leftover = -2*log10( ( E/(3.7*diameter(m)) ) + ( 2.51/(R_N(m)*sqrt(f_guess)) ) ) - (1/sqrt(f_guess));
            
            if leftover < 0
                
                f_min = f_guess;
                
            else
                
                f_max = f_guess;
                
            end
            
            f_f(m) = f_guess;
            
        end
        
    end
    
end

%Ask for total length per diameter

for k = 1:n;
    
    length(k) = -1;
    
    while length(k) < 0
        
        prompt = ['Enter total length for diameter',num2str(k),' in feet:'];
        
        Output = inputdlg(prompt,'Pumped',[1,20]);
        
        length(k) = str2num(Output{:});
        
    end
    
    k = k+1;
    
end

clear len


%MAJOR HEAD LOSSES

%Find h_L_major

maj_n = n-1;

h_L_major = 0;

for k = 1:maj_n;
    
    if length(k) > 0
        
        h_L_maj(k) = (f_f(k)*(length(k)/diameter(k))*(Velocity(k).^2)/(2*g));
        
        h_L_major = h_L_major+h_L_maj(k);
        
    end
    
    k = k+1;
    
end

clear h_L_maj


%MINOR HEAD LOSSES

%assume all globe valves in first diameter

n_valve = -10;

while n_valve < 0;
    
    valve = inputdlg('Input the number of globe valves: ','Pumped',[1,20]);
    
    n_valve = str2num(valve{:});
    
    clc
    
end

clear valve

h_L_minor_globe = n_valve*((K_L_globe*(Velocity(1).^2))/(2*g));

%number of contraction fittings
%contactions happen between two diameters and the last contaction is a nozzle
%contractions are only assumed to be 45 degrees so the program can be formatted
%to include more options

n_contrac = n-3;

h_L_minor_contraction = 0;

for k = 1:n_contrac
    
    h_L_minor_contrac(k) = ((K_L_contraction*(Velocity(k).^2))/(2*g));
    
    h_L_minor_contraction = h_L_minor_contraction+h_L_minor_contrac(k);
    
    k = k+1;
    
end

clear h_L_minor_contrac

%number of flanged elbows that are 90 degrees

for k = 1:n;
    
    n_elbow(k) = -1;
    
    while n_elbow(k) < 0;
        
        prompt = ['Enter number of flanged elbows for diameter',num2str(k),':'];
        
        elbow = inputdlg(prompt,'Pumped',[1,20]);
        
        n_elbow(k) = str2num(elbow{:});
        
    end
    
    k = k+1;
    
    clc
    
end

clear elbow

h_L_minor_90elbow = 0;

for k = 1:n
    
    h_L_minor_90_fl(k) = (((K_L_90_flanged*(Velocity(k).^2))/(2*g))*n_elbow(k));
    
    h_L_minor_90elbow = h_L_minor_90elbow+ h_L_minor_90_fl(k);
    
    k = k+1;
    
end

clear h_L_minor_90_fl n_elbow

%number of threaded elbows that are 90 degrees

for k = 1:n;
    
    n_elbow(k) = -1;
    
    while n_elbow(k) < 0;
        
        prompt = ['Enter number of threaded elbows for diameter',num2str(k),':'];
        
        elbow = inputdlg(prompt,'Pumped',[1,20]);
        
        n_elbow(k) = str2num(elbow{:});
        
    end
    
    k = k+1;
    
    clc
    
end

clear elbow

for k = 1:n
    
    h_L_minor_90_th(k) = (((K_L_90_threaded*(Velocity(k).^2))/(2*g))*n_elbow(k));
    
    h_L_minor_90elbow = h_L_minor_90elbow+ h_L_minor_90_th(k);
    
    k = k+1;
    
end

clear h_L_minor_90_th n_elbow

%last diameter for nozzle fitting
%only one nozzle fitting

h_L_minor_nozzle = ((K_L_nozzle*(Velocity(n-1).^2))/(2*g));

h_L_minor = h_L_minor_globe+h_L_minor_nozzle+h_L_minor_90elbow+h_L_minor_contraction;

h_L = h_L_minor+h_L_major;


%CAVITATION AND NSPH CHECK

for k = 1:2;
    
    Pressure(k) = -10;
    
    while Pressure(k) < 0;
        
        prompt = ['Enter Pressure ',num2str(k),' in pounds per square inch gage:'];
        
        Pres = inputdlg(prompt,'Pumped',[1,20]);
        
        Pressure(k) = (str2num(Pres{:}))*144; %change to pounds per square feet gage
        
    end
    
    k = k+1;
    
    clc
    
end

clear Pres

deltaz1 = -1;

while deltaz1 < 0;
    
    prompt = ['Enter change in height between the reservoir and the nozzle (z1-z2) in feet'];
    
    delta = inputdlg(prompt,'Pumped',[1,20]);
    
    deltaz1 = (str2num(delta{:})); %change to pounds per square feet gage
    
end

clear delta

velocity_no = n-1;

%Energy equation

h_P = (h_L+((Pressure(2)-Pressure(1))/(density*g))+(((Velocity(velocity_no)^2)-(Velocity(1)^2))/(2*g))+deltaz1);

deltaz2 = -1;

while deltaz2 < 0;
    
    prompt = ['Enter change in height between the reservoir and the pump (z1-z2) in feet'];
    
    delta = inputdlg(prompt,'Pumped',[1,20]);
    
    deltaz2 = (str2num(delta{:})); %change to pounds per square feet gage
    
end

clear delta

P_liquid = ((((Pressure(1)/(density*g))+deltaz2+h_P - h_L)+15.6959494)*144); %converting to psfa


%Check for cavitation

if P_Vapour <= P_liquid
    
    Cavitation_pipe = 0;
    
    prompt_pipe = ('NO cavitation in pipe');
    
else
    
    Cavitation_pipe = 1;
    
    prompt_pipe = ('YES, cavitation in pipe');
    
end


%Check for NPSH

NPSH_F = -1;

while NPSH_F <= 0;
    
    prompt = ('Enter the factory Net Positive Suction Head: ');
    
    NPSH2 = inputdlg(prompt,'Pumped',[1,20]);
    
    NPSH_F = (str2num(NPSH2{:}));
    
end

NPSH_C = (P_liquid/(density*g))+(((Velocity(1)).^2)/(2*g))-((P_Vapour)/(density*g));

if NPSH_C > NPSH_F
    
    Cavitation_pump = 0;
    
    prompt_pump = ('NO cavitation in pump');
    
else
    
    Cavitation_pump = 1;
    
    prompt_pump = ('YES, cavitation in pump');
    
end


%CLOSING MESSAGE

waitfor(msgbox('Thank you for using Pump Calculator 2016! The results from this calculation have been printed to both the command window and an Excel sheet located on the desktop.', 'Pump Calculator 2016'));


%PRINT RESULTS TO DISPLAY

%User Inputs

fprintf('USER INPUTS : \n \n');

fprintf('Diameters : \n');

for k = 1:n
    
    prompt = sprintf('Diameter %d (ft)',k);
    
    dia(k,1) = {prompt};
    
    dia(k,2) = {diameter(k)};
    
    k=k+1;
    
end

disp(dia);

fprintf('\n Lengths for each Diameter:  \n')

for k = 1:n
    
    prompt = sprintf('Length of Diameter %d (ft)',k);
    
    len(k,1) = {prompt};
    
    len(k,2) = {length(k)};
    
    k = k+1;
    
end

disp(len);

fprintf('\n Pressure initial and final: \n')

for k = 1:2
    
    prompt = sprintf('Pressure %d (psfg)',k);
    
    pres(k,1) = {prompt};
    
    pres(k,2) = {Pressure(k)};
    
    k = k+1;
    
end

disp(pres);

fprintf('\n Height change: \n');

prompt = sprintf('Reservoir and nozzle(ft)');

delta(1,1) = {prompt};

delta(1,2) = {deltaz1};

prompt = sprintf('Reservoir and pump(ft)');

delta(2,1) = {prompt};

delta(2,2) = {deltaz2};

disp(delta);

fprintf('\n \n');

hi={'User Inputs',' '};

Input = [hi;dia;len;pres;delta];

clear hi dia len pres delta

%OUTPUTS

fprintf('OUTPUTS: \n ');

fprintf('Reynold''s number for each velocity with cavitation concern:  \n')

prompt = sprintf('Diameter No');
Output(1,1) = {prompt};

prompt = sprintf('Velocity (ft/s)');
Output(1,2) = {prompt};

prompt = sprintf('Reynold''s Number');
Output(1,3) = {prompt};

prompt = sprintf('Laminar');
Output(1,4) = {prompt};

prompt = sprintf('Turbulent');
Output(1,5) = {prompt};

prompt = sprintf('f_f');
Output(1,6) = {prompt};

m = n+1;
l = 1;

for k = 2:m
    
    Output(k,1) = {l};
    
    Output(k,2) = {Velocity(l)};
    
    Output(k,3) = {R_N(l)};
    
    Output(k,4) = {category_laminar{l}};
    
    Output(k,5) = {category_turbulent{l}};
    
    Output(k,6) = {f_f(l)};
    
    k = k+1;
    
    l = l+1;
      
    end
    
prompt = sprintf('Calculated NPSH'); 

Output{1,7} = prompt;

Output{2,7} = NPSH_C;

prompt = sprintf('Factory NPSH'); %put inputs here for comparision

Output{1,8} = prompt;

Output{2,8} = NPSH_F;

prompt = sprintf('Cavitation in pipe'); 

Output{3,7} = prompt;

Output{3,7} = prompt_pipe;

prompt = sprintf('Cavitation in pump'); 

Output{3,8} = prompt;

Output{3,8} =  prompt_pump;

prompt = sprintf('Total Head Loss (ft)'); 

Output{4,7} = prompt;

Output{4,8} = h_L;

prompt = sprintf('Total Head of Pump (ft)'); 

Output{5,7} = prompt;

Output{5,8} = h_P;

prompt = sprintf('Volumetric Flow Rate (gpm)');

Output{1,9} = prompt;

Output{2,9} = VFL;

disp(Output);

%PRINT RESULTS TO FILE
%Volumetric Flowrate in gpm
%Reynolds Number
%Laminar or Turbulent Flow
%Friction Factor
%Pump Head (Print it as just head, not pump head)
%Cavitation concern (in pipe and in pump)
%NPSH at the pump inlet

xlswrite(filename,Input,inputsheet);

xlswrite(filename,Output,outputsheet);