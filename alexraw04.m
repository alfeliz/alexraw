###########################################################################################
#
#
#  OCTAVE SCRIPT TO:
#       - READ ALEX DATA FROM TEKTRONIK SCOPES,
#	- TRANSFORM THESE DATA INTO ELECTRICAL ENERGY DATA.
#    Made by Gonzalo Rodríguez Prieto
#       (gonzalo#rprietoATuclm#es)
#       (Mail: change "#" by "." and "AT" by "@")
#              Version 3.60
#
#
#########################################################
#
#  It uses the functions: 
#                        baseline	To find the zero on experimental data that could be biased. 
#                        chan		To extract voltage and time data from Tektronik made *.CVS files.
#                        supsmu		To smooth scope data through a filter.
#  They must be in the same directory.
#
###########################################################################################

more off; %To make the lines display when they appear in the script, not at the end of it.

clear; %Just in case there is some data in memory.

tic; %Total time of the script.




#########################################################################################################
#
#	SCOPE DATA INTERPRETATION PART.
#
#########################################################################################################



###
# Reading information over the shot features:
###

#Choosing the shot:
shot = "ALEX098"; %TO make it unique

#Reading the TXT file with info over the shot:
#Shot name transformation:
infoshot = horzcat(shot,".txt");
#Open it:
[fileinfo, msg] = fopen(infoshot, "r");
if (fileinfo == -1) 
   error ("Script alexraw: Unable to open file name: %s, %s",infoshot, msg); 
endif;
#Putting the information into a cell structure:
txtinfo = textscan(fileinfo,"%s", "treatasempty","/");
fclose(fileinfo); #Closing the file
txtinfo = cell2mat(txtinfo); #To convert the cell structure in a matrix of strings.

#Putting the channels info in matrices of strings:
#(This is valid if the format of file SHOT.txt follows the procedure of ALEX data)
channels = ""; #Channel name
chantype = ""; #Channel probe type
for i=1 : length(txtinfo)
 if (cell2mat(strfind(txtinfo(i),"CH")) ==1) %We have the channels here, with the type following.
   % Standard format of the *.txt file for ALEX.
   channels = strvcat(channels,txtinfo(i));
   chantype = strvcat(chantype,txtinfo(i+1));
 elseif(cell2mat(strfind(txtinfo(i),"Del")) ==1) %We should stop to avoid channels confusion.
   #After channels info cames other shot info, not stored or used in this version of ALEX.m
   break;
 endif;
endfor;


###
#Taking all channel data in one matrix (vol): 
###

for i=1 : rows(channels)
 channame = strcat(shot,".",channels(i,:),".CSV");
 fprintf(stdout,"Taking data from %s\n", channame);
 [t(:,i),vol(:,i),points(i)] = chan(channame); %chan is the function to take the channel file CHANNAME info.
 # It is already smoothed.
endfor;


###
#Transforming scope data from (t,volts) into (t,dimensional units) taking into account the type:
###

#Selecting the variable type from the file information:
selection = []; %Defining the variable
for i=1 : rows(channels)
#Choosing which adaptation is necessary from the channel type characters:
if strcmp(chantype(i,1:4),"Div.")==1
 selection = [selection,1]; %Number 1: My resistive divider
endif;
if strcmp(chantype(i,1:5),"Sonda")==1
 selection = [selection,2]; %Number 2: Tektronik voltage probe
endif;
if strcmp(chantype(i,1:3),"Rog")==1
 selection = [selection,3]; %Number 3: Rogowsky coil probe
endif;
if ( length(chantype(i,:))>8 )
    if strcmp(chantype(i,1:10),"Photodiode")==1
     selection = [selection,4]; %Number 4: Photodiode signal
    endif;
endif;
endfor;

# Data transformation of the channels from previous selection:
for i=1 : columns(selection)
 switch(selection(i))
    case 1 %My resistive divider data
      volt(:,2) = vol(:,i) ./ 0.00032; %Wire volts
      volt(:,1) = t(:,i).*1e6; %time in us
      volt(:,2) = supsmu(volt(:,1),volt(:,2),"span",0.01); %Smoothing the voltage data.
	#It allows for a better treatment of the data, without losing too much information.
      disp("Transforming voltage data from resistive divider");
    case 2 %The tektronik probe (if this exists, not the one up, and viceversa)
      volt(:,2) = vol(:,i); %Wire volts. No transformation, but not clear if correct.
      volt(:,1) = t(:,i).*1e6; %time in us
      volt(:,2) = supsmu(volt(:,1),volt(:,2),"span",0.01); %Smoothing
      disp("Transforming voltage data from Tektronik probe");
    case 3 %Rogowsky coil
      #Due to integration problems, I will try now to make the conversion AFTER the integration
      dint(:,2) = vol(:,i); 
        %Following 
      dint(:,1) = t(:,i).*1e6; %time in us
      dint(:,2) = supsmu(dint(:,1),dint(:,2),"span",0.01); %Smoothing
      disp("Transforming Rogowsky signal");
    case 4 %Photodiode signal (Normalized to 1)
      phot(:,2) = vol(:,i);
      phot(:,2) = phot(:,2)./max(phot(:,2)); 
      phot(:,1) = t(:,i).*1e6; %time in us
      phot(:,2) = supsmu(phot(:,1),phot(:,2),"span",0.02); %Smoothing
      disp("Transforming photodiode signal");
 endswitch
endfor;

#Checking and warning for missed data on the shot:
if (exist("volt","var")==0)
  warning("Script alex: There is no voltage signal in shot %s",shot);
elseif (exist("dint","var")==0)
  warning("Script alex: There is no intensity signal in shot %s",shot);
elseif (exist("phot","var")==0)
  warning("Script alex: There is no photodiode signal in shot %s",shot);
endif;


###
#Removing values on the baseline of signals for integration or further operations of the data:
###

#Putting the voltage baseline to zero:
disp("Removing the baselines from the scope data.");
if (exist("volt","var")==1)
  [basev,mediav,puntosv] = baseline(volt(1:250,2)); %Voltage signal
  if (puntosv<=10)
    warning("Script alexraw: only %u points in voltage baseline",puntosv);
  endif;
  volt(:,2) = volt(:,2) - mediav;
endif;

#Photodiode signal to zero:
disp("Removing the baselines from the photodiode data.");
if (exist("phot","var")==1)
  [basep,mediap,puntosp] = baseline(phot(1:200,2)); %Photodiode signal
  if (puntosp<=10)
    warning("Script alexraw: only %u points in photodiode baseline",puntosp);
  endif;
  phot(:,2) = phot(:,2) - mediap;
endif;

#Derivative intensity signal to zero:
disp("Removing the baselines from the d(i)/dt data.");
if (exist("dint","var")==1)
  [basep,mediap,puntosp] = baseline(dint(1:200,2)); %d(i)/dt signal
  if (puntosp<=10)
    warning("Script alexraw: only %u points in d(i)/dt baseline",puntosp);
  endif;
  dint(:,2) = dint(:,2) - mediap;
endif;


###
# Integrating the intensity derivative:
###
if (exist("dint","var")==1)
  int = - cumtrapz(dint(:,1)./1e6,dint(:,2)) .*(63.095/5.85e-9) ; %The minus because of the signal in the voltage and the derivative.
  #The division of 1e6 is because we need to use ALL unit in I.S!!!!!
  # The other factor (63.095/5.85e-9) follows the calibration page 33 of "Logbook of diverse projects 01" for the d(i)/dt signal.
  dint(:,2) = dint(:,2) .*(63.095/5.85e-9); %To transform d(i)/dt to use it later.
endif;


###
# Removing the inductive part on the voltage to obtain only the resistive voltage:
###

#First it is necessary to obtain the induction part of the wire, by checking the difference on the period on the voltage signal from short circuit case.
C = 2.27e-6; %Farads. Total capacity of the circuit. Measured and written in page 21 of "Logbook Exploding Wire".

[vma,posma] = max(volt(:,2)); %Find the value and position of the first voltage peak.
[vmi,posmi] = min(volt(:,2)); %Find the value and position of the second voltage peak.
period = abs(posma-posmi)*abs(volt(1,1)-volt(2,1))*2; %(µs) Period of the wire and the other part of the circuit.
Lt = ( (period*1e-6)/(2 * pi))**2 * (1/C) - 142e-9 %Henrios. Total inductance of the circuit minus short circuit inductance: Wire inductance.

Ltot = Lt + 65e-9 ;
Ltot * 1e9 %To check in the screen the vaue of the inductance


Vres(:,1) = volt(:,1); %Time vector is the same.
Vres(:,2) = volt(:,2) + ( Ltot .* dint(:,2) );%Volts, resistive part of the wire.
#Positive sign on the inductive part because of the Rogowsky sensor alignment!!!!

plot(Vres(:,1),Vres(:,2)) %To see the voltage through the wire as purely resistive.

####
## Finding the electrical power and energy delivered by time: 
####

if ( (exist("int","var")==1) && (exist("volt","var")==1) && (exist("Vres","var")==1) ) 
  elpowr2 = Vres(:,2).*int; %Wire electrical power (Watts)
  elpowr = volt(:,2).*int; %Total circuit electrical power (Watts)
  electricenergia(:,1) = dint(:,1); %Time (microseconds)
  electricenergia2(:,1) = dint(:,1); %Time (microseconds)
  electricenergia(:,2) = cumtrapz(elpowr).*(abs(dint(1,1)-dint(2,1)).*1e-6); %Energy (Joules) delivered to the total circuit
  electricenergia2(:,2) = cumtrapz(elpowr2).*(abs(dint(1,1)-dint(2,1)).*1e-6); %Energy (Joules) delivered to the wire
endif;





#############################################################################################################
#####
#####	SAVING AND PLOTTING CALCULATED DATA
#####
#############################################################################################################


#####
###Saving the data in files:
#####

disp("Creating and saving data files.");
#Voltage (if there is)
if (exist("Vres","var"))
  #Output file name:
  name = horzcat(shot,"_voltage.txt"); %Adding the right sufix to the shot name.
  output = fopen(name,"w"); %Opening the file.
  #First line:
  fdisp(output,"descriptor `time(micros)`  `voltage(V)` `tot. volt(V)`"); %Veusz format
  redond = [2 3 3]; %Saved precision 
  vol = [Vres, volt(:,2)];
  display_rounded_matrix(vol, redond, output); %This function is not made by my.
  fclose(output); %Closing the file.
  disp("Voltage saved.");
endif;

#Intensity (if there is)
if (exist("dint","var"))
  #Output file name:
  name = horzcat(shot,"_intensity.txt"); %Adding the right sufix to the shot name.
  output = fopen(name,"w"); %Opening the file.
  #First line:
  fdisp(output,"descriptor `time(micros)`  `intensity(A)`"); %Veusz format
  redond = [2 3]; %Saved precision 
  intensity = [dint(:,1),int];
  display_rounded_matrix(intensity, redond, output); 
  fclose(output); %Closign the file.
  disp("Intensity saved.");
endif;

#Photodiode signal (if there is)
if (exist("phot","var"))
  #Output file name:
  name = horzcat(shot,"_photodiode.txt"); %Adding the right sufix to the shot name.
  output = fopen(name,"w"); %Opening the file.
  #First line:
  fdisp(output,"descriptor `time(micros)`  `luminosity(A.U.)`"); %Veusz format
  redond = [2 2]; %Saved precision 
  display_rounded_matrix(phot, redond, output); %This function is not made by my.
  fclose(output); %Closing the file.
  disp("Relative luminosity saved.");
endif;

#Electrical energy (if there is)
if (exist("electricenergia2","var"))
  #Output file name:
  name = horzcat(shot,"_wire_energy.txt"); %Adding the right sufix to the shot name.
  output = fopen(name,"w"); %Opening the file.
  #First line:
  fdisp(output,"descriptor `time(micros)`  `wire energy(J)` `total energy(J)`"); %Veusz format
  redond = [2 3 3]; %Saved precision 
  energia = [electricenergia2, electricenergia(:,2)]; %Saving total and wire energies.
  display_rounded_matrix(energia, redond, output); %This function is not made by my.
  fclose(output);%Closing file
  disp("Electrical energy saved.");
endif;

#Resistance
if (exist("int","var") && exist("volt","var"))
  #Output file name:
  name = horzcat(shot,"_resistance.txt"); %Adding the right sufix to the shot name.
  output = fopen(name,"w"); %Opening the file.
  #First line:
  fdisp(output,"descriptor `time(micros)`  `res(Ohm)`"); %Veusz format
  redond = [2 3]; %Saved precision
  resistance = [dint(:,1),Vres(:,2)./int]; 
  display_rounded_matrix(resistance, redond, output); %This function is not made by my.
  fclose(output);%Closing file
  disp("Resistance saved.");
endif;


more on; #Revert more normal comments behaviour.

###
# Total processing time
###
timing = toc;
disp("Script alexraw execution time:")
disp(timing)
disp(" seconds")

#That's...that's all, folks! 
