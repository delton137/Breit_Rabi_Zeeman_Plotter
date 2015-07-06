%Calculates the energy levels and transition frequencies from the Breit-Rabi
%formula for Rb-85 and Rb-87. (this only applies to j=1/2 states)
%     Copyright (C) 2012 Daniel C. Elton 
% 
%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation; either version 2 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License along
%     with this program; if not, write to the Free Software Foundation, Inc.,
%     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

%Parameters
Z = 87; %Isotope 5/2 for Rb-85 and 3/2 for Rb-87 

mu = 9.274009*10^(-24); %Bohr magneton 
h = 6.626068*10^(-34); 

if Z == 85
   I = 5/2;
   W = h*3035.735*10^6;
   g_I = -.000293640;
end
if Z == 87
   I = 3/2;
   W = h*6834.685*10^6;
   g_I = -.0009951414;
end

g_l = 1;
g_s = 2.002;
S = 1/2;
L = 0;
J = 1/2;

g_J = g_l*(J*(J+1) + L*(L+1) - S*(S+1))/(2*J*(J+1)) + g_s*(J*(J+1) - L*(L+1) + S*(S+1))/(2*J*(J+1));


%%% Make plot of the Breit-Rabi equation for different H for both I+1/2 and I-1/2 
%%% Not the most elegant code, but it works
Bfield = linspace(0,.8,5000);

cases = [1/2,-1/2];

for j = 1:2
    F = I + cases(j);
    g_F = g_J*(F*(F+1) + J*(J+1) - I*(I+1))/(2*F*(F+1)) + g_I*(F*(F+1) + I*(I+1) - J*(J+1))/(2*F*(F+1));
   
    Energies = zeros(2*F+1,length(Bfield));
    m = linspace(-F,F,2*F+1);
          
    for i = 1:length(Bfield)
    
        x = (g_J - g_I)*mu*Bfield(i)/W;
    
        first = -W/(2*(2*I+1)) + mu*g_I*Bfield(i)*m; %first part of Breit-Rabi
        second = (W/2)*sqrt(1 + 4*m*x/(2*I + 1) + x^2); %second part of Breit-Rabi
        
        if j == 1
        second(1) = (W/2)*(1-x);
        end
        
        if j == 1
        E2_plus(:,i) = first + second;
        end
        if j == 2
        E2_minus(:,i) = first - second;
        end
    end
end

freq_shift_plus = E2_plus/h*10^-6;
freq_shift_minus = E2_minus/h*10^-6;

clf;
     plot(Bfield,freq_shift_plus,'k','LineWidth',1.5)
     hold on;
     plot(Bfield,freq_shift_minus,'k','LineWidth',1.5)
     set(gca,'YTick',linspace(-1.5*10^4,1.5*10^4,13),'LineWidth',1.5,'FontSize',15)
     set(gca,'XTick',linspace(0,.8,17),'LineWidth',1.5,'FontSize',15)
     xlabel('Magnetic Field (Tesla)','FontSize',25);
     ylabel('Freq Shift (MHz)','FontSize',25);
     title('Zeeman Splitting vs. Magnetic Field for Rb-87','FontSize',30)
 
