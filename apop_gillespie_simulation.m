%Apoptosis Extrinsic Pathway (Gillespie Stochastic Simulation)

P8 = 30; %(inactive) procaspase 8 molecules
C8 = 30; %(active) caspase 8 molecules
P3 = 200; % (inactive) procaspase 3 molecules
C3 = 0; %(active) caspase 3 molecules
IAP = 0; % XIAP molecules (inhibitor)
iC3I= 0;

%reaction rates (sec^-1)

k1 = 5.8*10^(-4);
k2 = 10^(-4);
k3= 5*10^(-3);
k4= 0.035;
k5 = 3*10^(-3);

tgraph = zeros(1,501);
C8graph = zeros(1,501);
C3graph = zeros(1,501);
IAPgraph = zeros(1,501);

r=1;

for i=1:30

t = 0;
while t<=r %run time up to r sec; 1 simulation
 
    
a1 = k1*C8*P3 ; %reaction 1 of 2 possible reactions
a2 = k2*C3*P8;
a3 = k3*C3*IAP;
a4 = k4*iC3I;
a5 = k5*C3*IAP;

a0 = a1 + a2 + a3 + a4 + a5;

r1 = rand;
r2 = rand;
tau = (-1/a0)*log(r1);
t = t + tau;

if (r2 >= 0) && (r2 < (a1/a0))
    P3 = P3 -1;
    C3 = C3 +1;
    
elseif (r2 >= (a1/a0)) && (r2 < ((a1 + a2)/a0))
    P8 = P8 -1;
    C8 = C8 + 1;
    
elseif (r2 >= ((a1 + a2)/a0)) && (r2 < ((a1 + a2+ a3)/a0))
    C3 = C3 - 1;
    IAP = IAP - 1;
    iC3I = iC3I + 1;

elseif (r2 >= ((a1 + a2+ a3)/a0)) && (r2 < ((a1 + a2+ a3 + a4)/a0))
    iC3I = iC3I -1;
    IAP = IAP + 1;
    C3 = C3 + 1;
     
else 
   IAP = IAP - 1; 
      
end


end
r= r+ 1;
C3graph(i+1) = C3;
C8graph(i+1) = C8;
IAPgraph(i+1) = IAP;
tgraph(i+1) = t;
end
plot(C3graph)
xlabel('time (sec)')
xlim([0 30])
ylabel('Caspase3 concentration (nM)')



