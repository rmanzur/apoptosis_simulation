%Euler Method: Extrinsic Pathway


%reaction rates (sec^-1)
k1 = 5.8*10^(-4);
k2 = 10^(-4);
k3= 5*10^(-3);
k4= 0.035;
k5 = 3*10^(-3);

%initial conditions
C8= [30,30,30,30];
P8= [30,45,60,90];
C3= [0,0,0,0];
P3= [200,200,200,200] ;
IAP =[30,30,30,30];
iC3I = [0,0,0,0];

dt = 1;
tgraph = zeros(4,1001);
C8graph = zeros(4,1001);
P8graph = zeros(4,1001);
P3graph = zeros(4,1001);
C3graph = zeros(4,1001);
IAPgraph = zeros(4,1001);

for i = 1:4
for t = 0:1000
    

C8(i) = C8(i) + (dt*k2*P8(i)*C3(i));
P8(i) = P8(i) + (dt*-k2*P8(i)*C3(i));
C3(i) = C3(i) + (dt*(k1*C8(i)*P3(i)-k3*C3(i)*IAP(i)+ k4*iC3I(i)));
P3(i) = P3(i) + (dt*-k1*C8(i)*P3(i));
IAP(i) = IAP(i)+ (dt*(-k3*IAP(i)*C3(i)+ k4*iC3I(i)-k5*IAP(i)*C3(i)));

C8graph(i,t+1) = C8(i);
P8graph(i,t+1) = P8(i);
P3graph(i,t+1) = P3(i);
C3graph(i,t+1) = C3(i);
IAPgraph(i,t+1) = IAP(i);
tgraph(i,t+1) = t;

end

end

figure
x =plot (C3graph(1,:))
hold on;
y=plot (C3graph(2,:))
hold on;
z=plot (C3graph(3,:))
hold on;
m=plot (C3graph(4,:))
xlabel('time (sec)')
xlim([0 250])
ylim([0 200])
set([x y z m],'LineWidth',1)
ylabel('Caspase 3 concentration (nM)')
legend('[ProCaspase 8] = baseline','[ProCaspase 8] = + 50%','[ProCaspase 8] = + 100%','[ProCaspase 8] = + 200%')
figure
x =plot (P3graph(1,:))
hold on;
y=plot (P3graph(2,:))
hold on;
z=plot (P3graph(3,:))
hold on;
m=plot (P3graph(4,:))
xlabel('time (sec)')
xlim([0 250])
ylim([0 200])
set([x y z m],'LineWidth',1)
ylabel('ProCaspase 3 concentration (nM)')
legend('[ProCaspase 8] = baseline','[ProCaspase 8] = + 50%','[ProCaspase 8] = + 100%','[ProCaspase 8] = + 200%')

