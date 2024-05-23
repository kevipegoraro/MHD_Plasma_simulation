function ip_control

%%
dtcharge = 40;
dtup = 20;
dttop = 60;
dtdown = 20;
Ip = 80e3;
dtdecharge = 40;

%%
t1= -dtcharge*1e-3;
t2 = 0;
t3 = dtup*1e-3;
t4 = (dtup + dttop)*1e-3;
t5 = (dtup + dttop + dtdown)*1e-3;
t6 = (dtup + dttop + dtdown + dtdecharge)*1e-3;

%% Machine parameters
Roh = 0.123e-3;
Loh = 11e-3;
M = 51e-6;
Lp = 1.4e-6;

% Solving equations
t = linspace(t1-10e-3,t6+10e-3,1001).';
Ip = interp1([t1 t2 t3 t4 t5 t6],[0 0 Ip Ip 0 0],t);
Ip = smooth(Ip,5);
Rp = 42./(Ip+1e3).^1.25;

dIpdt = dFdx(t,Ip);
dIohdt = -1/M*(Rp.*Ip + Lp*dIpdt);
Ioh = integration(t,dIohdt,0);
Ioh = Ioh + (max(Ioh) - min(Ioh))/2;
Ioh(t <= 0) = max(Ioh) + max(Ioh)/(t2 - t1)*t(t<=0);
Ioh(t >= t5) = min(Ioh(t>=t5 & t<=t5+1e-3)) - min(Ioh(t>=t5 & t<=t5+1e-3))/(t6-t5)*(t(t>=t5) - t5);
Ioh(t <= t1 | t >= t6) = 0;
dIohdt = dFdx(t,Ioh);
V = Roh*Ioh + Loh*dIohdt + M*dIpdt;
Rp(Rp>0.0009) = 1; % Just for plotting reasons

% Plotting
figure(1)
clf
plot(t*1000,Ip/100e3,'b','linewidth',2)
hold on
plot(t*1000,Ioh/1e3,'r','linewidth',2)
plot(t*1000,V/1e3,'k','linewidth',2)
plot(t*1000,Rp*1e3,'color',[0 0.5 0],'linewidth',2)
%plot(t*1000,M*dIohdt,'color',[1 1 1]*0.4,'linewidth',2)
hold off
xlabel('Time ( ms )')
legend('I_{Plasma} ( 100 kA )','I_{OH} ( kA )','V ( kV )','R_{plasma} ( m\Omega )')%,'V_{Loop}')
axis([t1*1000-10 t6*1000+10 -4 4])

figure(2)
clf
plot(Ioh/1e3,V/1e3,'b','linewidth',2)
hold on
plot([-6 6],[0 0],'k')
plot([0 0],[-4 4],'k')
hold off
xlabel('I_{OH} ( kA )')
ylabel('V ( kV )')
axis([-4 4 -3 3])

end