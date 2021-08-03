% script_cooperativity
%
% This script is to investigate a question asked by one of the reviewers:
% what are the effects of cooperativity between X and W or Y and W?
%
% We will create a separate m-file function that looks at cooperativity
% between X and W, and Y and W. First: the work on C13 (cooperativity
% between X and W) and then on C23 (cooperativity between Y and W).
%
% These cooperativities, at most, can change the width of the NPA domain in
% the FF/FB case. So all we have to do is re-calculate Delta C_FFFB and
% produce a set of heatmaps for different values of C13 and C23 that are
% similar to Fig 3E and F and put these in the supplement.
%
% Here is the goal (do this for C13 first, then C23 second): Set C13 = [0.5
% 0.1 0.05 0.01] and plot the corresp. Fig3E,F. Do this again for C23.
% 


clear
close all
options = ddeset('RelTol',1e-6);


%% ========================================================================
% First, we will calculate the Delta(K3)/K3 for the situation with negative
% feedback. This is happening first, even though it's like part E,F, or
% something, because A,B are hand-drawn, and C,D will use parameters that
% come from this K1-K2 sweep.
% =========================================================================
% {

clear
load Mat/DeltaK_FFLonly_dense K11_dense K22_dense DeltaK
load Mat/C_pointone_contour XC YC VC
load Mat/Peaks_pointone_contours_FFFB XP YP Xzmax Yzmax VP Vzmax
nK = 500;
K11 = logspace(-2,6,nK)'; K111 = repmat(K11',nK,1);
K22 = logspace(-6,2,nK)';
K4 = 0.1; K8 = 0.1;

%
% creating K3_NPAplus and minus arrays
%
% script_Fig3_FBL_K3NPA_K12sweep
load Mat/Fig3_FBL_K3NPA_K12sweep


% %
% % Plotting DeltaK3^{+}/K3
% %
Y = (real(K3_NPA_FBL(:,:,1)) - K3_PA)./K3_PA; 
% pcolor_contour(K11,K22,Y,[0.15 0.3 1 5 10],'r',true,0,5)
% plot(K1_I_FBL,K2_I_FBL,'r','linewidth',2)
% plot(K1_II_FBL,K2_II_FBL,'r','linewidth',2)
% title('\Delta K_3^{NPA,FBL,+}/K_3^{PA}')
% 
% 
% %
% % Plotting DeltaK3^{-}/K3
% %
Yminus = (K3_PA - real(K3_NPA_FBL(:,:,2)))./K3_PA;
% % pcolor_contour(K11,K22,log10(K3_NPA_FBL(:,:,2)),[0 0.15 0.3 1],'r',true,-6,-0.7)
% pcolor_contour(K11,K22,Yminus,[0.15 0.3 1 5 10],'r',true)
% plot(K1_I_FBL,K2_I_FBL,'r','linewidth',2)
% plot(K1_II_FBL,K2_II_FBL,'r','linewidth',2)
% title('\Delta K_3^{NPA,FBL,-}/K_3^{PA}')
% 
% 
% %
% % Plotting DeltaK3^{tot}/K3
% %
DeltaK_FBL = Y + Yminus;
pcolor_contour(K11,K22,DeltaK_FBL,[0.15 0.3 1 5],'r',true,0,5)
set(gca,'xtick',10.^([-2 0 2 4 6]),'ytick',10.^([-6 -4 -2 0 2]))
% plot([K1_I;K1upper],[K2_I;max(K22)],'r','linewidth',2)
% plot([K1lower;K1_II],[min(K22);K2_II],'r','linewidth',2)
plot(K1_I_FBL,K2_I_FBL,'r','linewidth',2)
plot(K1_II_FBL,K2_II_FBL,'r','linewidth',2)
plot([0.1 0.1 10 10 0.1],[0.1 0.5 0.5 0.1 0.1],':k','linewidth',2)
plot(XC,YC,'w',XP,YP,'w',Xzmax,Yzmax,'w')
% savepcolor(gcf,'Figs/DeltaK3_NPA_FFFB')
xlabel('K1')
ylabel('K2')
title('\Delta K_3^{NPA,FBL}/K_3^{PA}')


DeltaK1 = interp2(repmat(K11_dense',length(K22_dense),1),...
	repmat(K22_dense,1,length(K11_dense)),DeltaK,...
	repmat(K11',nK,1),repmat(K22,1,nK));
% pcolor_contour(K11,K22,DeltaK1,[0.15 0.3 1 5 10],'r',true,0,5)
% set(gca,'xtick',10.^([-2 0 2 4 6]),'ytick',10.^([-6 -4 -2 0 2]))
% plot([K1_I;K1upper],[K2_I;max(K22)],'r','linewidth',2)
% plot([K1lower;K1_II],[min(K22);K2_II],'r','linewidth',2)
% plot([0.1 0.1 10 10 0.1],[0.1 0.5 0.5 0.1 0.1],':k','linewidth',2)
% xlabel('K1')
% ylabel('K2')
% title('\Delta K_3^{NPA}/K_3^{PA}')



%
% pcolor and contour of the ratio of the two DeltaK's. Note that in this
% contour, there is a weird curve for the ratio equal to 1.5.
%
RDK = DeltaK_FBL./DeltaK1; % ratio of DeltaK's
pcolor_contour(K11,K22,RDK,[1 1.1 1.25 1.5 1.75 2],'r',true,0.99,2)
set(gca,'xtick',10.^([-2 0 2 4 6]),'ytick',10.^([-6 -4 -2 0 2]))
plot(K1_I_FBL,K2_I_FBL,'r','linewidth',2)
plot(K1_II_FBL,K2_II_FBL,'r','linewidth',2)
% plot([K1_I;K1upper],[K2_I;max(K22)],'r','linewidth',2)
% plot([K1lower;K1_II],[min(K22);K2_II],'r','linewidth',2)
plot([0.1 0.1 10 10 0.1],[0.1 0.5 0.5 0.1 0.1],':k','linewidth',2)
plot(XC,YC,'w',XP,YP,'w',Xzmax,Yzmax,'w')
h = get(gcf,'children');
for i = 1:length(h)
	if strcmp(h(i).Tag,'Colorbar')
		cbh = h(i);
		set(h(i),'YTick',1:0.25:2)
	end
end
% savepcolor(gcf,'Figs/Ratio_of_DeltaK3')
xlabel('K1')
ylabel('K2')
title('\Delta K_3^{NPA,FBL}/\Delta K_3^{NPA}')

%
% Create mathematical definition of realistic region
%
V = VP & VC & Vzmax;
min(RDK(V))
max(RDK(V))

%}



%% ========================================================================
% Now that we have found the values of K3 for which f  = eps = 0.05 (and
% for f = -0.05, we will run a set of simulations for a particular choice
% of K1 and K2 to show what it looks like.
% =========================================================================
%{

clear
options = ddeset('RelTol',1e-6);
load Mat/Fig3_FBL_K3NPA_K12sweep
K1 = 1; K2 = 0.1; 

K3_NPAplus = interp2(repmat(K11',nK,1),repmat(K22,1,nK),K3_NPA_FBL(:,:,1),K1,K2);
K3_NPAminus = interp2(repmat(K11',nK,1),repmat(K22,1,nK),K3_NPA_FBL(:,:,2),K1,K2);

%
% K3_PA
%
F0 = 1;
F1 = 10;
y0 = F0./(K1 + F0);
y1 = F1./(K1 + F1);
RHS = y1 - y0;
LHS = (1/F0-1/F1) - (1./K2)*(1./(K1 + F1) - 1./(K1 + F0));
K3 = RHS/LHS;



%
% Disturbance function
%
t0 = 0;
disthand = @(t)ftn_unitstep(t,t0,F0,F1);

%
% Param values
%
n = 1;
thetaz = 0.5; thetaw = 0.5;
theta = [0.5 thetaz thetaw];
tauz = 1; tauw = 1; tauy = 1;
tspan = [0 20];
K4 = 0.1; K8 = 0.1;
p = [tauz K1 K2 K3 K4 K8 tauw tauy];
K4 = 1; K8 = Inf;
p1 = [tauz K1 K2 K3 K4 K8 tauw tauy];


%
% Perfect adaptation FF/FB sim
%
ftnhand = @ftn_goentoro_dde;
[y0,z0,w0] = ftn_goentoro_ss(p,F0,n);
Y0 = [y0; z0; w0];
soln = dde23(ftnhand,theta,Y0,tspan,options,theta,p,n,disthand);
t = soln.x';
Y = soln.y';
z = Y(:,2);

% Then FF only sim
[y0,z0FF,w0] = ftn_goentoro_ss(p1,F0,n);
Y0 = [y0; z0FF; w0];
soln = dde23(ftnhand,theta,Y0,tspan,options,theta,p1,n,disthand);
tFF = soln.x';
Y = soln.y';
zFF = Y(:,2);


%
% Near perfect adaptation (plus), FF/FB sim
%
p(4) = K3_NPAplus;
[y0,z0plus,w0] = ftn_goentoro_ss(p,F0,n);
Y0 = [y0; z0plus; w0];
soln = dde23(ftnhand,theta,Y0,tspan,options,theta,p,n,disthand);
tplus = soln.x';
Y = soln.y';
zplus = Y(:,2);

% Then FF only sim
p1(4) = K3_NPAplus;
[y0,z0FFplus,w0] = ftn_goentoro_ss(p1,F0,n);
Y0 = [y0; z0FFplus; w0];
soln = dde23(ftnhand,theta,Y0,tspan,options,theta,p1,n,disthand);
tFFplus = soln.x';
Y = soln.y';
zFFplus = Y(:,2);


%
% Near perfect adaptation (plus), FF only sim
%
p(4) = K3_NPAminus;
[y0,z0minus,w0] = ftn_goentoro_ss(p,F0,n);
Y0 = [y0; z0minus; w0];
soln = dde23(ftnhand,theta,Y0,tspan,options,theta,p,n,disthand);
tminus = soln.x';
Y = soln.y';
zminus = Y(:,2);

% Then FF only sim
p1(4) = K3_NPAminus;
[y0,z0FFminus,w0] = ftn_goentoro_ss(p1,F0,n);
Y0 = [y0; z0FFminus; w0];
soln = dde23(ftnhand,theta,Y0,tspan,options,theta,p1,n,disthand);
tFFminus = soln.x';
Y = soln.y';
zFFminus = Y(:,2);


%
% Plotting
%
figure
plot(t,z,tplus,zplus,tminus,zminus,'linewidth',2)
set(gca,'fontsize',24)
print(gcf,'Figs/NPA_simulation_FFFB.eps','-depsc')
print(gcf,'Figs/NPA_simulation_FFFB.jpg','-djpeg','-r150')

figure
plot(t,z/z0,tplus,zplus/z0plus,tminus,zminus/z0minus,'linewidth',2)
hold on
plot(xlim,[1 1]*1.05,'k:',xlim,[1 1]*0.95,'k:')
set(gca,'fontsize',24)
print(gcf,'Figs/NPA_simulation_FFFB_norm.eps','-depsc')
print(gcf,'Figs/NPA_simulation_FFFB_norm.jpg','-djpeg','-r150')

plot(tFF,zFF/z0FF,tFFplus,zFFplus/z0FFplus,tFFminus,zFFminus/z0FFminus,'linewidth',2)
print(gcf,'Figs/NPA_simulation_FFFB_norm1.eps','-depsc')
print(gcf,'Figs/NPA_simulation_FFFB_norm1.jpg','-djpeg','-r150')



%}



%% ========================================================================
% Now we run a parameter sweep to get the peak value parameter P while
% varying values of K1,K2, at K3_PA, for FFL only. It turns out that this
% parameter is very weak for Region I. So that's a no-go region.
% =========================================================================
% {
nK = 200;
K11 = logspace(-2,6,nK)';
K22 = logspace(-6,2,nK)';
F0 = 1;
F1 = 10;

%
% K3_PA
%
y0 = F0./(K11 + F0);
y1 = F1./(K11 + F1);
y00 = repmat(y0',nK,1);
y11 = repmat(y1',nK,1);
% RHS = y11 - y00;
% LHS = (1 + (1./K22)*y0')/F0 - (1 + (1./K22)*y1')/F1;
% K3_PA = RHS./LHS;

%
% Disturbance function
%
t0 = 0;
disthand = @(t)ftn_unitstep(t,t0,F0,F1);

%
% Param values
%
n = 1;
thetaz = 0.5; thetaw = 0.5;
theta = [0.5 thetaz thetaw];
K1 = 1; K2 = 0.3; K3 = 0.3; K4 = 0.1; K8 = 0.1;
K5 = K1*K4; K6 = K2*K4; K7 = K3*K4;
tauz = 1; tauw = 1; tauy = 1;
tspan = [0 100];

%
% All (slow) simulation commands in the following script, which is
% commented out now, since the sims were already run and stored in the Mat
% file below.
%
% script_Fig3_FFFB_dde_K12sweep
load Mat/Fig3_peak_K12sweep_FFFB

%
% Boundary of Region I (and II, too)
%
load Mat/Fig3_FBL_K3NPA_K12sweep K1_I_FBL K2_I_FBL K1_II_FBL K2_II_FBL 

%
% Contour of C_PA = 0.1
%
load Mat/C_pointone_contour XC YC VC


%
% Plotting pcolor and contour of P.
%
P = (Zpeak - Zinitial)./Zinitial;
pcolor_contour(K11,K22,P,[0.01 0.1 1 1.95 2.05 3],'r',true)
set(gca,'xtick',10.^([-2 0 2 4 6]),'ytick',10.^([-6 -4 -2 0 2]))
plot(K1_I_FBL,K2_I_FBL,'r','linewidth',2)
plot(K1_II_FBL,K2_II_FBL,'r','linewidth',2)
plot([0.1 0.1 10 10 0.1],[0.1 0.5 0.5 0.1 0.1],':k','linewidth',2)
% savepcolor(gcf,'Figs/Peak_surface_FFFB')
xlabel('K1')
ylabel('K2')
title('P = (Z_{max}-Z_0)/Z_0')

%
% Get the 0.1 contour
%
h = get(gca,'children');
for i = 1:length(h)
	if strcmp(h(i).DisplayName,'0.1')
		XP = get(h(i),'Xdata')';
		YP = get(h(i),'Ydata')';
		break
	end
end
nK2 = 500;
K11_dense2 = logspace(-2,6,nK2)';
K22_dense2 = logspace(-6,2,nK2)';
P2 = interp2(repmat(K11',nK,1),repmat(K22,1,nK),P,...
	repmat(K11_dense2',nK2,1),repmat(K22_dense2,1,nK2));
VP = P2 >= 0.1;

%
% Plotting pcolor and contour of Zpeak. This is for Fig. S2.
%
pcolor_contour(K11,K22,Zpeak,[0.01 0.1 1],'r',true)
set(gca,'xtick',10.^([-2 0 2 4 6]),'ytick',10.^([-6 -4 -2 0 2]))
plot(XP,YP,'w')
plot(XC,YC,'w')
plot(K1_I_FBL,K2_I_FBL,'k','linewidth',1)
plot(K1_II_FBL,K2_II_FBL,'k','linewidth',1)
plot([0.1 0.1 10 10 0.1],[0.1 0.5 0.5 0.1 0.1],':k','linewidth',2)
% savepcolor(gcf,'Figs/Zmax_surface_FFFB')
xlabel('K1')
ylabel('K2')
title('Z_{max}')

%
% Get the zmax = 0.1 contour
%
h = get(gca,'children');
for i = 1:length(h)
	if strcmp(h(i).DisplayName,'0.01')
		Xzmax = get(h(i),'Xdata')';
		Yzmax = get(h(i),'Ydata')';
		break
	end
end
nK2 = 500;
K11_dense2 = logspace(-2,6,nK2)';
K22_dense2 = logspace(-6,2,nK2)';
Zpeak2 = interp2(repmat(K11',nK,1),repmat(K22,1,nK),Zpeak,...
	repmat(K11_dense2',nK2,1),repmat(K22_dense2,1,nK2));
Vzmax = Zpeak2 >= 0.01;

%
% Plotting pcolor and contour of Z0.
%
pcolor_contour(K11,K22,Zinitial,[0.01 0.1 1],'r',true)
set(gca,'xtick',10.^([-2 0 2 4 6]),'ytick',10.^([-6 -4 -2 0 2]))
plot(XP,YP,'w')
plot(K1_I_FBL,K2_I_FBL,'r','linewidth',2)
plot(K1_II_FBL,K2_II_FBL,'r','linewidth',2)
plot([0.1 0.1 10 10 0.1],[0.1 0.5 0.5 0.1 0.1],':k','linewidth',2)
% savepcolor(gcf,'Figs/Z0_surface_FFFB')
xlabel('K1')
ylabel('K2')
title('Z_{0}')


% save Mat/Peaks_pointone_contours_FFFB XP YP Xzmax Yzmax VP Vzmax

%}









