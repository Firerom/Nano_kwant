%% 
clear all
close all
Squarre_var_E_fixed_W_L=dlmread('/Users/YOLO/Desktop/1_master_complementaire/cours/Q1/nanoelectronics/TP/Squarre_var_E_fixed_W_L.txt');
Squarre_var_L_fixed_E02_W=dlmread('/Users/YOLO/Desktop/1_master_complementaire/cours/Q1/nanoelectronics/TP/Squarre_var_L_fixed_E02_W.txt');
Squarre_var_L_fixed_E05_W=dlmread('/Users/YOLO/Desktop/1_master_complementaire/cours/Q1/nanoelectronics/TP/Squarre_var_L_fixed_E05_W.txt');
Squarre_var_L_fixed_E08_W=dlmread('/Users/YOLO/Desktop/1_master_complementaire/cours/Q1/nanoelectronics/TP/Squarre_var_L_fixed_E08_W.txt');

Squarre_var_W_fixed_E02_L=dlmread('/Users/YOLO/Desktop/1_master_complementaire/cours/Q1/nanoelectronics/TP/Squarre_var_W_fixed_E02_L.txt');
Squarre_var_W_fixed_E05_L=dlmread('/Users/YOLO/Desktop/1_master_complementaire/cours/Q1/nanoelectronics/TP/Squarre_var_W_fixed_E05_L.txt');
Squarre_var_W_fixed_E08_L=dlmread('/Users/YOLO/Desktop/1_master_complementaire/cours/Q1/nanoelectronics/TP/Squarre_var_W_fixed_E08_L.txt');
honey_var_E=dlmread('/Users/YOLO/Desktop/1_master_complementaire/cours/Q1/nanoelectronics/TP/honey_var_E.txt');


figure('units','normalized','outerposition',[0 0 1 1],'visible','on')
plot(Squarre_var_E_fixed_W_L(:,1),Squarre_var_E_fixed_W_L(:,2), 'linewidth', 5)
hold on 
xlabel('E','interpreter', 'latex','FontSize',30)
ylabel('T', 'interpreter', 'latex','FontSize',30)
set(gca,'fontsize',30)
set(gcf,'color','w')
grid on


figure('units','normalized','outerposition',[0 0 1 1],'visible','on')
plot(Squarre_var_L_fixed_E02_W(:,1),Squarre_var_L_fixed_E02_W(:,2), 'linewidth', 5,'color','r')
hold on 
plot(Squarre_var_L_fixed_E05_W(:,1),Squarre_var_L_fixed_E05_W(:,2), 'linewidth', 5,'color','g')

plot(Squarre_var_L_fixed_E08_W(:,1),Squarre_var_L_fixed_E08_W(:,2), 'linewidth', 5)

xlabel('L','interpreter', 'latex','FontSize',30)
ylabel('T', 'interpreter', 'latex','FontSize',30)
set(gca,'fontsize',30)
set(gcf,'color','w')
grid on
h=legend('$E=0.2$','$E=0.5$','$E=0.8$');
set(h,'interpreter', 'latex')

figure('units','normalized','outerposition',[0 0 1 1],'visible','on')
plot(Squarre_var_W_fixed_E02_L(:,1),Squarre_var_W_fixed_E02_L(:,2), 'linewidth', 5,'color','r')
hold on 
plot(Squarre_var_W_fixed_E05_L(:,1),Squarre_var_W_fixed_E05_L(:,2), 'linewidth', 5,'color','g')

plot(Squarre_var_W_fixed_E08_L(:,1),Squarre_var_W_fixed_E08_L(:,2), 'linewidth', 5)

xlabel('W','interpreter', 'latex','FontSize',30)
ylabel('T', 'interpreter', 'latex','FontSize',30)
set(gca,'fontsize',30)
set(gcf,'color','w')
grid on
h=legend('$E=0.2$','$E=0.5$','$E=0.8$');
set(h,'interpreter', 'latex')



figure('units','normalized','outerposition',[0 0 1 1],'visible','on')
plot(honey_var_E(:,1),honey_var_E(:,2), 'linewidth', 5)
hold on 
xlabel('E','interpreter', 'latex','FontSize',23)
ylabel('T', 'interpreter', 'latex','FontSize',23)
set(gca,'fontsize',20)
set(gcf,'color','w')
grid on
