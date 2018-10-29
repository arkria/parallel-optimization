% clear
% clc
% % close all;
% load('asy_res_5_10.mat')
% figure(2)
% hold on
% for i = 2:2:10
%     target = res{1,i};
%     plot(1:2:target.ADMM_step_sy-1, log10(target.test_save(1:2:end-1,5)),'linewidth', 2);
% end
% set(gca, 'FontName', 'Times New Roman','Fontsize',15);
% xlabel('Iteration step','Fontname', 'Times New Roman','FontSize',15);
% ylabel('Maximum error [log(%)]','Fontname', 'Times New Roman','FontSize',15);
% legend('2','4','6','8','10','fontsize',50,'fontname','Times New Roman','location','northeast')
% box on

clear
clc
close all;
load('asy_res_20_47.mat')

for i = 7:8:47
    target = res{1,i};
    semilogy(1:2:target.ADMM_step_sy-1, target.test_save(1:2:end-1,5),'linewidth', 2);
    hold on
end
set(gca, 'FontName', 'Times New Roman','Fontsize',15);
xlabel('Iteration step','Fontname', 'Times New Roman','FontSize',15);
ylabel('Maximum error [%]','Fontname', 'Times New Roman','FontSize',15);
legend('7','15','23','31','39','47','fontsize',50,'fontname','Times New Roman','location','northeast')
box on