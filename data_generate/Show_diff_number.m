clear
clc
close all;
load('asy_diff_num.mat')
% figure(1)

for i = 2:6
    target = res{1,i};
%     plot(1:2:target.ADMM_step_sy-1, log10(target.test_save(1:2:end-1,5)),'linewidth', 2);
    semilogy(1:2:target.ADMM_step_sy-1, target.test_save(1:2:end-1,5),'linewidth', 2);
    hold on
end
set(gca, 'FontName', 'Times New Roman','Fontsize',15);
xlabel('Global iteration step','Fontname', 'Times New Roman','FontSize',15);
ylabel('Maximum error [%]','Fontname', 'Times New Roman','FontSize',15);
legend('10','15','20','25','30', 'fontsize',50,'fontname','Times New Roman','location','northeast')
box on