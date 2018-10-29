clc
close all
figure(1)
hold on
time_mat = [];
for i = 24:47
    plot(i, res{1,i}.test_save(end,5),'s');
end

figure(2)
hold on
for i = 24:1:47
    plot(i, res{1,i}.max_time,'*');
    time_mat = [time_mat;i, res{1,i}.max_time];
end
time_mat(size(time_mat,1), 3) = 1;
time_mat(size(time_mat,1), 4) = 0.88;
for i = 1:size(time_mat,1)-1
    time_mat(i, 3) = time_mat(size(time_mat,1), 2)/time_mat(i,2);
    time_mat(i, 4) = time_mat(size(time_mat,1), 4)/time_mat(i,3);
end