clear
clc
load('Compare_n324s127.mat')
Neigh_dist = 25;

index_vect = [10, 3, 6,7, 22, 2, 18, 8, 34, 1, 30, 9, 11, 15, 5, 19, 23, ...
    14, 17, 20, 35, 13, 29, 21, 12, 27, 4, 31, 24, 26, 16, 32, 36, 25, 28, 33];
index_vect = [index_vect, index_vect+36, index_vect+72,...
    index_vect+108, index_vect+144, index_vect+180,...
    index_vect+216, index_vect+252, index_vect+288];
cen_time = [];
den_time = [];
step_save = [];
net_save = []; 
end_num = 324;
epsi_rel = 0.01;
epsi_abs = 0.01;
rho = 20000;
max_step = 10000;

for veh_n = 10
    veh_cell = cell(1,veh_n);
    for i = 1:veh_n
        veh_cell{1,i} = Veh_cell{1,index_vect(i)};
    end

    veh_couple_mat = ADMM_coupleCheck(veh_cell, Neigh_dist);
    [f,g,rel_mat] = ADMM_transfer2(veh_cell, veh_couple_mat, obstacle);
    [f1,g1,rel_mat1] = ADMM_transfer(veh_cell, veh_couple_mat, obstacle);
    [delta_temp1, empty_tag, MPC_time] = Central_write(f1, g1, rel_mat1);
    w_num = size(f, 2);
    m_num = size(g, 2);
    file_w = fopen('worker.txt','w');
    for i = 1:w_num
        [m,n] = size(f{1,i}.A_o);
        for j = 1:n
            for k = 1:n
                fprintf(file_w, '%20.4f', f{1,i}.H_o(j,k));
            end
        end
        fprintf(file_w, '\n');
        for k = 1:n
            fprintf(file_w, '%20.4f', f{1,i}.f_o(k));
        end
        fprintf(file_w, '\n');
        for j = 1:m
            for k = 1:n
                fprintf(file_w, '%20.4f', f{1,i}.A_o(j,k));
            end
        end
        fprintf(file_w, '\n');
        for k = 1:m
            fprintf(file_w, '%20.4f', f{1,i}.b_o(k));
        end
        fprintf(file_w, '\n');
        for k = 1:n
            fprintf(file_w, '%20.4f', f{1,i}.lb_o(k));
        end
        fprintf(file_w, '\n');
        for k = 1:n
            fprintf(file_w, '%20.4f', f{1,i}.ub_o(k));
        end
        for k = 1:n
            fprintf(file_w, '%20.4f', f{1,i}.z_temp(k));
        end
        fprintf(file_w, '\n');
    end

    fclose(file_w);
    file_m = fopen('master.txt','w');
    for i = 1:m_num
        m = size(f{1,i}.H_o, 1);
        for k = 1:m
            fprintf(file_m, '%20.4f', g{1,i}.z_temp(k));
        end
        fprintf(file_w, '\n');
    end
    fclose(file_m);
    
    file_o = fopen('other.txt', 'w');
    [W, M] = size(rel_mat);
    fprintf(file_o, '%10d', M);
    fprintf(file_o, '%10d', W);
    fprintf(file_o, '\n');
    for i = 1:W
        c = find(rel_mat(i,:) == 1);
        for j = 1:size(c, 2)
            fprintf(file_o, '%10d', c(j)-1);
        end
    end
    fprintf(file_o, '\n');
    for i = 1:M
        s = size(find(rel_mat(:,1:i) == 1), 1);
        fprintf(file_o, '%10d', s);
    end
    fprintf(file_o, '\n');
    temp = zeros(1, sum(sum(rel_mat)));
    tag = 1;
    for i = 1:M
        c = find(rel_mat(:,i) == 1);
        for j = 1:size(c, 1)
            fprintf(file_o, '%10d', c(j)-1);
            temp(tag) = c(j)-1;
            tag = tag + 1;
        end
    end
    fprintf(file_o, '\n');
    dict = zeros(W, 1);
    for i = 1:size(temp,2)
        if dict(temp(i)+1) == 0
            dict(temp(i)+1) = 1;
            fprintf(file_o, '%10d', 0);
        else
            fprintf(file_o, '%10d', 1);
        end
    end
    fprintf(file_o, '\n');
    fclose(file_o);
    
end