function veh_couple_mat = ADMM_coupleCheck( veh_cell, Neigh_dist )
%ADMM_COUPLECHECK �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
    veh_n = size(veh_cell,2);
    veh_couple_mat = zeros(veh_n);
    for i = 1:veh_n-1
        for j = i+1:veh_n
            veh_M = veh_cell{1,i}; veh_N = veh_cell{1,j};
            if pdist([veh_M.x_now, veh_M.y_now; veh_N.x_now, veh_N.y_now]) <= Neigh_dist
                veh_couple_mat(i,j) = 1;
            end
        end
    end
    veh_couple_mat = veh_couple_mat+veh_couple_mat';
end

