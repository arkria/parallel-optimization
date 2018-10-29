function out_mat = Other_boundFun( in_mat, bound )
%BOUNDFUN �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
    [m,n] = size(in_mat);
    out_mat = zeros(m,n);
    for i = 1:m
        for j = 1:n
            if in_mat(i,j) >= bound(i,2)
                out_mat(i,j) = bound(i,2);
            elseif in_mat(i,j) <= bound(i,1)
                out_mat(i,j) = bound(i,1);
            else
                out_mat(i,j) = in_mat(i,j);
            end
        end
    end


end

