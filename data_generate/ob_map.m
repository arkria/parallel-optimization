function ob_near = ob_map( x, y, obstacle )
%MAP_CON �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
    ob_near.x = 0.5*(obstacle(1,1)+obstacle(2,1));
    ob_near.y = 0.5*(obstacle(1,2)+obstacle(4,2));

end

