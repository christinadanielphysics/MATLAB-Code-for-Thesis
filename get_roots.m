function roots = get_roots(min,max,step,function_handle)
    roots = [];
    for x0 = min:step:max
        fun = function_handle; % e.g. @d 
        root = fzero(fun,x0);
        if ~ismember(round(root,4),round(roots,4))
            roots(end+1) = round(root,4);
        end
    end
end