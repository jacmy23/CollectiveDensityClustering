function [id, x] = unionRoot(id, x)
        while x ~= id(x)
            id(x) = id(id(x));
            x = id(x);
        end
end