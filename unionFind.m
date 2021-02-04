function [id, sz] = unionFind(id, sz, a, b)
    [id, ra] = unionRoot(id, a);
    [id, rb] = unionRoot(id, b);
    
    if ra == rb
        return;
    end

    if sz(ra) < sz(rb)
        id(ra) = id(rb);
        sz(rb) = sz(rb) + sz(ra);
    else
        id(rb) = id(ra);
        sz(ra) = sz(ra) + sz(rb);
    end
end
