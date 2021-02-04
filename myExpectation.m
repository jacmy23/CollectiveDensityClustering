function [deltaE] = myExpectation(indi, indj, v)

diri = v(indi, :);
dirj = v(indj, :);

diri = diri ./ norm(diri);
dirj = dirj ./ norm(dirj);

if size(diri, 1) > 1
    diri = sum(diri) / length(diri);
end
Ei = acos(diri(1) / norm(diri(:)));
if (diri(2) < 0)
    Ei = 2 * pi - Ei;
end

if size(dirj, 1) > 1
    dirj = sum(dirj) / length(dirj);
end
Ej = acos(dirj(1) / norm(dirj(:)));
if (dirj(2) < 0)
    Ej = 2 * pi - Ej;
end

result = abs(Ei - Ej);
if result > pi
    result = 2 * pi - result;
end

% result
if result <= pi / 4
    deltaE = 1;
else
    deltaE = 0;
end
end