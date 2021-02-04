function [NCLUST, cl, union_id] = CDC(ND, dist, distd, xx, v)

percent=2.0;
position=round(size(xx, 1)*percent/100);

if position == 0
    NCLUST = 0; cl = 0; union_id = 0;
    return;
end

sda=sort(xx);
dc=sda(position);
dc_d = pi / 3;

% theta
% theta = exp(-(distd / dc_d).^2);
% theta(distd <= dc_d) = theta(distd <= dc_d) + 1;
% theta = theta.*(~eye(size(theta)));
% 
% rho
% adjM = dist < theta * dc;
% rho = sum(adjM)';
% 
% AdjMatrix
% AdjMatrix = dist < dc * 2;
% AdjMatrix = AdjMatrix.*(~eye(size(AdjMatrix)));

rho = zeros(ND, 1);
delta = zeros(ND, 1);
AdjMatrix = zeros(ND, ND);
% "Cut off" kernel
for i=1:ND-1
    for j=i+1:ND
        if(distd(i,j)>dc_d)
            tmp = exp( - ( distd(i,j) / dc_d )^2);
        else
            tmp = 1 + exp( - ( distd(i,j) / dc_d )^2);
        end
        
        %         tmp = 1 + exp( - ( distd(i,j) / dc_d )^2);
        
        if (dist(i,j) < dc * tmp)
            %rho(i)=rho(i)+1.;
            %rho(j)=rho(j)+1.;
            rho(i)=rho(i)+exp(-(dist(i,j)/dc)*(dist(i,j)/dc)/tmp);
            rho(j)=rho(j)+exp(-(dist(i,j)/dc)*(dist(i,j)/dc)/tmp);
        end
        
        if (dist(i,j) < dc * 2)
            AdjMatrix(i, j) = 1;
        end
    end
end

AdjMatrix = AdjMatrix + AdjMatrix';

maxd=max(max(dist));
[rho_sorted,ordrho]=sort(rho,'descend');
delta(ordrho(1))= -1.;
nneigh(ordrho(1))= ordrho(1);

for ii=2:ND
    delta(ordrho(ii))=maxd;
    for jj=1:ii-1
        if(dist(ordrho(ii),ordrho(jj))<delta(ordrho(ii)))
            delta(ordrho(ii))=dist(ordrho(ii),ordrho(jj));
            nneigh(ordrho(ii))=ordrho(jj);
        end
    end
end
delta(ordrho(1))=max(delta(:));
gamma = rho + delta;

cl = ones(ND, 1) * -1;
ind = find(gamma == max(gamma));
cl(ind) = ind;

%assignation
for i = 1 : ND
    if (cl(ordrho(i)) == -1 && distd(ordrho(i), nneigh(ordrho(i))) <= dc_d) %&& delta(ordrho(i)) < dc * theta(ordrho(i), nneigh(ordrho(i))))
        cl(ordrho(i)) = cl(nneigh(ordrho(i)));
    else
        cl(ordrho(i)) = ordrho(i);
        nneigh(ordrho(i)) = ordrho(i);
    end
end

u_cl = unique(cl);
for i = 1 : size(u_cl, 1)
    ind = find(cl == u_cl(i));
    if size(ind, 1) < 10
        u_cl(i) = -1;
        cl(ind) = -1;
    end
end

u_cl(u_cl == -1) = [];

NCLUST = size(u_cl, 1);
for i = 1 : NCLUST
    cl(cl == u_cl(i)) = i;
end

union_id = 1:NCLUST;
union_sz = ones(NCLUST, 1);

for i = 1 : NCLUST - 1
    for j = i + 1 : NCLUST
        adjIndx = find(cl == i);
        adjIndy = find(cl == j);
        [x, y] = find(AdjMatrix(cl == i, cl == j) > 0);
        if ~isempty(x)
            if union_id(i) == union_id(j)
                continue;
            end
            
            if myExpectation(adjIndx(unique(x)), adjIndy(unique(y)), v)
                [union_id, union_sz] = unionFind(union_id, union_sz, i, j);
            end
        end
    end
end

for i = 1 : NCLUST
    [union_id, union_id(i)] = unionRoot(union_id, union_id(i));
end

u_uid = unique(union_id);
for i = 1 : length(u_uid)
    union_id(union_id == u_uid(i)) = i;
end
end