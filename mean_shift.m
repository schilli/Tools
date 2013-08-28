
function mean_shift()

close all;

% parameters
nPoints = 2000;
scale   = 1;
SD      = 0.1;
eps     = 0.001;
maxIter = 1000000000;
normEps = 1e-12;
h       = SD;

% generate means
a = 1/sqrt(2);
means   = [[1 -a -a];
           [0  a -a]];
means = scale * means;

% generate points
points     = SD * randn(2,nPoints);
assignment = randi(size(means,2), 1, nPoints);
for i = 1:size(means,2)
    points(:,assignment == i) = bsxfun(@plus, ...
                                points(:,assignment == i), ...
                                means(:,i));
end


endPoints = zeros(size(points));

for i = 1:size(points,2)
    p = points(:,i);
    M = m(p, points, h);

    %fprintf('i = %d\n', i)
    %fflush(stdout);

    %fprintf('M = [%d, %d]\n', M(1), M(2));

    iter = 0;
    oldNorm = 2*norm(M);
    
    while norm(M) > eps && iter < maxIter && norm(M) ~= oldNorm && abs(norm(M) - oldNorm) > normEps
        iter = iter + 1;
        p = p + M;
        oldNorm = norm(M);
        M = m(p, points, h);

        %fprintf('norm(M) = %d, iter = %d\n', norm(M), iter);
        %fflush(stdout); 
    end
    
    if iter >= maxIter
        fprintf('Finished with iter reaching maxIter\n');
    end
    
    if norm(M) == oldNorm
        fprintf('Finished with norm(M) == oldNorm\n');
    end
    
    if abs(norm(M) - oldNorm) <= normEps
        fprintf('Finished with abs(norm(M) - oldNorm) <= normEps\n');
    end

    endPoints(:,i) = p;
end

P = find_convergence_points(endPoints, 1e-1)

% plot points
plot(points(1,:), points(2,:), 'bo', 'LineWidth', 2);
hold on;
plot(endPoints(1,:), endPoints(2,:), 'r*', 'LineWidth', 2);

end % mean_shift()


% ================================ %

function M = m(x, points, h)
    A = zeros(size(x));
    B = 0;
    for p = points
        A = A + p * g(norm(x-p)^2/h^2);
        B = B +     g(norm(x-p)^2/h^2);
    end

    M = (A/B - x);
end % m(x)

% ================================ %

function G = g(x)
    G = -2*x*exp(-x^2);
end

% ================================ %

function P = find_convergence_points(points, eps)

    P = [];
    done = zeros(size(points,2));
    
    for i = 1:size(points,2)
        if done(i) ~= 1
            p1 = points(:,i);
            P(:,end+1) = p1;
            
            for j = i+1:size(points,2)
                p2 = points(:,j);
                
                if norm(p1-p2) <= eps
                    done(j) = 1;
                end
            end
        end
    end
end
