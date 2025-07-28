function [ result ] = objj( Chrom )

% AckleysFunction 10
% for i = 1:size(Chrom,1)
%     x = Chrom;
%     result(i) =  -20 * exp(-0.2 * sqrt(1/size(x,2) * sum(x(i,:).^2)))-exp(1/size(x,2) * sum(cos(2 * pi * x(i,:))))+20+exp(1);
% end


% % Axis parallel hyper-ellipsoid function
% for i = 1:size(Chrom,1)
%     x = Chrom;
%     result(i) = sum([1:size(x,2)].*x(i,:).^2);
% end

% % % Rastrigin's function 6
for i = 1:size(Chrom,1)
    x = Chrom;
    n = size(x,2);
    result(i) =  10 * n + sum(x(i,:).^2-10 * cos(2 * pi * x(i,:)));

end



end


