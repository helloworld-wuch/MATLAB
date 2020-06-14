function Endmember = test_p(x)
[Lines, Columns, Bands]=size(x); 
N = Lines * Columns;    % 像元数
Data = reshape(x, N, Bands);
Data = permute(Data, [2,1]);
numEndmember = 145;   % 端元数
numSkewers = 145;   % 投影向量数
skewers = randn(Bands, numSkewers);  % 投影向量
votes = zeros(N, 1);    % 纯像元指数
for i = 1:numSkewers
    r = skewers(:,i).' * Data;
    [val, id] = max(r);
    votes(id) = votes(id) + 1;
    [val, id] = min(r);
    votes(id) = votes(id) + 1;
end
[val, id] = sort(votes, 'descend');

Endmember = Data(:, id(1:numEndmember));
Endmember =sum(Endmember);
return;