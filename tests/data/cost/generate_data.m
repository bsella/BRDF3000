name_file = 'costSet3';

num_rows = 5;
d = 60;

#num_rows
dlmwrite(name_file, num_rows);

#num_cols
dlmwrite(name_file, d, ' ', 1, "-append");

#Z
Z = rand(num_rows, d);
dlmwrite(name_file, Z, ' ', 1, "-append", "precision", "%.100f");

#K_minus1
do
	random = rand(num_rows, num_rows) ./ 10;
	K = (random + random') .* 5;
	detK = det(K);
until (abs(detK) > 0.00000000001 &&  !iscomplex(log(detK)));

K_minus1 = inverse(K);

dlmwrite(name_file, K_minus1, ' ', 1, "-append", "precision", "%.100f");

#detK
dlmwrite(name_file, detK, ' ', 1, "-append", "precision", "%.100f");

name_file = [name_file, "_output"];

#cost
cost = 0.5 * d * log(detK) + 0.5 * trace(K_minus1 * Z * transpose(Z));

dlmwrite(name_file, cost, ' ', 1, "precision", "%.100f");