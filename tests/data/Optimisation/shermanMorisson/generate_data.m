pkg load symbolic;
digits(200);

name_file = 'shermanMorissonUpdateSet2';

num_rows = 100;
index_changed = 10;

#num_rows
dlmwrite(name_file, num_rows);

#lv_num
dlmwrite(name_file, index_changed - 1, ' ', 1, "-append");

#old_K_minus1
do
	random = rand(num_rows, num_rows) ./ 10;
	K = (random + random') .* 5;
	detK = det(K);
until (abs(detK) > 0.00000000001);

dlmwrite(name_file, inverse(K), ' ', 1, "-append", "precision", "%.100f");

#old_detK
dlmwrite(name_file, detK, ' ', 1, "-append", "precision", "%.100f");

#diff_cov_vector
new_values = rand(1, num_rows) ./ 10;
diff = new_values - K(index_changed, :);

dlmwrite(name_file, diff, ' ', 1, "-append", "precision", "%.100f");

name_file = [name_file, "_output"];

#new_K_minus1
K(:, index_changed) = new_values';
K(index_changed, :) = new_values;

dlmwrite(name_file, inverse(K), ' ', "precision", "%.100f");

#new_detK
detK = det(K);

dlmwrite(name_file, detK, ' ', 1, "-append", "precision", "%.100f");
