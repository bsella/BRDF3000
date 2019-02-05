num_rows = 5;
index_changed = 3;

#num_rows
dlmwrite('shermanMorissonUpdateSet1', num_rows, ' ', "precision", "%.30f");

#lv_num
dlmwrite('shermanMorissonUpdateSet1', index_changed, ' ', 1, "-append", "precision", "%.30f");

#old_K_minus1
do
	random = rand(num_rows, 1) .* 10;
	K_minus1 = random * random';
	detK = det(K_minus1);
until (detK > 0.0001);

dlmwrite('shermanMorissonUpdateSet1', K_minus1, ' ', 1, "-append", "precision", "%.30f");

#old_detK
detK = det(K_minus1);

dlmwrite('shermanMorissonUpdateSet1', detK, ' ', 1, "-append", "precision", "%.30f");

#diff_cov_vector
new_values = rand(1, num_rows) .* 10;
diff = new_values - K_minus1(index_changed, :);

dlmwrite('shermanMorissonUpdateSet1', diff, ' ', 1, "-append", "precision", "%.30f");

#new_K_minus1
K_minus1(:, index_changed) = new_values';
K_minus1(index_changed, :) = new_values;

dlmwrite('shermanMorissonUpdateSet1', K_minus1, ' ', 1, "-append", "precision", "%.30f");

#new_detK
detK = det(K_minus1);

dlmwrite('shermanMorissonUpdateSet1', detK, ' ', 1, "-append", "precision", "%.30f");