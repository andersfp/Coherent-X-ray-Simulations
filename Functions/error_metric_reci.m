function [chi_square_last]=error_metric_reci(diff_avant,diff_apres,chi_square)

chi_square_new=(1/(length(diff_avant)).^2).*sum(sum(sum((abs(diff_avant-diff_apres)).^2)));
chi_square_last=[chi_square chi_square_new];
