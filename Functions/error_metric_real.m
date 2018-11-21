function [chi_square_R_last]=error_metric_real(objet_avant,objet_apres,chi_square)

chi_square_R=sum(sum(sum((abs(objet_avant-objet_apres)).^2)));
chi_square_R_last=[chi_square chi_square_R];
