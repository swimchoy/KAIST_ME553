ga_true = [-39.5612
-5.39841
 6.32672
 63.7426
 25.6895
 13.9361];
gf = [0.1, 0.2, 0.3, 0.4, 0.4, 0.3]';
a_raisim_0 = zeros(6,1);
a_raisim_1 = [        0
          0
          0
8.18217e-12
          0
    39.5612];
a_raisim_2 = [-0.0632979
   -1.6e-05
1.30915e-14
      -0.02
    5.39841
    39.5612];
a_raisim_3 =[ -2.26764
  1.60297
-0.204648
     0.01
  11.7251
  39.5612];
a_raisim_4 =[ 0.601424
   2.4177
-0.449368
 -6.39345
  11.7291
  -23.859];
a_raisim_5 =[ 1.81162
   2.8254
-0.571846
 -6.64245
 -13.9603
  -23.839];
a_raisim_6 =[1.30398
0.715958
0.795293
-19.7526
-13.6853
-28.5716];
M_6 = [       1.337            0            0           -0   -0.0262576 -7.66874e-14
           0        1.337            0    0.0262576           -0   -0.0768585
           0            0        1.337  7.66874e-14    0.0768585           -0
           0    0.0262576  7.66874e-14    0.0112362  4.22781e-15  -0.00144763
  -0.0262576            0    0.0768585  4.22792e-15    0.0154735  1.44442e-15
-7.66874e-14   -0.0768585            0  -0.00144763  1.44462e-15     0.014979];
b_6 = [   0.0380576
   -0.025951
  -0.0184222
-0.000488788
 -0.00173247
  0.00143073];
S6 = [          0           0           0     -0.9463 9.44193e-13    -0.32329]';

M_5 = [       0.463            0            0           -0  -0.00744237   0.00132246
           0        0.463            0   0.00744237           -0   -0.0217846
           0            0        0.463  -0.00132246    0.0217846           -0
           0   0.00744237  -0.00132246  0.000251494 -6.22227e-05 -0.000454043
 -0.00744237            0    0.0217846 -6.22227e-05   0.00157674 -2.12575e-05
  0.00132246   -0.0217846            0 -0.000454043 -2.12575e-05    0.0014254];
b_5 = [  0.00871909
 -0.00150041
 0.000505987
-3.42914e-05
-0.000147282
 0.000121048];
S5 = [           0            0            0 -1.18451e-12           -1 -1.45995e-11]';

X65 = [           1            0            0            0            0            0
            0            1            0            0            0            0
           0            0            1            0            0            0
           0   -0.0335413 -6.05981e-13            1            0            0
   0.0335413            0   -0.0981786            0            1            0
 6.05981e-13    0.0981786            0            0            0            1];