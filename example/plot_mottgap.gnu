# Copyright (c) 2013, Robert Rueger <rueger@itp.uni-frankfurt.de>
#
# This file is part of hVMC.
#
# hVMC is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# hVMC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with hVMC.  If not, see <http://www.gnu.org/licenses/>.


set key top left

set xlabel '|q|'
set xrange [0:]

set ylabel 'E(|q|)'
set yrange [0:]

E(x) = exp(a)*x**2 + exp(b)*x + exp(logD)
fit E(x) 'ana_res_ssfac.txt' using (sqrt($1**2+$2**2)):(($1**2+$2**2)/$3):(($1**2+$2**2)**2) via a,b,logD

plot \
  'ana_res_ssfac.txt' using (sqrt($1**2+$2**2)):(($1**2+$2**2)/$3) title 'Data: |q|^2 / N(q)', \
   E(x) title 'Fit: E(|q|) = a |q|^2 + b |q| + D'

pause -1
