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


set grid x y

set xlabel "Stochastic reconfiguration iterations"

set ylabel "Variational parameters"

plot \
'opt_vpar_hist.txt' using 2 with lines title "t'", \
'opt_vpar_hist.txt' using 3 with lines title "t''", \
'opt_vpar_hist.txt' using 4 with lines title "Delta_onsite", \
'opt_vpar_hist.txt' using 5 with lines title "Delta", \
'opt_vpar_hist.txt' using 6 with lines title "Delta'", \
'opt_vpar_hist.txt' using 7 with lines title "Delta''", \
'opt_vpar_hist.txt' using 8 with lines title "mu", \
'opt_vpar_hist.txt' using 9 with lines title "mu_m"

pause -1

plot \
'opt_vpar_hist.txt' using 10 with lines title "onsite Jastrow"

pause -1

set ylabel "Energy"

plot \
'opt_E_hist.txt' using 2 with lines title "E"

pause -1
