/*
 * Copyright (c) 2013, Robert Rueger <rueger@itp.uni-frankfurt.de>
 *
 * This file is part of hVMC.
 *
 * hVMC is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * hVMC is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with hVMC.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef MSGTAGS_H_INCLUDED
#define MSGTAGS_H_INCLUDED

enum msgtag_t {
  // convention:
  // MSGTAG_SENDER_RECEIVER_CONTENT
  MSGTAG_S_M_REQUEST_BINS,
  MSGTAG_M_S_DISPATCHED_BINS,
  MSGTAG_S_M_FINISHED_BINS
};

enum schedmsg_t {
  SCHEDMSG_START_MCC,
  SCHEDMSG_EXIT
};

#endif // MSGTAGS_H_INCLUDED
