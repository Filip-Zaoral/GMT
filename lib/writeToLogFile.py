
# GMT - Copyright (C) 2022 Filip Zaoral, IT4Innovations,
#                          VSB-Technical University of Ostrava, Czech Republic

# This file is a part of GMT.

# See the LICENSE.txt file in the GMT root directory for license information.

# GMT is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.

# GMT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

def write(Log,name):
    
    # Write the message console outputs to the Log File:
    logFile = open(name + '.log','wt')
    logFile.writelines('\n'.join(Log))
    logFile.close()                                                            # The saved file can be located in the working directory.
    return
