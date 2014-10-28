/*
# Tailor, a BWT-based aligner for non-templated RNA tailing
# Copyright (C) 2014 Min-Te Chou, Bo W Han, Jui-Hung Hung
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#include <bitset>
#include <iostream>
#include <deque>
#include <string>
#include <algorithm>
#include <fstream>
#include "cstring"
#include "boost/dynamic_bitset.hpp"
#include "boost/utility/binary.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/algorithm/string.hpp"
#include "boost/random.hpp"
#include "../src/abwt.hpp"
#include "../src/abwt_search.hpp"
#include "../src/abwt_table.hpp"
#include "../src/abwt_format.hpp"
#include "../src/compression/abit.hpp"
#include "../src/compression/jbit.hpp"
#include "../src/tailer.hpp"
#include "../src/tailor_build.hpp"
#include "../src/tailor_map.hpp"
