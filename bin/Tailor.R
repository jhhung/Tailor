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

.libPaths(c(paste(Sys.getenv ("PIPELINE_DIRECTORY"),"Rlib",sep='/'), .libPaths()))

# this function is modified from http://stackoverflow.com/questions/9341635/how-can-i-check-for-installed-r-packages-before-running-install-packages
pkgTest <- function(x)
{
	if (!require(x,character.only = TRUE))
	{
		install.packages(x,dep=TRUE, lib=paste(Sys.getenv ("PIPELINE_DIRECTORY"),"Rlib",sep="/"), repos='http://cran.us.r-project.org')
		if(!require(x,character.only = TRUE)) stop ("Failed to install the package. Please check the internet access or update your R if it is too old.")
	}
}


