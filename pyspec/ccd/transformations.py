#
# transformations.py (c) Stuart B. Wilkins 2010
#
# $Id$
# $HeadURL$
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Part of the "pyspec" package
#

import numpy

class CCDExperiment:
    detectorSampleDistance = 100 #mm

    def __init__():
        return

    def getTransformationMatrix():
        return numpy.eye(3)

class CCDTardis(CCDExperiment):
    detectorSampleDistance = 300 # mm


class CCDTransformation:
    def __init__(ccd = None):
        self.ccd = ccd
