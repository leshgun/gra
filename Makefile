#
# This file is part of a graduate thesis to modernize the Galbraith-Rupray 
# algorithm on equivalence classes instead of the exponential part of 
# the Gaudry-Schost algorithm for counting points of a curve of genus 2.
#
# Developed for the Immanuel Kant Baltic Federal University.
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#


FLAGS = -O4 ##-O4 -Wall
#FLAGS = -g

LIBS = NTLJac2/libNTLJac2.a -lntl -lgmp -lpthread #../pml/pml.a
#LIBS = NTLJac2/libNTLJac2.a /usr/local/lib/libntl.a -lgmp -lpthread
#LIBS = NTLJac2/libNTLJac2.a ../NTL/libNTLext.a -lntl -lgmp -lpthread
#LIBS = NTLJac2/libNTLJac2.a ../NTL/libNTLext.a /localdisk/gaudry/lib/libntl.a /localdisk/gaudry/lib/libgmp.a

INC = -I. -I NTLJac2 -I.. #-I../pml/include

CXX = g++ -std=c++11

objects = elemQuo.o resultantsumroots.o ZZ_pXResultant.o elltorsion.o schoof.o fastinterp.o vec_pair_long_long.o ZZ_pXCantorpoly.o 2ktorsion.o 4torsion.o memusage.o geometric_interpolation.o 2ktorsion2.o 4torsion2.o s1s2_heuristic.o bsgs_mct.o

path = ""
pd = ""
nt = "t50"
hs = "hs1000"
qs = "qs64"

all: $(objects) libschoof.a gen_curve GRA GRAT

%.o: %.c %.h
	$(CXX) $(FLAGS) $(INC) -c $<

libschoof.a: $(objects)
	ar r libschoof.a $(objects)

gen_curve: gen_curve.c $(objects)
	$(CXX) -o gen_curve gen_curve.c $(FLAGS) $(INC) libschoof.a $(LIBS)

GRA: GRA.c $(objects)
	$(CXX) -o $(path)GRA$(nt)$(hs)$(qs)$(pd) GRA.c $(FLAGS) $(INC) libschoof.a $(LIBS)

GRAT: GRAT.c $(objects)
	$(CXX) -o $(path)GRAT$(nt)$(hs)$(qs)$(pd) GRAT.c $(FLAGS) $(INC) libschoof.a $(LIBS)

clean:
	rm -f $(objects) libschoof.a main bsgs_mct_test LMPMCT


