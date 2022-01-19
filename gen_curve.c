/*
 * This file is part of a graduate thesis to modernize the Galbraith-Rupray 
 * algorithm on equivalence classes instead of the exponential part of 
 * the Gaudry-Schost algorithm for counting points of a curve of genus 2.
 *
 * Developed for the Immanuel Kant Baltic Federal University.
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "NTL/ZZ_p.h"
#include "NTL/ZZ_pX.h"

NTL_CLIENT

void randomizeF2(ZZ_pX& f)
{
    while(true)
    {
        ZZ_pX h;
        random(h, 6);
        if (deg(h) != 5) continue;
        MakeMonic(h);

        ZZ_p d = resultant(h, diff(h));
        if (d != 0)
        {
            SetCoeff(f, 5, coeff(h, 5));
            SetCoeff(f, 4, coeff(h, 4));
            SetCoeff(f, 3, coeff(h, 3));
            SetCoeff(f, 2, coeff(h, 2));
            SetCoeff(f, 1, coeff(h, 1));
            SetCoeff(f, 0, coeff(h, 0));
            break;
        }
    }
}

int main(int argc, char* argv[])
{
    long l;
    ZZ p;
    ZZ_pX f;
    
    if (argc == 1) 
    {
        cout << "Enter num bits of 'p' : ";
        cin >> l;
    }
    else l = stoi(argv[1]);

    GenPrime(p, l);
    cout << "p : " << p << endl;
    ZZ_p::init(p);

    randomizeF2(f);
    cout << "f : [ ";
    for (int i = 0; i < 6; ++i)
        cout << coeff(f, i) << ' ';
    cout << ']' << endl;

    return 0;
}
