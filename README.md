# Listing irreducible polynomials modulo n

`main.f90` is a Fortran 2008 program that lists all irreducible polynomials of degree at most `D` with coefficients modulo `n`. An irreducible polynomial is nonzero, has no multiplicative inverse, and cannot be written as a product of two polynomials that both lack multiplicative inverses. The listing includes irreducible constants and different polynomials that differ by multiplication by a unit (an invertible polynomial).

## Build and run

From the folder containing `main.f90`, compile with GNU Fortran and run the program. The following commands use Windows PowerShell.

```powershell
gfortran -std=f2008 -O3 main.f90 -o list_irreducibles.exe
.\list_irreducibles.exe
```

At the prompt, enter `n D` on one line. The modulus `n` must be between 2 and 2,147,483,647, and `D` must be nonnegative. For example, `4 1` produces one irreducible constant (`2`) and eight irreducible polynomials of degree 1.

The program writes the list and counts to `irreducibles_Z<n>_degree<D>.txt` in the current folder, replacing any file with that name. It also prints the count for each degree and the total in the terminal.

## How the algorithm works

1. **Split the modulus.** Factor `n` into prime powers `q = p^k`. The Chinese remainder theorem lets the program work modulo each `q` separately and then combine the results coefficient by coefficient.
2. **Screen monic polynomials modulo `p`.** A monic polynomial has leading coefficient 1. For each degree, the program enumerates monic polynomials modulo `p` and uses polynomial greatest common divisors to check whether each is a power `h^m` of a single irreducible polynomial `h`. A reduction with different irreducible factors can be split into two coprime factors, which lift to a factorization modulo `q`.
3. **Check the lifts modulo `q`.** A lift is a polynomial modulo `q` that reduces to the given polynomial modulo `p`. If `m = 1`, every monic lift is irreducible. Otherwise, the program tests monic lifts of `h^j` for `1 ≤ j ≤ floor(m/2)`. Any factorization has a smaller or equal degree factor in that range. If a candidate divides the polynomial with zero remainder modulo `q`, the program discards that polynomial.
4. **Generate the full listing.** The program multiplies the irreducible monic polynomials by units and, when `k > 1`, separately generates the irreducible zero divisors, which are `p` times a unit modulo `q`. It chooses one irreducible component and a unit in every other component, then recombines them to obtain an irreducible polynomial modulo `n`.

The program writes each result as it finds it. 
