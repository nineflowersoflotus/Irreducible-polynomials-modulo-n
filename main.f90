! List irreducible polynomials over Z/nZ through a given degree.
! Coefficients stored in ascending degree order.

module polynomials
  use iso_fortran_env, only: int64
  implicit none
contains
  ! The zero polynomial has degree -1 here.
  integer function degree(coefficients) result(d)
    integer(int64), intent(in) :: coefficients(0:)

    d = size(coefficients) - 1
    do while (d >= 0)
      if (coefficients(d) /= 0) return
      d = d - 1
    end do
  end function degree

  ! Increment a little-endian digit array, return false when it wraps to zero.
  logical function advance_digits(digits, radix) result(advanced)
    integer(int64), intent(inout) :: digits(0:)
    integer(int64), intent(in) :: radix
    integer :: j

    advanced = .true.
    do j = 0, size(digits) - 1
      if (digits(j) < radix - 1) then
        digits(j) = digits(j) + 1
        return
      end if
      digits(j) = 0
    end do
    advanced = .false.
  end function advance_digits

  integer(int64) function mod_power(base, exponent, modulus) result(value)
    integer(int64), intent(in) :: base, exponent, modulus
    integer(int64) :: factor, remaining

    value = 1
    factor = modulo(base, modulus)
    remaining = exponent
    do while (remaining > 0)
      if (btest(remaining, 0)) value = modulo(value * factor, modulus)
      remaining = remaining / 2
      factor = modulo(factor * factor, modulus)
    end do
  end function mod_power

  ! A monic divisor lets long division work even over a ring with zero divisors.
  subroutine monic_remainder(dividend, divisor, modulus, remainder)
    integer(int64), intent(in) :: dividend(0:), divisor(0:), modulus
    integer(int64), intent(out) :: remainder(0:)
    integer(int64) :: leading_coefficient
    integer :: power, j, divisor_degree

    remainder = 0
    remainder(0:size(dividend) - 1) = modulo(dividend, modulus)
    divisor_degree = degree(divisor)
    do power = degree(remainder), divisor_degree, -1
      leading_coefficient = remainder(power)
      do j = 0, divisor_degree
        remainder(power - divisor_degree + j) = modulo( &
          remainder(power - divisor_degree + j) - leading_coefficient * divisor(j), modulus)
      end do
    end do
  end subroutine monic_remainder

  subroutine polynomial_gcd(a, b, p, g)
    integer(int64), intent(in) :: a(0:), b(0:), p
    integer(int64), intent(out) :: g(0:)
    integer(int64) :: dividend(0:size(g) - 1), divisor(0:size(g) - 1), remainder(0:size(g) - 1), inverse
    integer :: leading_degree

    dividend = 0
    divisor = 0
    dividend(0:size(a) - 1) = a
    divisor(0:size(b) - 1) = b
    do
      leading_degree = degree(divisor)
      if (leading_degree < 0) exit
      inverse = mod_power(divisor(leading_degree), p - 2, p)
      divisor = modulo(divisor * inverse, p)
      call monic_remainder(dividend, divisor, p, remainder)
      dividend = divisor
      divisor = remainder
    end do
    leading_degree = degree(dividend)
    g = dividend
    if (leading_degree >= 0) g = modulo(dividend * mod_power(dividend(leading_degree), p - 2, p), p)
  end subroutine polynomial_gcd

  subroutine multiply_mod_polynomial(a, b, f, p, reduced_product)
    integer(int64), intent(in) :: a(0:), b(0:), f(0:), p
    integer(int64), intent(out) :: reduced_product(0:)
    integer(int64) :: product(0:2 * (size(f) - 1)), remainder(0:2 * (size(f) - 1))

    call polynomial_product(a, b, p, product)
    call monic_remainder(product, f, p, remainder)
    reduced_product = remainder(0:size(reduced_product) - 1)
  end subroutine multiply_mod_polynomial

  ! Raise a polynomial to the p-th power modulo the given divisor.
  subroutine frobenius(polynomial, divisor, p, power)
    integer(int64), intent(in) :: polynomial(0:), divisor(0:), p
    integer(int64), intent(out) :: power(0:)
    integer(int64) :: factor(0:size(polynomial) - 1), product(0:size(polynomial) - 1), remaining

    factor = polynomial
    power = 0
    power(0) = 1
    remaining = p
    do while (remaining > 0)
      if (btest(remaining, 0)) then
        call multiply_mod_polynomial(power, factor, divisor, p, product)
        power = product
      end if
      remaining = remaining / 2
      if (remaining > 0) then
        call multiply_mod_polynomial(factor, factor, divisor, p, product)
        factor = product
      end if
    end do
  end subroutine frobenius

  logical function is_field_irreducible(a, p) result(ok)
    integer(int64), intent(in) :: a(0:), p
    integer(int64) :: f(0:size(a) - 1), h(0:size(a) - 1), t(0:size(a) - 1), g(0:size(a) - 1)
    integer :: s, i

    f = modulo(a, p)
    s = degree(f)
    ok = .false.
    if (s < 1) return
    if (s == 1) then
      ok = .true.
      return
    end if
    f = modulo(f * mod_power(f(s), p - 2, p), p)
    h = 0
    h(1) = 1
    do i = 1, s / 2
      call frobenius(h(0:s), f(0:s), p, t(0:s))
      h(0:s) = t(0:s)
      t = h
      t(1) = modulo(t(1) - 1, p)
      call polynomial_gcd(f, t, p, g)
      if (degree(g) > 0) return
    end do
    ok = .true.
  end function is_field_irreducible

  ! Ordinary polynomial product without polynomial reduction. c must be large enough for a*b.
  subroutine polynomial_product(a, b, q, c)
    integer(int64), intent(in) :: a(0:), b(0:), q
    integer(int64), intent(out) :: c(0:)
    integer :: i, j

    c = 0
    do i = 0, degree(a)
      do j = 0, degree(b)
        c(i + j) = modulo(c(i + j) + a(i) * b(j), q)
      end do
    end do
  end subroutine polynomial_product

  ! Division by a monic polynomial, returning both quotient and remainder.
  subroutine divide_by_monic(dividend, divisor, modulus, quotient, remainder)
    integer(int64), intent(in) :: dividend(0:), divisor(0:), modulus
    integer(int64), intent(out) :: quotient(0:), remainder(0:)
    integer(int64) :: leading_coefficient
    integer :: power, j, divisor_degree

    quotient = 0
    remainder = 0
    remainder(0:size(dividend) - 1) = modulo(dividend, modulus)
    divisor_degree = degree(divisor)
    do power = degree(remainder), divisor_degree, -1
      leading_coefficient = remainder(power)
      quotient(power - divisor_degree) = leading_coefficient
      do j = 0, divisor_degree
        remainder(power - divisor_degree + j) = modulo( &
          remainder(power - divisor_degree + j) - leading_coefficient * divisor(j), modulus)
      end do
    end do
  end subroutine divide_by_monic

  ! Check whether the monic reduction is a power of one irreducible factor in F_p[x].
  ! In the first nontrivial GCD with x**(p**i)-x, all factors found
  ! have degree i and occur once. If degree(GCD) > i, there are
  ! different factors. Otherwise the GCD is the only candidate for the factor.
  ! If the search finds no such GCD through half the degree, the reduction is irreducible.
  subroutine primary_reduction(reduction, p, factor, exponent, is_primary)
    integer(int64), intent(in) :: reduction(0:), p
    integer(int64), intent(out) :: factor(0:)
    integer, intent(out) :: exponent
    logical, intent(out) :: is_primary
    integer(int64) :: x_power(0:size(reduction) - 1), difference(0:size(reduction) - 1)
    integer(int64) :: gcd(0:size(reduction) - 1), remaining(0:size(reduction) - 1)
    integer(int64) :: quotient(0:size(reduction) - 1), remainder(0:size(reduction) - 1)
    integer :: reduction_degree, factor_degree, gcd_degree

    is_primary = .false.
    exponent = 0
    factor = 0
    reduction_degree = degree(reduction)
    if (reduction_degree < 1) return
    x_power = 0
    x_power(1) = 1
    do factor_degree = 1, reduction_degree / 2
      call frobenius(x_power, reduction, p, difference)
      x_power = difference
      difference(1) = modulo(difference(1) - 1, p)
      call polynomial_gcd(reduction, difference, p, gcd)
      gcd_degree = degree(gcd)
      if (gcd_degree <= 0) cycle
      if (gcd_degree /= factor_degree) return
      remaining = reduction
      do while (degree(remaining) > 0)
        call divide_by_monic(remaining, gcd, p, quotient, remainder)
        if (any(remainder /= 0)) return
        exponent = exponent + 1
        remaining = quotient
      end do
      factor(0:size(gcd) - 1) = gcd
      is_primary = .true.
      return
    end do
    factor(0:size(reduction) - 1) = reduction
    exponent = 1
    is_primary = .true.
  end subroutine primary_reduction

  ! The input is monic and reduces modulo p to an irreducible factor power.
  ! In a proper factorization, the smaller monic divisor reduces to factor**j
  ! for some j <= exponent/2. Test all lifts of these powers modulo q.
  logical function is_monic_irreducible(polynomial, p, q, factor, exponent) result(irreducible)
    integer(int64), intent(in) :: polynomial(0:), p, q, factor(0:)
    integer, intent(in) :: exponent
    integer(int64) :: factor_power(0:size(polynomial) - 1), product(0:size(polynomial) - 1)
    integer(int64) :: candidate(0:size(polynomial) - 1)
    integer(int64) :: digits(0:size(polynomial) - 1), remainder(0:size(polynomial) - 1)
    integer :: j, candidate_degree, factor_degree

    irreducible = .true.
    if (exponent == 1) return
    if (q == p) then
      irreducible = .false.
      return
    end if
    factor_power = 0
    factor_power(0) = 1
    factor_degree = degree(factor)
    do j = 1, exponent / 2
      call polynomial_product(factor_power, factor, p, product)
      factor_power = product
      candidate_degree = j * factor_degree
      digits = 0
      do
        candidate = factor_power
        candidate(0:candidate_degree - 1) = factor_power(0:candidate_degree - 1) + p * digits(0:candidate_degree - 1)
        call monic_remainder(polynomial, candidate(0:candidate_degree), q, remainder)
        if (all(remainder == 0)) then
          irreducible = .false.
          return
        end if
        if (.not. advance_digits(digits(0:candidate_degree - 1), q / p)) exit
      end do
    end do
  end function is_monic_irreducible

  subroutine write_polynomial(f, line, file_unit)
    integer(int64), intent(in) :: f(0:)
    character(*), intent(inout) :: line
    character(64) :: term
    integer, intent(in) :: file_unit
    integer :: power, line_length, term_length, io_status

    line = ''
    line_length = 0
    do power = degree(f), 0, -1
      if (f(power) == 0) cycle
      if (line_length > 0) then
        line(line_length + 1:line_length + 3) = ' + '
        line_length = line_length + 3
      end if
      if (power == 0) then
        write (term, '(i0)') f(power)
      else if (power == 1) then
        write (term, '(i0,a)') f(power), '*x'
      else
        write (term, '(i0,a,i0)') f(power), '*x^', power
      end if
      term_length = len_trim(term)
      line(line_length + 1:line_length + term_length) = term(1:term_length)
      line_length = line_length + term_length
    end do
    write (file_unit, '(a)', iostat = io_status) line(1:line_length)
    if (io_status /= 0) error stop 'Failed to write the listing.'
  end subroutine write_polynomial
end module polynomials

program list_irreducibles
  use iso_fortran_env, only: error_unit, int64
  use polynomials, only: advance_digits, degree, is_monic_irreducible, mod_power, &
    polynomial_product, primary_reduction, write_polynomial
  implicit none
  ! This bound keeps products of modular residues within int64.
  integer(int64), parameter :: max_modulus = 2147483647_int64
  integer, parameter :: term_buffer_size = 64
  integer(int64) :: n, prime, prime_power, cofactor, total
  integer(int64) :: component_primes(32), component_moduli(32), crt_weights(32)
  integer(int64), allocatable :: monic_polynomial(:), reduction(:), factor(:), lift_digits(:)
  integer(int64), allocatable :: unit(:), associate_digits(:)
  integer(int64), allocatable :: choices(:, :), unit_digits(:, :), crt_sum(:, :), counts(:)
  character(:), allocatable :: line
  character(128) :: filename
  character(512) :: message
  integer :: max_degree, d, num_components, i, active_component, io_status, alloc_status, output_unit, exponent
  logical :: is_primary

  write (error_unit, '(a)') 'Enter n (2..2147483647) and the maximum degree (>=0):'
  read (*, *, iostat = io_status) n, max_degree
  if (io_status /= 0) error stop 'Invalid input.'
  if (n < 2 .or. n > max_modulus .or. max_degree < 0) error stop 'Input is out of range.'
  if (term_buffer_size * (int(max_degree, int64) + 1) > int(huge(max_degree), int64)) &
    error stop 'Maximum degree is too large for the output buffer.'
  call factor_modulus(n, component_primes, component_moduli, num_components)

  ! The Chinese remainder theorem (CRT) splits Z/nZ into prime-power rings.
  ! Each CRT weight is 1 modulo its own prime power and 0 modulo the others.
  do i = 1, num_components
    prime_power = component_moduli(i)
    cofactor = n / prime_power
    ! Euler's theorem gives the inverse of the cofactor modulo the prime power.
    crt_weights(i) = modulo(cofactor * &
      mod_power(cofactor, prime_power - prime_power / component_primes(i) - 1, prime_power), n)
  end do

  allocate(monic_polynomial(0:max_degree), reduction(0:max_degree), factor(0:max_degree), &
           lift_digits(0:max_degree), unit(0:max_degree), associate_digits(0:max_degree), &
           choices(0:max_degree, num_components), unit_digits(0:max_degree, num_components), &
           crt_sum(0:max_degree, 0:num_components), counts(0:max_degree), stat = alloc_status)
  if (alloc_status /= 0) error stop 'Not enough memory.'
  allocate(character(term_buffer_size * (max_degree + 1)) :: line, stat = alloc_status)
  if (alloc_status /= 0) error stop 'Not enough memory for output.'

  write (filename, '(a,i0,a,i0,a)') 'irreducibles_Z', n, '_degree', max_degree, '.txt'
  open (newunit = output_unit, file = trim(filename), status = 'replace', action = 'write', &
       form = 'formatted', iostat = io_status, iomsg = message)
  if (io_status /= 0) then
    write (error_unit, '(a)') trim(message)
    error stop 'Could not create the output file.'
  end if
  write (output_unit, '(a,i0,a,i0)', iostat = io_status) &
    'Irreducible polynomials over Z', n, ', maximum degree ', max_degree
  if (io_status /= 0) error stop 'Failed to write the header.'
  counts = 0
  total = 0
  crt_sum = 0

  ! Choose one irreducible component and fill the other components with units.
  do active_component = 1, num_components
    prime = component_primes(active_component)
    prime_power = component_moduli(active_component)
    write (error_unit, '(a,i0,a,i0)') 'Component ', active_component, ', modulus ', prime_power
    if (prime_power > prime) call generate_zero_divisors()
    do d = 1, max_degree
      ! Analyze each monic reduction in F_p[x] once.
      reduction = 0
      reduction(d) = 1
      do
        call primary_reduction(reduction(0:d), prime, factor(0:d), exponent, is_primary)
        ! Over the field, powers with exponent > 1 are reducible.
        if (is_primary .and. (prime_power > prime .or. exponent == 1)) then
          lift_digits = 0
          do
            monic_polynomial = 0
            monic_polynomial(0:d) = reduction(0:d) + prime * lift_digits(0:d)
            if (is_monic_irreducible(monic_polynomial(0:d), prime, prime_power, factor(0:d), exponent)) then
              call generate_associates(d)
            end if
            if (.not. advance_digits(lift_digits(0:d - 1), prime_power / prime)) exit
          end do
        end if
        if (.not. advance_digits(reduction(0:d - 1), prime)) exit
      end do
    end do
  end do

  do d = 0, max_degree
    write (*, '(a,i0,a,i0)') 'Degree ', d, ': ', counts(d)
    write (output_unit, '(a,i0,a,i0)', iostat = io_status) 'Degree ', d, ': ', counts(d)
    if (io_status /= 0) error stop 'Failed to write a degree count.'
  end do
  write (*, '(a,i0)') 'Total: ', total
  write (output_unit, '(a,i0)', iostat = io_status) 'Total: ', total
  if (io_status /= 0) error stop 'Failed to write the total.'
  write (output_unit, '(a)', iostat = io_status) 'Listing complete.'
  if (io_status /= 0) error stop 'Failed to finish the listing.'
  close (output_unit, iostat = io_status)
  if (io_status /= 0) error stop 'Failed to close the file.'
  write (*, '(a)') 'Listing saved to: ' // trim(filename)
contains

  subroutine factor_modulus(n, component_primes, component_moduli, num_components)
    integer(int64), intent(in) :: n
    integer(int64), intent(out) :: component_primes(:), component_moduli(:)
    integer, intent(out) :: num_components
    integer(int64) :: remaining, prime, prime_power

    ! Try 2, then odd divisors. Stop once the divisor exceeds sqrt(remaining).
    num_components = 0
    remaining = n
    prime = 2
    do while (prime <= remaining / prime)
      if (modulo(remaining, prime) == 0) then
        prime_power = 1
        do while (modulo(remaining, prime) == 0)
          prime_power = prime_power * prime
          remaining = remaining / prime
        end do
        num_components = num_components + 1
        component_primes(num_components) = prime
        component_moduli(num_components) = prime_power
      end if
      if (prime == 2) then
        prime = 3
      else
        prime = prime + 2
      end if
    end do
    ! Any factor left over is prime.
    if (remaining > 1) then
      num_components = num_components + 1
      component_primes(num_components) = remaining
      component_moduli(num_components) = remaining
    end if
  end subroutine factor_modulus

  subroutine generate_associates(monic_degree)
    integer, intent(in) :: monic_degree
    integer(int64) :: constant_term
    integer :: unit_degree

    ! Multiplying by units gives the associates of the monic polynomial.
    ! Their degrees add because the polynomial is monic, so cap the unit's degree.
    unit_degree = max_degree - monic_degree
    associate_digits = 0
    do
      unit = 0
      unit(0:unit_degree) = prime * associate_digits(0:unit_degree)
      do constant_term = 1, prime_power - 1
        if (modulo(constant_term, prime) == 0) cycle
        unit(0) = constant_term
        call polynomial_product(monic_polynomial(0:monic_degree), unit(0:unit_degree), &
          prime_power, choices(:, active_component))
        call combine_components(1)
      end do
      if (.not. advance_digits(associate_digits(1:unit_degree), prime_power / prime)) exit
    end do
  end subroutine generate_associates

  subroutine generate_zero_divisors()
    integer(int64) :: constant_term, unit_modulus, base

    ! An irreducible zero divisor is p times a unit.
    ! Multiplication by p identifies units that agree modulo q/p, so enumerate them once.
    unit_modulus = prime_power / prime
    base = unit_modulus / prime
    associate_digits = 0
    do
      choices(:, active_component) = (prime * prime) * associate_digits
      do constant_term = 1, unit_modulus - 1
        if (modulo(constant_term, prime) == 0) cycle
        choices(0, active_component) = prime * constant_term
        call combine_components(1)
      end do
      if (.not. advance_digits(associate_digits(1:max_degree), base)) exit
    end do
  end subroutine generate_zero_divisors

  recursive subroutine combine_components(component)
    integer, intent(in) :: component
    integer(int64) :: constant_term, base
    integer :: polynomial_degree

    if (component > num_components) then
      polynomial_degree = degree(crt_sum(:, num_components))
      if (polynomial_degree < 0) error stop 'Internal CRT error.'
      if (counts(polynomial_degree) == huge(total) .or. total == huge(total)) error stop 'Count exceeds the int64 limit.'
      counts(polynomial_degree) = counts(polynomial_degree) + 1
      total = total + 1
      call write_polynomial(crt_sum(:, num_components), line, output_unit)
      return
    end if
    ! The irreducible component already has a value.
    if (component == active_component) then
      crt_sum(:, component) = modulo(crt_sum(:, component - 1) + crt_weights(component) * choices(:, component), n)
      call combine_components(component + 1)
    else
      ! A unit's constant coefficient is coprime to p, its other coefficients are divisible by p.
      base = component_moduli(component) / component_primes(component)
      unit_digits(:, component) = 0
      do
        choices(:, component) = component_primes(component) * unit_digits(:, component)
        do constant_term = 1, component_moduli(component) - 1
          if (modulo(constant_term, component_primes(component)) == 0) cycle
          choices(0, component) = constant_term
          crt_sum(:, component) = modulo(crt_sum(:, component - 1) + crt_weights(component) * choices(:, component), n)
          call combine_components(component + 1)
        end do
        if (.not. advance_digits(unit_digits(1:max_degree, component), base)) exit
      end do
    end if
  end subroutine combine_components
end program list_irreducibles
