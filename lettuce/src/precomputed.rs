/// Stuff that should probably be pre-computed in production systems.
use std::collections::HashMap;

/// Greatest common denominator of two numbers.
fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        let tmp = b;
        b = a % b;
        a = tmp;
    }
    a
}

/// Return a prime factor of the input, if one exists.
/// Will overflow the stack if the input is prime.
pub fn pollard_rho(v: u128) -> u128 {
    pollard_rho_inner(v, 2, 1)
}

fn pollard_rho_inner(v: u128, x0: u128, c: u128) -> u128 {
    if v % 2 == 0 {
        return 2;
    }
    let (mut x, mut y, mut d) = (x0, x0, 1);
    let f = |x: u128| (x * x + c) % v;
    while d == 1 {
        x = f(x);
        y = f(f(y));
        d = gcd(x.abs_diff(y), v);
    }
    if d != v {
        d
    } else {
        pollard_rho_inner(v, (x0 + 1) % v, (x + 1) % v)
    }
}

#[test]
fn pollards_rho_factors() {
    let v_orig = 100u64;
    let mut prod = 1u64;
    let mut v = v_orig;

    while !const_primes::is_prime(v) {
        let prime_maybe = pollard_rho(v as u128) as u64;
        assert!(const_primes::is_prime(prime_maybe));
        prod *= prime_maybe;
        assert_eq!(v % prime_maybe, 0);
        v /= prime_maybe;
    }
    assert_eq!(v * prod, v_orig);
}

/// Compute the prime factorization of a number.
pub fn prime_factorize(mut v: u128) -> HashMap<u128, usize> {
    assert!(v > 1, "cannot factorize the identity values");
    let mut factors = HashMap::default();
    for p in const_primes::primes::<1024>() {
        if v == 1 {
            break;
        }
        let p = p as u128;
        let mut count = 0usize;
        // prime p might be a factor of v multiple times
        // determine the number
        loop {
            let quotient = v / p;
            if quotient * p != v {
                break;
            }
            v = quotient;
            count += 1;
        }
        if count > 0 {
            factors.insert(p, count);
        }
    }
    // we've cycled the first few primes and haven't returned.
    // use pollards rho to find the remaining prime factors
    while v > 1 {
        if const_primes::is_prime(v as u64) {
            *factors.entry(v).or_default() += 1;
            break;
        }
        let factor = pollard_rho(v);
        let quotient = v / factor;
        assert_eq!(factor * quotient, v, "inexact divisor");
        v = quotient;
        *factors.entry(v).or_default() += 1;
    }
    factors
}

#[test]
fn prime_factorize_test() {
    for v in 2..100_000 {
        let mut prod = 1u128;
        for (factor, count) in prime_factorize(v) {
            assert!(const_primes::is_prime(factor.try_into().unwrap()));
            for _ in 0..count {
                prod *= factor;
            }
        }
        assert_eq!(prod, v);
    }
}

pub fn find_unity_root(len: usize, q: u128, generator: u128) -> Option<u128> {
    let len = len as u128;
    let div = (q - 1) / len;
    if div * len != q - 1 {
        return None;
    }
    let omega = modpow(generator, div, q);

    #[cfg(debug_assertions)]
    {
        assert_eq!(modpow(omega, len as u128, q), 1);
    }

    Some(omega)
}

pub fn find_generator(q: u128, factorization: &HashMap<u128, usize>) -> u128 {
    let mut g = 2u128;
    loop {
        assert!(g < q);
        if modpow(g, q - 1, q) != 1 {
            g += 1;
            continue;
        }
        // check that the order of g is q
        let mut is_found = true;
        for (factor, _) in factorization {
            if modpow(g, (q - 1) / factor, q) == 1 {
                is_found = false;
                break;
            }
        }
        if is_found {
            return g;
        }
        g += 1;
    }
}

fn modpow(mut v: u128, mut exp: u128, q: u128) -> u128 {
    let mut out = 1u128;
    loop {
        if exp % 2 == 1 {
            out *= v;
            out %= q;
        }
        exp >>= 1;
        if exp == 0 {
            break;
        }
        v = (v * v) % q;
    }
    out
}
