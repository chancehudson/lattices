use crate::*;

/// Negacyclic number theoretic transform. Applies a number theoretic transform using a
/// 2N'th root of unity to yield the negacyclic convolution of two vectors. This is equivalent to
/// polynomial multiplication modulo X^N + 1
///
/// Panics if the not root of unity exists with order `2N`.
pub fn ntt_negacyclic<const N: usize, E: FieldScalar>(input: &mut [E; N]) {
    // TODO: unity root iterators, run in parallel with rayon
    // TODO: cache a root lookup table
    let psi = match E::unity_root(2 * N) {
        Some(root) => root,
        None => panic!(
            "Z/{} does not have unity root cycle of len: {}",
            E::Q,
            input.len()
        ),
    };
    debug_assert!(psi.modpow((2 * N) as u128) == E::one());
    debug_assert!(psi.modpow(N as u128) == E::negone());

    let psi_data = E::unity_root_powers(psi, N);
    let (_psi_inv, psi_powers, _psi_inv_powers) = psi_data.as_ref();

    for (j, x) in input.iter_mut().enumerate() {
        *x *= psi_powers[j];
    }

    let omega = psi * psi;
    let omega_data = E::unity_root_powers(omega, N);
    let (_omega_inv, omega_powers, _omega_inv_powers) = omega_data.as_ref();

    ntt_inplace::<N, E>(input, omega_powers);
}

/// Inverse negacyclic number theoretic transform. Applies an inverse number theoretic transform using a
/// 2N'th root of unity to convert the evaluation form of a polynomial modulo X^N + 1 back to
/// coefficient form.
///
/// Panics if the not root of unity exists with order `2N`.
pub fn intt_negacyclic<const N: usize, E: FieldScalar>(input: &mut [E; N]) {
    let psi = match E::unity_root(2 * N) {
        Some(root) => root,
        None => panic!(
            "Z/{} does not have unity root cycle of len: {}",
            E::Q,
            input.len()
        ),
    };
    let psi_data = E::unity_root_powers(psi, N);
    let (_psi_inv, _psi_powers, psi_inv_powers) = psi_data.as_ref();
    let omega = psi * psi;
    let omega_data = E::unity_root_powers(omega, N);
    let (_omega_inv, _omega_powers, omega_inv_powers) = omega_data.as_ref();

    ntt_inplace::<N, E>(input, omega_inv_powers);

    // post-multiply by psi^(-j) / N
    let n_inv = E::from(N).inverse();
    for (j, x) in input.iter_mut().enumerate() {
        *x *= psi_inv_powers[j] * n_inv;
    }
}

/// Bit-reversal permutation
fn bit_reverse_copy<E: Copy>(a: &mut [E]) {
    let n = a.len();
    let log_n = n.trailing_zeros() as usize;

    for i in 0..n {
        let j = reverse_bits(i, log_n);
        if i < j {
            a.swap(i, j);
        }
    }
}

/// Reverse the bottom `bits` bits of `n`
fn reverse_bits(mut n: usize, bits: usize) -> usize {
    let mut result = 0;
    for _ in 0..bits {
        result = (result << 1) | (n & 1);
        n >>= 1;
    }
    result
}

fn ntt_inplace<const N: usize, E: FieldScalar>(input: &mut [E], powers: &Vec<E>) {
    bit_reverse_copy(input);
    let mut len = 2;
    while len <= N {
        // let w_len = root.modpow((N / len) as u128);
        let w_len = powers[N / len];

        for i in (0..N).step_by(len) {
            let mut w = E::one();
            for j in 0..len / 2 {
                let u = input[i + j];
                let v = input[i + j + len / 2] * w;
                input[i + j] = u + v;
                input[i + j + len / 2] = u - v;
                w = w * w_len;
            }
        }
        len *= 2;
    }
}

/// Number theoretic transform.
///
/// Given a vector of ring elements v and a root of unity u in a ring of
/// cardinality q. Output a vector of evaluations at points in the root of unity.
///
/// Panics if the not root of unity exists with order `N`.
pub fn ntt<const N: usize, E: FieldScalar>(input: &mut [E; N]) {
    // TODO: unity root iterators, run in parallel with rayon
    // TODO: cache a root lookup table
    let root = match E::unity_root(N) {
        Some(root) => root,
        None => panic!(
            "Z/{} does not have unity root cycle of len: {}",
            E::Q,
            input.len()
        ),
    };
    debug_assert!(
        root.modpow(input.len() as u128) == E::one(),
        "lettuce::ntt root of unity incorrect cycle length"
    );
    let root_data = E::unity_root_powers(root, N);
    let (_root_inv, root_powers, _root_inv_powers) = root_data.as_ref();
    ntt_inplace::<N, E>(input, root_powers);
}

/// Inverse number theoretic transform.
///
/// Given a vector of ring elements in evaluation form convert to a vector of elements in
/// coefficient form.
///
/// Panics if the not root of unity exists with order `N`.
pub fn intt<const N: usize, E: FieldScalar>(input: &mut [E; N]) {
    let root = match E::unity_root(N) {
        Some(root) => root,
        None => panic!(
            "Z/{} does not have unity root cycle of len: {}",
            E::Q,
            input.len()
        ),
    };
    debug_assert!(
        root.modpow(input.len() as u128) == E::one(),
        "lettuce::intt root of unity incorrect cycle length"
    );

    let root_data = E::unity_root_powers(root, N);
    let (_root_inv, _root_powers, root_inv_powers) = root_data.as_ref();

    ntt_inplace::<N, E>(input, root_inv_powers);

    let n_inv = E::from(N).inverse();
    for x in input {
        *x *= n_inv;
    }
}

#[test]
fn ntt_negacyclic_roundtrip() {
    type E = MilliScalar;
    const N: usize = 8;

    let rng = &mut rand::rng();

    for _ in 0..10 {
        let mut orig = Polynomial::<N, E>::sample_uniform(rng);
        let orig_clone = orig.clone();

        ntt_negacyclic::<N, _>(orig.coefs_slice_mut());
        intt_negacyclic::<N, _>(orig.coefs_slice_mut());
        assert_eq!(orig, orig_clone);
    }
}

#[test]
fn ntt_roundtrip() -> Result<()> {
    type E = MilliScalar;
    const N: usize = 8;

    let rng = &mut rand::rng();

    for _ in 0..10 {
        let mut orig = Polynomial::<N, E>::sample_uniform(rng);
        let orig_clone = orig.clone();

        ntt::<N, _>(orig.coefs_slice_mut());
        intt::<N, _>(orig.coefs_slice_mut());
        assert_eq!(orig, orig_clone);
    }

    Ok(())
}

#[test]
fn ntt_field() -> Result<()> {
    let rng = &mut rand::rng();
    type E = MilliScalar;
    const N: usize = 64;
    for _ in 0..100 {
        let mut poly = Polynomial::<N, _>::from(&Vector::<E>::sample_uniform(N, rng));
        let poly_clone = poly.clone();
        let root = E::unity_root(N).unwrap();
        ntt::<N, _>(poly.coefs_slice_mut());

        for (i, eval) in poly.coefs().enumerate() {
            // eval should be coefs evaluated at x = unity root
            let out = poly_clone.evaluate(root.modpow(i as u128));
            assert_eq!(eval, out);
        }
    }
    Ok(())
}

#[test]
fn ntt_negacyclic_mul() -> Result<()> {
    let rng = &mut rand::rng();
    type E = MilliScalar;
    const N: usize = 64;
    for _ in 0..100 {
        let mut a = Polynomial::<N, E>::sample_uniform(rng);
        let mut b = Polynomial::<N, E>::sample_uniform(rng);
        let c = a * b;

        ntt_negacyclic::<N, _>(a.coefs_slice_mut());
        ntt_negacyclic::<N, _>(b.coefs_slice_mut());
        let c_ntt = a
            .coefs()
            .zip(b.coefs())
            .map(|(a_v, b_v)| a_v * b_v)
            .collect::<Vector<_>>();
        let mut c_computed = Polynomial::<N, _>::from(&c_ntt);
        intt_negacyclic::<N, _>(c_computed.coefs_slice_mut());
        assert_eq!(c, c_computed);
    }
    Ok(())
}

#[test]
fn intt_field() -> Result<()> {
    let rng = &mut rand::rng();
    type E = MilliScalar;
    const N: usize = 64;
    let root = E::unity_root(N).unwrap();
    for _ in 0..100 {
        let poly = Polynomial::<N, _>::from(&Vector::<E>::sample_uniform(N, rng));
        let evals = (0..N)
            .map(|i| poly.evaluate(root.modpow(i as u128)))
            .collect::<Vector<_>>();
        let mut evals_poly = Polynomial::from(&evals);
        intt::<N, _>(evals_poly.coefs_slice_mut());
        assert_eq!(evals_poly, poly);
    }
    Ok(())
}
