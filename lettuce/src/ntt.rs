use crate::*;

pub fn intt_negacyclic<const N: usize, E: FieldScalar>(
    input: impl Iterator<Item = E> + ExactSizeIterator,
) -> Result<impl Iterator<Item = E> + ExactSizeIterator> {
    assert!(input.len() == N);
    let outer_root_inv = match E::unity_root(2 * N) {
        Some(root) => root.inverse(),
        None => anyhow::bail!(
            "Z/{} does not have unity root cycle of len: {}",
            E::Q,
            input.len()
        ),
    };
    let n_inv = E::from(N).inverse();
    let input_vec = input.collect::<Vector<_>>();
    let o = (0..N).map(move |j| {
        (0..N)
            .map(|k| {
                let exp = (2 * k + 1) * j;
                let root = outer_root_inv.modpow(exp as u128);
                input_vec[k] * root
            })
            .sum::<E>()
            * n_inv
    });
    Ok(o)
}

pub fn ntt_negacyclic<const N: usize, E: FieldScalar>(
    input: impl Iterator<Item = E> + ExactSizeIterator,
) -> Result<impl Iterator<Item = E> + ExactSizeIterator> {
    assert!(input.len() == N);
    // TODO: unity root iterators, run in parallel with rayon
    // TODO: cache a root lookup table
    let outer_root = match E::unity_root(2 * N) {
        Some(root) => root,
        None => anyhow::bail!(
            "Z/{} does not have unity root cycle of len: {}",
            E::Q,
            input.len()
        ),
    };
    assert!(outer_root.modpow((2 * N) as u128) == E::one());
    assert!(outer_root.modpow(N as u128) == E::negone());
    let input_vec = input.collect::<Vector<_>>();
    let o = (0..N).map(move |k| {
        (0..N)
            .map(|j| {
                let exp = (2 * k + 1) * j;
                let root = outer_root.modpow(exp as u128);
                input_vec[j] * root
            })
            .sum::<E>()
    });
    Ok(o)
}

/// Number theoretic transform.
///
/// Given a vector of ring elements v and a root of unity u in a ring
/// cardinality q. Output a vector of evaluations at points in the root of unity.
pub fn ntt<const N: usize, E: FieldScalar>(
    input: impl Iterator<Item = E> + ExactSizeIterator + Clone,
) -> Result<impl Iterator<Item = E> + ExactSizeIterator + Clone> {
    assert!(input.len() == N);
    // TODO: unity root iterators, run in parallel with rayon
    // TODO: cache a root lookup table
    let root = match E::unity_root(N) {
        Some(root) => root,
        None => anyhow::bail!(
            "Z/{} does not have unity root cycle of len: {}",
            E::Q,
            input.len()
        ),
    };
    debug_assert!(
        root.modpow(input.len() as u128) == E::one(),
        "lettuce::ntt root of unity incorrect cycle length"
    );
    Ok((0..N).map(move |k| {
        input
            .clone()
            .enumerate()
            .map(|(j, v)| v * root.modpow((k * j) as u128))
            .sum::<E>()
    }))
}

pub fn intt<const N: usize, E: FieldScalar>(
    input: impl Iterator<Item = E> + ExactSizeIterator + Clone,
) -> Result<impl Iterator<Item = E> + ExactSizeIterator + Clone> {
    assert!(input.len() == N);
    // const N_INV =
    let root = match E::unity_root(N) {
        Some(root) => root,
        None => anyhow::bail!(
            "Z/{} does not have unity root cycle of len: {}",
            E::Q,
            input.len()
        ),
    };
    debug_assert!(
        root.modpow(input.len() as u128) == E::one(),
        "lettuce::intt root of unity incorrect cycle length"
    );
    let input = input.collect::<Vector<_>>();
    Ok((0..N).map(move |j| {
        E::from(N).inverse()
            * (0..N)
                .map(|k| input[k] * root.modpow((j * k) as u128).inverse())
                .sum::<E>()
    }))
}

#[test]
fn ntt_field() -> Result<()> {
    let rng = &mut rand::rng();
    type E = MilliScalar;
    const N: usize = 64;
    for _ in 0..100 {
        let coefs = Vector::<E>::sample_uniform(N, rng);
        let root = E::unity_root(N).unwrap();
        for (i, eval) in ntt::<N, _>(coefs.iter().copied()).unwrap().enumerate() {
            // eval should be coefs evaluated at x = unity root
            let out = Polynomial::<N, _>::from(&coefs).evaluate(root.modpow(i as u128));
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
        let a = Polynomial::<N, E>::sample_uniform(rng);
        let b = Polynomial::<N, E>::sample_uniform(rng);
        let c = a * b;

        let a_ntt = ntt_negacyclic::<N, _>(a.coefs())?;
        let b_ntt = ntt_negacyclic::<N, _>(b.coefs())?;
        let c_ntt = a_ntt.zip(b_ntt).map(|(a_v, b_v)| a_v * b_v);
        let c_computed =
            Polynomial::<N, E>::from(&intt_negacyclic::<N, _>(c_ntt)?.collect::<Vector<_>>());
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
        let coefs = Vector::<E>::sample_uniform(N, rng);
        let evals =
            (0..N).map(|i| Polynomial::<N, E>::from(&coefs).evaluate(root.modpow(i as u128)));
        let out = intt::<N, _>(evals).unwrap().collect::<Vector<_>>();
        assert_eq!(out, coefs);
    }
    Ok(())
}
