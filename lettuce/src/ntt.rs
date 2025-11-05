use crate::*;

/// Number theoretic transform.
///
/// Given a vector of ring elements v and a root of unity u in a ring
/// cardinality q. Output a vector of evaluations at points in the root of unity.
pub fn ntt<const N: usize, E: FieldScalar>(input: &Vector<E>) -> Result<impl Iterator<Item = E>> {
    assert!(input.len() <= N);
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
        (0..N)
            .map(|j| input[j] * root.modpow((k * j) as u128))
            .sum()
    }))
}

pub fn intt<const N: usize, E: FieldScalar>(input: &Vector<E>) -> Result<impl Iterator<Item = E>> {
    assert!(input.len() <= N);
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
        for (i, eval) in ntt::<N, _>(&coefs).unwrap().enumerate() {
            // eval should be coefs evaluated at x = unity root
            let out = Polynomial::<N, _>::from(&coefs).evaluate(root.modpow(i as u128));
            assert_eq!(eval, out);
        }
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
        let evals = (0..N)
            .map(|i| Polynomial::<N, E>::from(&coefs).evaluate(root.modpow(i as u128)))
            .collect::<Vector<_>>();
        let out = intt::<N, _>(&evals).unwrap().collect::<Vector<_>>();
        assert_eq!(out, coefs);
    }
    Ok(())
}
