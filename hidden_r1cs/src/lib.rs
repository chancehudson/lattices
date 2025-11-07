use anyhow::Result;
use lettuce::*;
use rand::Rng;
use rand::SeedableRng;

/// Given an r1cs instance and a witness generate a hidden witness that fulfills a relaxed R1CS
/// relation and lattice NIZK arguments of consistency.
pub struct HiddenR1CS<const N: usize, E: FieldScalar> {
    /// The base r1cs constraint system
    r1cs: R1CS<E>,
    /// The final s value in the RR1CS
    wtns_s_final: E,
    /// The final masked and challenged witness that fulfills the RR1CS
    wtns_final: Vector<E>,
    /// The e vector in the final RR1CS relation
    wtns_e_final: Vector<E>,
    /// Argument that witness vector was multiplied by challenge
    wtns_arg: BDLOPLinearNIZKArg<N, E>,
    /// Commitment to witness mask
    wtns_mask_commit: BDLOP<N, E>,
    /// E vector in the relaxed r1cs for wtns_mask
    wtns_mask_e_commit: BDLOP<N, E>,
    /// Linear zk arg of crossterm challenge
    crossterm_arg: BDLOPLinearNIZKArg<N, E>,
}

impl<const N: usize, E: FieldScalar> HiddenR1CS<N, E> {
    /// Commit to a witness and generate a hidden witness that fulfills a relaxed R1CS instance
    /// formed from the input r1cs.
    ///
    /// The input witness is a vector of scalars. This vector will be encoded to a vector of
    /// polynomials for operation on the lattice. During verification we'll decode the hidden
    /// witness from a vector of polynomials to a vector of scalars.
    pub fn commit<R: Rng>(wtns: Vector<E>, r1cs: R1CS<E>, rng: &mut R) -> Result<Self> {
        assert!(
            r1cs.eval(&wtns)?.is_zero(),
            "hidden_r1cs::commit provided r1cs is not solved"
        );
        let wtns_mask = Vector::<E>::sample_uniform(wtns.len(), rng);
        let wtns_mask_e = r1cs.eval(&wtns_mask)?;
        // build the cross term
        let crossterm = (&r1cs.a * &wtns_mask) * &(&r1cs.b * &wtns)
            + &((&r1cs.a * &wtns) * &(&r1cs.b * &wtns_mask))
            - (&r1cs.c * &wtns)
            - (&r1cs.c * &wtns_mask);

        // check that our cross terms match the combined eval
        assert!(
            crossterm
                == &r1cs.eval(&(wtns.clone() + &wtns_mask))?
                    - (&r1cs.a * &wtns) * &(&r1cs.b * &wtns)
                    - (&r1cs.a * &wtns_mask) * &(&r1cs.b * &wtns_mask)
        );

        let wtns_polys = wtns
            .chunks(N)
            .map(|v| v.iter().copied().collect::<Vector<E>>())
            .map(|v| Polynomial::<N, _>::from(&v))
            .collect::<Vector<_>>();
        let wtns_mask_polys = wtns_mask
            .chunks(N)
            .map(|v| v.iter().copied().collect::<Vector<E>>())
            .map(|v| Polynomial::<N, _>::from(&v))
            .collect::<Vector<_>>();
        let crossterm_polys = crossterm
            .chunks(N)
            .map(|v| v.iter().copied().collect::<Vector<E>>())
            .map(|v| Polynomial::<N, _>::from(&v))
            .collect::<Vector<_>>();

        let wtns_mask_e_polys = wtns_mask_e
            .chunks(N)
            .map(|v| v.iter().copied().collect::<Vector<E>>())
            .map(|v| Polynomial::<N, _>::from(&v))
            .collect::<Vector<_>>();

        let wtns_lattice = BDLOP::lattice_for(wtns_mask_polys.len(), rng);
        let crossterm_lattice = BDLOP::lattice_for(crossterm_polys.len(), rng);

        let (wtns_mask_commit, wtns_mask_secret) =
            BDLOP::commit(wtns_mask_polys.clone(), &wtns_lattice, rng);
        let (wtns_mask_e_commit, wtns_mask_e_secret) =
            BDLOP::commit(wtns_mask_e_polys.clone(), &crossterm_lattice, rng);

        let (wtns_commit, wtns_secret) = BDLOP::commit(wtns_polys.clone(), &wtns_lattice, rng);
        let (crossterm_commit, crossterm_secret) =
            BDLOP::commit(crossterm_polys.clone(), &crossterm_lattice, rng);

        // need to hash wtns_commit, wtns_mask_commit, wtns_mask_e_commit, crossterm_commit
        //
        // sample a challenge
        let hash = blake3::hash(
            &wtns_commit
                .as_le_bytes()
                .chain(wtns_mask_commit.as_le_bytes())
                .chain(wtns_mask_e_commit.as_le_bytes())
                .chain(crossterm_commit.as_le_bytes())
                .collect::<Vec<_>>(),
        );
        let mut csprng = rand_chacha::ChaCha20Rng::from_seed(hash.into());
        let challenge = E::sample_uniform(&mut csprng);

        // we have our challenge, generate zk arguments of linear relations
        let wtns_challenged = wtns_polys
            .iter()
            .map(|w| *w * challenge)
            .collect::<Vector<_>>();
        let crossterm_challenged = crossterm_polys
            .iter()
            .map(|v| *v * challenge)
            .collect::<Vector<_>>();

        let wtns_final = wtns_polys
            .iter()
            .zip(wtns_mask_polys.iter())
            .map(|(w, w_m)| *w_m + *w * challenge)
            .collect::<Vector<_>>();
        let crossterm_final = crossterm_polys
            .iter()
            .map(|v| *v * challenge)
            .collect::<Vector<_>>()
            + &wtns_mask_e_polys;

        // generate scalar values that will be stored on the struct
        let wtns_final_scalar = wtns_final
            .iter()
            .flat_map(|poly| poly.coefs().collect::<Vector<_>>())
            .take(r1cs.width())
            .collect::<Vector<E>>();
        let crossterm_final_scalar = crossterm_final
            .iter()
            .flat_map(|poly| poly.coefs().collect::<Vector<_>>())
            .take(r1cs.height())
            .collect::<Vector<E>>();

        let (wtns_final_commit, wtns_final_secret) =
            BDLOP::commit(wtns_challenged, &wtns_lattice, rng);
        let (crossterm_final_commit, crossterm_final_secret) =
            BDLOP::commit(crossterm_challenged, &crossterm_lattice, rng);

        let wtns_arg = BDLOP::try_open_linear_zk(
            (wtns_commit, &wtns_secret),
            (wtns_final_commit, &wtns_final_secret),
            Polynomial::from(challenge),
            rng,
        )?;

        let crossterm_arg = BDLOP::try_open_linear_zk(
            (crossterm_commit, &crossterm_secret),
            (crossterm_final_commit, &crossterm_final_secret),
            Polynomial::from(challenge),
            rng,
        )?;

        debug_assert!(
            crossterm_final_scalar
                == (&r1cs.a * &wtns_final_scalar) * &(&r1cs.b * &wtns_final_scalar)
                    - (&r1cs.c * &wtns_final_scalar) * (E::one() + challenge),
            "hidden_r1cs final relation mismatch"
        );

        Ok(Self {
            wtns_final: wtns_final_scalar,
            wtns_e_final: crossterm_final_scalar,
            wtns_s_final: E::one() + challenge,
            wtns_mask_e_commit,
            wtns_arg,
            crossterm_arg,
            r1cs,
            wtns_mask_commit,
        })
    }

    pub fn verify(&self) -> Result<()> {
        let hash = blake3::hash(
            &self
                .wtns_arg
                .c_a
                .as_le_bytes()
                .chain(self.wtns_mask_commit.as_le_bytes())
                .chain(self.wtns_mask_e_commit.as_le_bytes())
                .chain(self.crossterm_arg.c_a.as_le_bytes())
                .collect::<Vec<_>>(),
        );
        let mut csprng = rand_chacha::ChaCha20Rng::from_seed(hash.into());
        let challenge = E::sample_uniform(&mut csprng);

        self.wtns_arg.verify()?;
        self.crossterm_arg.verify()?;

        let wtns_final_commit = self.wtns_mask_commit.clone() + &self.wtns_arg.c_b;
        let wtns_e_final_commit = self.wtns_mask_e_commit.clone() + &self.crossterm_arg.c_b;
        // TODO: open these

        if self.wtns_e_final
            != (&self.r1cs.a * &self.wtns_final) * &(&self.r1cs.b * &self.wtns_final)
                - (&self.r1cs.c * &self.wtns_final) * (E::one() + challenge)
        {
            anyhow::bail!("hidden_r1cs: verification failed constraint matrices are not fulfilled");
        }

        Ok(())
    }
}

#[test]
fn commit_rand() -> Result<()> {
    let rng = &mut rand::rng();
    type E = MilliScalar;
    const N: usize = 1024;
    let (r1cs, wtns) = R1CS::<E>::sample_uniform(120, 192, rng);
    let arg = HiddenR1CS::<64, _>::commit(wtns, r1cs, rng)?;
    arg.verify()?;

    Ok(())
}
