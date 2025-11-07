use std::sync::Arc;

use anyhow::Result;
use rand::SeedableRng;

use crate::*;

const SIGMA: f64 = 27000f64;

/// Argues knowledge of two BDLOP commitments such that
/// c_a = g * c_b
pub struct BDLOPLinearNIZKArg<const N: usize, E: FieldScalar> {
    pub g: Polynomial<N, E>,
    pub c_a: BDLOP<N, E>,
    pub c_b: BDLOP<N, E>,
    pub t_1: Vector<Polynomial<N, E>>,
    pub t_2: Vector<Polynomial<N, E>>,
    pub u: Vector<Polynomial<N, E>>,
    pub z_1: Vector<Polynomial<N, E>>,
    pub z_2: Vector<Polynomial<N, E>>,
}

impl<const N: usize, E: FieldScalar> BDLOPLinearNIZKArg<N, E> {
    pub fn verify(&self) -> Result<()> {
        // check z terms norm bounds
        let l2_norm_bound = 4.0 * SIGMA * (N as f64).sqrt();
        for (z_1, z_2) in self.z_1.iter().zip(self.z_2.iter()) {
            if z_1.norm_l2() >= l2_norm_bound {
                anyhow::bail!("BDLOPLinearNIZKArg z_1 is beyond bounds");
            }
            if z_2.norm_l2() >= l2_norm_bound {
                anyhow::bail!("BDLOPLinearNIZKArg z_2 is beyond bounds");
            }
        }
        // compute the d value
        let hash = blake3::hash(
            &self
                .t_1
                .iter()
                .chain(self.t_2.iter())
                .chain(self.u.iter())
                .flat_map(|p| p.as_le_bytes())
                .collect::<Vec<_>>(),
        );
        let mut csprng = rand_chacha::ChaCha20Rng::from_seed(hash.into());

        let d: Polynomial<N, E> =
            std::array::from_fn(|_| E::at_displacement(csprng.random_range(0..2) - 1)).into();

        // check equalities
        if self.c_a.a_1.as_ref() * &self.z_1 != self.t_1.clone() + &(self.c_a.c_1.clone() * d) {
            anyhow::bail!("BDLOPLinearNIZKArg z_1 equality mismatch");
        }
        if self.c_a.a_1.as_ref() * &self.z_2 != self.t_2.clone() + &(self.c_b.c_1.clone() * d) {
            anyhow::bail!("BDLOPLinearNIZKArg z_2 equality mismatch");
        }

        let lhs = self.c_a.a_2.as_ref() * &self.z_1 * self.g - self.c_a.a_2.as_ref() * &self.z_2;
        let rhs = (self.c_a.c_2.clone() * self.g - &self.c_b.c_2) * d + &self.u;
        if lhs != rhs {
            anyhow::bail!("BDLOPLinearNIZKArg final equality mismatch");
        }

        Ok(())
    }
}

/// An implementation of Baum et. al. commitments.
/// https://eprint.iacr.org/2016/997.pdf
///
#[derive(Clone, Debug)]
pub struct BDLOP<const N: usize, E: FieldScalar> {
    a_1: Arc<Matrix<Polynomial<N, E>>>,
    a_2: Arc<Matrix<Polynomial<N, E>>>,
    c_1: Vector<Polynomial<N, E>>,
    c_2: Vector<Polynomial<N, E>>,
}

impl<const N: usize, E: FieldScalar> BDLOP<N, E> {
    /// Return a byte representation of the commitments.
    pub fn as_le_bytes(&self) -> impl Iterator<Item = u8> {
        self.c_1
            .iter()
            .flat_map(|v| v.as_le_bytes())
            .chain(self.c_2.iter().flat_map(|v| v.as_le_bytes()))
    }

    /// Check the parameters for safe operation.
    fn check_params() -> Result<()> {
        // our N value must be a power of 2
        if 2usize.pow(N.ilog2()) != N {
            anyhow::bail!("BDLOP degree of polynomial ring must be a power of 2");
        }
        // we don't check the field for invertibility. Instead we'll use a computational statement
        // about invertibility during proving/verification.
        Ok(())
    }

    /// Given a message length determine the dimension of a commitment matrix.
    ///
    /// BDLOP commitments vertically compose an SIS commitment to 0 with an SIS commitment
    /// to a message vector. We need two lattice bases, A_1 and A_2.
    ///
    /// Pages 11 and 10 of https://eprint.iacr.org/2016/997.pdf
    fn dimension(msg_len: usize) -> (usize, usize) {
        // approx n*log(q) for the message length only
        // we want our message to be mixed with at least this many vectors of random elements
        //
        // we subtract msg_len because the A_2 matrix will provide an additional msg_len
        // vectors of mixing
        let a_1_height = msg_len;
        // the A_2 matrix always has height equal to the message length
        let width = a_1_height + msg_len;
        (a_1_height, width)
    }

    pub fn lattice_for<R: Rng>(
        msg_len: usize,
        rng: &mut R,
    ) -> (Arc<Matrix<Polynomial<N, E>>>, Arc<Matrix<Polynomial<N, E>>>) {
        let (a_1_height, width) = Self::dimension(msg_len);
        // the A_1 lattice base
        let a_1 = Matrix::identity(a_1_height).compose_horizontal(Matrix::random(
            a_1_height,
            width - a_1_height,
            rng,
        ));

        // the A_2 lattice base
        let a_2 = Matrix::zero(msg_len, a_1_height)
            .compose_horizontal(Matrix::identity(msg_len))
            .compose_horizontal(Matrix::random(msg_len, width - a_1_height - msg_len, rng));
        (Arc::new(a_1), Arc::new(a_2))
    }

    /// Generate a BDLOP commitment to a vector of scalar elements.
    ///
    /// Elements will be encoded into polynomials.
    pub fn commit<R: Rng>(
        val: Vector<Polynomial<N, E>>,
        lattice: &(Arc<Matrix<Polynomial<N, E>>>, Arc<Matrix<Polynomial<N, E>>>),
        rng: &mut R,
    ) -> (Self, Vector<Polynomial<N, E>>) {
        let (a_1, a_2) = lattice.clone();

        // the secret committing to the zero component
        let r = (0..a_1.width())
            .map(|_| {
                let mut p = Polynomial::default();
                for coef in p.coefs_mut() {
                    let sample: i32 = rng.random_range(0..=2) - 1;
                    *coef = E::at_displacement(sample);
                }
                p
            })
            .collect::<Vector<_>>();

        let c_1 = a_1.as_ref() * &r;
        let c_2 = a_2.as_ref() * &r + &val;

        (Self { a_1, a_2, c_1, c_2 }, r)
    }

    /// Attempt to open a commitment directly using the stored `r` value.
    /// First attempts to open c_1 to the zero vector. If this succeeds c_2 is opened to whatever
    /// value is committed.
    pub fn try_open(&self, r: &Vector<Polynomial<N, E>>) -> Result<Vector<Polynomial<N, E>>> {
        if self.a_1.as_ref() * r != self.c_1 {
            anyhow::bail!("Failed to open commitment, secret is incorrect");
        }
        Ok(&self.c_2 - self.a_2.as_ref() * &r)
    }

    /// Attempt to generate a non-interactive ZK argument of opening.
    ///
    /// Described on page 15 of https://eprint.iacr.org/2016/997.pdf
    pub fn try_open_zk<R: Rng>(
        &self,
        r: &Vector<Polynomial<N, E>>,
        rng: &mut R,
    ) -> Result<(
        Polynomial<N, E>,
        Vector<Polynomial<N, E>>,
        Vector<Polynomial<N, E>>,
    )> {
        if self.a_1.as_ref() * r != self.c_1 {
            anyhow::bail!("Failed to open commitment, secret is incorrect");
        }
        loop {
            let y = (0..self.a_1.width())
                .map(|_| Polynomial::sample_gaussian(SIGMA, rng))
                .collect::<Vector<_>>();
            let t = self.a_1.as_ref() * &y;
            let hash = blake3::hash(&t.iter().flat_map(|v| v.as_le_bytes()).collect::<Vec<u8>>());
            let mut csprng = rand_chacha::ChaCha20Rng::from_seed(hash.into());

            let d: Polynomial<N, E> =
                std::array::from_fn(|_| E::at_displacement(csprng.random_range(0..2) - 1)).into();
            let z = y.clone() + &(r.clone() * d);
            // TODO: rejection sampling
            // if p.norm_l2() >= 2.0 * sigma * (self.a_1.width() as f64).sqrt() {
            //     // reject and try again
            //     continue;
            // }
            return Ok((d, t, z));
        }
    }

    /// Generate a NIZK argument of linear relation between two commitments.
    /// Argues `c_a = g * c_b`.
    pub fn try_open_linear_zk<R: Rng>(
        (c_a, r_a): (Self, &Vector<Polynomial<N, E>>),
        (c_b, r_b): (Self, &Vector<Polynomial<N, E>>),
        g: Polynomial<N, E>,
        rng: &mut R,
    ) -> Result<BDLOPLinearNIZKArg<N, E>> {
        assert!(c_a.a_1 == c_b.a_1);
        assert!(c_a.a_2 == c_b.a_2);
        let a_1 = &c_a.a_1;
        let a_2 = &c_a.a_2;
        loop {
            let y_1 = (0..a_1.width())
                .map(|_| Polynomial::sample_gaussian(SIGMA, rng))
                .collect::<Vector<_>>();
            let y_2 = (0..a_1.width())
                .map(|_| Polynomial::sample_gaussian(SIGMA, rng))
                .collect::<Vector<_>>();
            let t_1 = a_1.as_ref() * &y_1;
            let t_2 = a_1.as_ref() * &y_2;

            let u = (a_2.as_ref() * &y_1) * g - a_2.as_ref() * &y_2;

            let hash = blake3::hash(
                &t_1.iter()
                    .chain(t_2.iter())
                    .chain(u.iter())
                    .flat_map(|p| p.as_le_bytes())
                    .collect::<Vec<_>>(),
            );
            let mut csprng = rand_chacha::ChaCha20Rng::from_seed(hash.into());

            let d: Polynomial<N, E> =
                std::array::from_fn(|_| E::at_displacement(csprng.random_range(0..2) - 1)).into();

            let z_1 = y_1 + &(r_a.clone() * d);
            let z_2 = y_2 + &(r_b.clone() * d);
            // TODO: rejection sampling
            // if p.norm_l2() >= 2.0 * sigma * (self.a_1.width() as f64).sqrt() {
            //     // reject and try again
            //     continue;
            // }
            return Ok(BDLOPLinearNIZKArg {
                c_a,
                c_b,
                t_1,
                t_2,
                u,
                g,
                z_1,
                z_2,
            });
        }
    }
}

impl<const N: usize, E: FieldScalar> Add<&Self> for BDLOP<N, E> {
    type Output = Self;
    fn add(mut self, rhs: &Self) -> Self::Output {
        self += rhs;
        self
    }
}

impl<const N: usize, E: FieldScalar> AddAssign<&Self> for BDLOP<N, E> {
    fn add_assign(&mut self, rhs: &Self) {
        assert!(
            self.a_1 == rhs.a_1 && self.a_2 == rhs.a_2,
            "BDLOP: refusing to add commitments on different lattices"
        );
        self.c_1 += &rhs.c_1;
        self.c_2 += &rhs.c_2;
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn bdlop_commit_var_dimension() -> Result<()> {
        type Field = Seven753Scalar;
        const RING_DEGREE: usize = 16;
        let rng = &mut rand::rng();
        // just make sure our dimensions match in matrix/vector ops
        BDLOP::<RING_DEGREE, Field>::check_params()?;
        for i in 1..10 {
            let lattice = BDLOP::<64, Field>::lattice_for(i, rng);
            let (commitment, r) = BDLOP::commit(Vector::sample_uniform(i, rng), &lattice, rng);
            commitment
                .try_open(&r)
                .expect("failed to open BDLOP commitment");
            commitment
                .try_open(&(r.clone() + &r))
                .expect_err("should fail to open BDLOP commitment to bad r value");
        }
        Ok(())
    }

    #[test]
    fn bdlop_open_zk() -> Result<()> {
        type Field = MilliScalar;
        const RING_DEGREE: usize = 128;
        let rng = &mut rand::rng();
        let msg_len = 1;
        let lattice = BDLOP::<RING_DEGREE, Field>::lattice_for(msg_len, rng);
        let msg = Vector::sample_uniform(msg_len, rng);

        let (c, r) = BDLOP::commit(msg, &lattice, rng);
        let (d, t, z) = c.try_open_zk(&r, rng)?;

        for p in &z {
            let max_l2 = 4.0 * SIGMA * (RING_DEGREE as f64).sqrt();
            assert!(p.norm_l2() < max_l2);
        }
        assert_eq!(c.a_1.as_ref() * &z, t + &(c.c_1 * d));
        Ok(())
    }

    #[test]
    fn bdlop_open_linear_zk() -> Result<()> {
        type Field = MilliScalar;
        const RING_DEGREE: usize = 512;
        let rng = &mut rand::rng();
        let msg_len = 1;
        let lattice = BDLOP::<RING_DEGREE, Field>::lattice_for(msg_len, rng);
        let g = Polynomial::<RING_DEGREE, Field>::sample_uniform(rng);
        let a = Vector::sample_uniform(msg_len, rng);
        let b = a.clone() * g;

        let (c_a, r_a) = BDLOP::commit(a, &lattice, rng);
        let (c_b, r_b) = BDLOP::commit(b, &lattice, rng);

        let zk_arg = BDLOP::try_open_linear_zk((c_a, &r_a), (c_b, &r_b), g, rng)?;
        zk_arg.verify()?;

        Ok(())
    }
}
