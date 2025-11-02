use crate::*;

use anyhow::Result;
use rand::SeedableRng;

pub trait ElementHasher<E: Element> {
    fn finish(&self) -> E;
    fn write(&mut self, bytes: &[u8]);
}

/// An implementation of Baum et. al. commitments over a scalar field.
/// https://eprint.iacr.org/2016/997.pdf
///
#[derive(Clone, Debug)]
pub struct BDLOPScalar<E: Element> {
    a_1: Matrix<E>,
    a_2: Matrix<E>,
    c_1: Vector<E>,
    c_2: Vector<E>,
}

impl<E: Element> BDLOPScalar<E> {
    /// Given a message length determine the dimension of a commitment matrix.
    ///
    /// BDLOP commitments essentially vertically compose an SIS commitment to 0 with an SIS
    /// commitment to a message vector. We need two lattice bases, A_1 and A_2. The contents of A_2
    /// requires that the total commitment width is > msg_len + 0_len, otherwise the final portion
    /// of the commitment (c_2) will output the message vector. To account for this we expand the
    /// commitment width to 3 * msg_len to fully shift the message vector away from the identity
    /// matrix in A_2.
    ///
    /// Pages 11 and 10 of https://eprint.iacr.org/2016/997.pdf
    fn dimension(msg_len: usize) -> (usize, usize) {
        // approx n*log(q) for the message length only
        // we want our message to be mixed with at least this many vectors of random elements
        //
        // we subtract msg_len because the A_2 matrix will provide an additional msg_len
        // vectors of mixing
        let a_1_height = msg_len * E::BIT_WIDTH - msg_len;
        // the A_2 matrix always has height equal to the message length
        //
        // our width is equal to the height of the A_1 matrix plus 2 * msg_len
        // we need to shift our message by 2x to create mixing elements in the A_2 matrix
        // otherwise the message is output in the plain by the identity component of A_2
        let width = a_1_height + 2 * msg_len;
        (a_1_height, width)
    }

    pub fn lattice_for<R: Rng>(msg_len: usize, rng: &mut R) -> (Matrix<E>, Matrix<E>) {
        let (a_1_height, width) = Self::dimension(msg_len);
        // the A_1 lattice base
        let a_1 = Matrix::<E>::identity(a_1_height).compose_horizontal(Matrix::random(
            a_1_height,
            width - a_1_height,
            rng,
        ));

        // the A_2 lattice base
        let a_2 = Matrix::<E>::zero(msg_len, a_1_height)
            .compose_horizontal(Matrix::identity(msg_len))
            .compose_horizontal(Matrix::random(msg_len, width - a_1_height - msg_len, rng));
        (a_1, a_2)
    }

    /// Generate a BDLOP commitment to a vector of scalar elements.
    pub fn commit<R: Rng>(
        val: Vector<E>,
        lattice: (Matrix<E>, Matrix<E>),
        rng: &mut R,
    ) -> (Vector<E>, Self) {
        let (a_1, a_2) = lattice;

        // the secret committing to the zero component
        let r = (0..a_1.width())
            .map(|_| {
                let sample: i32 = rng.random_range(0..=2) - 1;
                E::at_displacement(sample)
            })
            .collect::<Vec<_>>()
            .into();

        let c_1 = &a_1 * &r;
        let c_2 = &a_2 * &r + &val;

        (r, Self { a_1, a_2, c_1, c_2 })
    }

    /// Attempt to open a commitment directly using the stored `r` value.
    /// First attempts to open c_1 to the zero vector. If this succeeds c_2 is opened to whatever
    /// value is committed. We assume that if r_1 is correct, r_2 is correct.
    pub fn try_open(&self, secret: (&Vector<E>, &Vector<E>)) -> Result<Vector<E>> {
        let (r_1, r_2) = secret;
        if &self.a_1 * r_1 != self.c_1 {
            anyhow::bail!("Failed to open commitment, secret is incorrect");
        }
        Ok(&self.c_2 - &self.a_2 * r_2)
    }

    /// Attempt to generate a non-interactive ZK proof of opening.
    ///
    /// Described on page 15 of https://eprint.iacr.org/2016/997.pdf
    ///
    /// This implementation modifies the d value to be a single ternary element. This challenge is
    /// not meaningfully sound.
    pub fn try_open_zk<R: Rng>(
        &self,
        r: &Vector<E>,
        rng: &mut R,
    ) -> Result<(E, Vector<E>, Vector<E>)> {
        let sigma = 5f64;
        let y = GaussianCDT::cache_or_init::<E>(sigma).sample_vec(self.a_1.width(), rng);
        let t = &self.a_1 * &y;
        let hash = blake3::hash(
            &t.iter()
                .flat_map(|v: &E| v.as_le_bytes())
                .collect::<Vec<u8>>(),
        );
        let mut csprng = rand_chacha::ChaCha20Rng::from_seed(hash.into());

        let z;
        let d;
        loop {
            let d_maybe = E::at_displacement(csprng.random_range(0..2) - 1);
            let z_maybe = y.clone() + &(r.clone() * d_maybe);
            if z_maybe.norm_l2() >= 2.0 * sigma * (self.a_1.width() as f64).sqrt() {
                // reject and try again
                continue;
            }
            z = z_maybe;
            d = d_maybe;
            break;
        }
        Ok((d, t, z))
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn bdlop_commit_var_dimension() {
        type Field = OxfoiScalar;
        let rng = &mut rand::rng();
        // just make sure our dimensions match in matrix/vector ops
        for i in 1..10 {
            let lattice = BDLOPScalar::lattice_for(i, rng);
            let _ = BDLOPScalar::<Field>::commit(Vector::random(i, rng), lattice, rng);
        }
    }

    #[test]
    fn bdlop_open_zk() -> Result<()> {
        type Field = OxfoiScalar;
        let rng = &mut rand::rng();
        let msg_len = 5;
        let lattice = BDLOPScalar::lattice_for(msg_len, rng);
        let msg = Vector::random(msg_len, rng);

        let (r, c) = BDLOPScalar::<Field>::commit(msg, lattice, rng);
        let (d, t, z) = c.try_open_zk(&r, rng)?;

        assert_eq!(&c.a_1 * &z, t + &(c.c_1 * d));
        Ok(())
    }
}
