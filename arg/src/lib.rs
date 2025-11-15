mod witness;

use std::marker::PhantomData;
use std::sync::Arc;

use anyhow::Result;
use lettuce::FieldScalar;
use serde::Deserialize;
use serde::Serialize;
use zkpo::*;

use lettuce::MilliScalarMont;

#[derive(Clone, Copy, Serialize, Deserialize, Debug)]
pub struct Pow5Arg<E: FieldScalar = MilliScalarMont> {
    input: u64,
    _phantom: PhantomData<E>,
}
//
// impl<E: FieldScalar> zkpo::ZKIO<E> for Pow5Arg<E> {
//     fn serialize(&self) -> impl Iterator<Item = E> {
//         unimplemented!()
//     }
//
//     fn deserialize(from: impl Iterator<Item = E>) -> Self {
//         unimplemented!()
//     }
// }
//
// impl<E: FieldScalar> zkpo::ZKArg for Pow5Arg<E> {
//     fn program_id(&self) -> Arc<[u8; 32]> {}
//     fn program(&self) -> Option<impl ZKProgram<E>> {}
//     fn cipher_bytes(&self) -> &[u8] {}
// }
//
// impl<E: FieldScalar> ZKProgram<E> for Pow5Arg<E> {
//     fn elf(&self) -> Arc<&[u8]> {
//         unimplemented!()
//     }
//
//     fn execute(&self, input: impl ZKIO<E>, program: impl ZKProgram<E>) -> Result<impl ZKArg<E>> {
//         // let mut csprng = rand_chacha::ChaCha20Rng::from_seed(rand::random::<[u8; 32]>());
//         // let rng = &mut csprng;
//         // type E = MilliScalarMont;
//         // let (r1cs, wtns) = R1CS::<E>::sample_uniform(2 * 8192, 2 * 8192, rng);
//         // let arg = HiddenR1CS::<1024, _>::commit(wtns, r1cs, rng)?;
//         unimplemented!()
//     }
//
//     fn verify(&self, arg: impl ZKArg<E>) -> Result<impl ZKIO<E>> {
//         arg.verify()?;
//         unimplemented!()
//     }
// }
