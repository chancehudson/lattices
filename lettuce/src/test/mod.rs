use crate::*;

use anyhow::Result;

#[test]
fn identity_r1cs() -> Result<()> {
    let r1cs = R1CS::<OxfoiScalar>::identity(10, 10);
    let witness = Vector::new(r1cs.dimension().0);
    assert_eq!(r1cs.eval(&witness)?, witness);
    Ok(())
}
