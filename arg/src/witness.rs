use std::collections::BTreeMap;
use std::sync::Arc;

use anyhow::Result;

use lettuce::*;

/// a_coef * a_member * b_coef * b_member = c_coef * c_member
pub struct WitnessConstraint<E: FieldScalar> {
    pub relation: MultiPolynomial<E>,
    pub x_wtns: Option<Arc<WitnessMember>>,
    pub y_wtns: Option<Arc<WitnessMember>>,
    pub z_wtns: Option<Arc<WitnessMember>>,
}

impl<E: FieldScalar> WitnessConstraint<E> {
    pub fn from_wtns_members(
        x_wtns: Option<Arc<WitnessMember>>,
        y_wtns: Option<Arc<WitnessMember>>,
        z_wtns: Option<Arc<WitnessMember>>,
        relation: MultiPolynomial<E>,
    ) -> Self {
        if [&x_wtns, &y_wtns, &z_wtns]
            .iter()
            .filter(|v| v.is_some())
            .count()
            < 2
        {
            panic!("WitnessConstraint cannot exist without two variables");
        }
        Self {
            x_wtns,
            y_wtns,
            z_wtns,
            relation,
        }
    }
}

#[derive(Default)]
pub struct LettuceProgram<E: FieldScalar> {
    witness_map: BTreeMap<usize, Arc<WitnessMember>>,
    /// quadratic constraints between max 3 variables
    /// and 4 constants
    constraints: Vec<WitnessConstraint<E>>,
}

impl<E: FieldScalar> LettuceProgram<E> {
    pub fn new() -> Result<Self> {
        let mut out = Self::default();
        // zero element
        out.private_input()?;
        // one element
        out.private_input()?;
        Ok(out)
    }

    pub fn add_constraint(&mut self, constraint: WitnessConstraint<E>) -> Result<()> {
        self.constraints.push(constraint);
        Ok(())
    }

    pub fn private_intermediate(&mut self) -> Result<Arc<WitnessMember>> {
        self.add_witness(WitnessMember {
            idx: self.witness_map.len(),
            is_input: false,
            is_public: false,
            is_committed: false,
        })
    }

    pub fn public_input(&mut self) -> Result<Arc<WitnessMember>> {
        self.add_witness(WitnessMember {
            idx: self.witness_map.len(),
            is_input: true,
            is_public: true,
            is_committed: true,
        })
    }

    pub fn private_input(&mut self) -> Result<Arc<WitnessMember>> {
        self.add_witness(WitnessMember {
            idx: self.witness_map.len(),
            is_input: true,
            is_public: false,
            is_committed: true,
        })
    }

    pub fn public_output(&mut self) -> Result<Arc<WitnessMember>> {
        self.add_witness(WitnessMember {
            idx: self.witness_map.len(),
            is_input: false,
            is_public: true,
            is_committed: true,
        })
    }

    pub fn add_witness(&mut self, w: WitnessMember) -> Result<Arc<WitnessMember>> {
        let w = Arc::new(w);
        if self.witness_map.contains_key(&w.idx) {
            anyhow::bail!("Witness error: entry {} already exists!", w.idx);
        }
        self.witness_map.insert(w.idx, w.clone());
        Ok(w)
    }

    /// Get the element in the witness that is 0.
    pub fn zero(&self) -> Arc<WitnessMember> {
        self.witness_map
            .get(&0)
            .expect("First witness entry did not exist.")
            .clone()
    }

    /// Get the element in the witness that is 1.
    pub fn one(&self) -> Arc<WitnessMember> {
        self.witness_map
            .get(&1)
            .expect("Second witness entry did not exist.")
            .clone()
    }
}

#[derive(Clone, Default)]
pub struct WitnessMember {
    /// index of the element in the witness vector
    pub idx: usize,
    /// is the value an input ?
    pub is_input: bool,
    /// is the value publicly revealed
    pub is_public: bool,
    /// is the value committed for some other reason?
    pub is_committed: bool,
}

#[test]
fn build_circuit() -> Result<()> {
    type E = MilliScalarMont;
    let mut program = LettuceProgram::<MilliScalarMont>::new()?;
    let x = program.public_input()?;
    let y = program.private_input()?;
    let mut relation = MultiPolynomial::<E>::zero();
    relation.a = 20.into();
    relation.b = 50.into();
    relation.c = 100.into();

    let constraint = WitnessConstraint::from_wtns_members(Some(x), Some(y), None, relation);

    Ok(())
}
