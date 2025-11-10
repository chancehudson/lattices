# hidden\_r1cs

First consider blindfold: [blindfoldzkp.github.io/slide22.html](https://blindfoldzkp.github.io/slide22.html)

The notation in the link is slightly inconsistent, but we have the following structures:

- `AwBw = Cw`: R1CS. `A`, `B`, `C` are matrices with `m` rows and `n` columns. `m =` number of constraints. `n =` number of variables. `w = ` witness (aka program variables).
- `AwBw = sCw - E`: Relaxed R1CS (aka RR1CS aka R3CS). As above. `s =` scalar variable, `E = ` vector of dimension `m`.

Relaxed R1CS adds a constraint to the consistency of the system itself. A scalar constraint on the witness is, in some cases, an effective way to constrain a witness. Consider the case where variables are relatively independent, breaking a 32 bit scalar constraint may be < 2^50 work. However, this can be accounted for at the arithmetization level by adding non-sparse safety constraints.

A vector of safety constraints may be introduced. Each constraint is a vector containing non-identity scalars. The vector is appended to the R1CS `A` and `C` matrices, and a zero vector is appended to the `B` matrix.

Indeed this can be applied to sparse R1CS instances by ensuring some minimum number of variables are constrained in non-additive ways.

=== RKCS digression

Is it possible to put some set of constraints `k` on the system? At first look the dimensions of the structures necessitate the use of a scalar and a vector (`s` and `E`).

An R1CS could be segmented into `k` groups of `m / k` constraints each. Each segment would be argued independently with final relation shown behind a lattice construction. Each zk argument of relation requires linear data in the witness length, so more outer (non-program?) variables = bigger proofs/longer proving time.

Proofs and proving refer to construction of arguments as information. Constructions are arguments unless stated otherwise.

A similar question I've been pondering: we can make arguments that a program was executed in the past.

Can we make arguments that a program was executed in the future?

Provided time is not totally ordered.

===

`A, B, C: Matrix<Scalar>`

`w: Vector<Scalar>`

A vector `w_m =` witness mask is randomly sampled from the field.

An R3CS instance is created using `w_m` and `s = 1`:

`E <- Aw_mBw_m - sCw_m`

An R3CS is created by joining the R1CS and the R3CS above:

`E == A(w + w_m)B(w + w_m) - sC(w + w_m)`

A scalar `c =` challenge is randomly sampled from commitments to above. The system is algebraically transformed to the following relation:

- `w_f <-? w_m + c * w`
- `x <- Aw_mBw + AwBw_m - Cw - Cw_m`
- `c*x == Aw_f + Bw_f - (1 + c)Cw_f`

Commitments to `w` and `w_m` must argue existence of `w_f` such that `w_f == w_m + c * w`.

`w_f` is public knowledge.

This requires linear homomorphism only. Recursive operations require binding only (not hiding). Statistical hiding is possible using large lattice bases.

This system transforms an R1CS into an R3CS over the same constraint matrices with a hidden witness.
