# hidden\_r1cs

First consider blindfold: [blindfoldzkp.github.io/slide22.html](https://blindfoldzkp.github.io/slide22.html)

The notation in the link is slightly inconsistent, but we have the following structures:

- R1CS: $\textbf{AwBw} = \textbf{Cw}$ with $\textbf{A}, \textbf{B}, \textbf{C} \in \mathbb{Z}_q^{m \times n}$ and $\textbf{w} \in \mathbb{Z}_q^n$ with $m$ constraints and $n$ variables.
- Relaxed R1CS (aka R3CS): $\textbf{AwBw} = s\textbf{Cw} + \textbf{e}$. As above. $s \in \mathbb{Z}_q$, $\textbf{e} \in \mathbb{Z}_q^m$.

Relaxed R1CS adds a constraint to the consistency of the system itself. A scalar constraint on the witness is, in some cases, an effective way to constrain a witness. However consider the case where variables are relatively independent. Breaking a 32 bit scalar constraint may be < 2^50 work. This can be accounted for at the arithmetization level by adding non-sparse safety constraints.

Each safety constraint is a vector containing non-identity scalars. The vector is appended as a row to the R1CS $\textbf{A}$ and $\textbf{C}$ matrices, and a vector with a single $1$ followed by $m - 1$ trailing $0$'s is appended to the $\textbf{B}$ matrix.

Indeed this can be applied to sparse R1CS instances by ensuring some minimum number of variables are constrained in non-additive ways.

#### blindfold relation

$\textbf{w}' \overset{\\$}{\leftarrow} \mathbb{Z}_q^n$ sample a random witness mask vector.

An R3CS instance is created using $\textbf{w}'$ (and $s = 1$ omitted):

$\textbf{e}' \leftarrow \textbf{Aw}'\textbf{Bw}' - \textbf{Cw}' \in \mathbb{Z}_q^m$

A final R3CS relation is formed using the above R1CS and R3CS:

$\textbf{e} \leftarrow \textbf{Aw}' \textbf{Bw} + \textbf{AwBw}' - \textbf{Cw} - \textbf{Cw}'$

$\textbf{e} + \textbf{e}' = \textbf{A}(\textbf{w} + \textbf{w}')\textbf{B}(\textbf{w} + \textbf{w}') - \textbf{C}(\textbf{w} + \textbf{w}')$

$c \overset{\\$}{\leftarrow} \mathbb{Z}_q$ sample a challenge scalar from commitments to $\textbf{w}$, $\textbf{w}'$, $\textbf{e}$ and $\textbf{e}'$.

A final public witness is formed:

$\textbf{w}_f \overset{\pi}{\leftarrow} \textbf{w}' + c\textbf{w}$

A final public crossterm is formed:

$\textbf{e}_f \overset{\pi}{\leftarrow} \textbf{e}' + c \textbf{e}$

Commitments to $\textbf{e}$, $\textbf{e}'$, $\textbf{w}$ and $\textbf{w}'$ must argue existence of $\textbf{e}_f$ and $\textbf{w}_f$ as above.

Verification requires checking the following:

- $c \overset{\\$}{\leftarrow} \textbf{w}, \textbf{w}', \textbf{e}, \textbf{e}'$
- $\textbf{w}_f^{\pi}$
- $\textbf{e}_f^{\pi}$
- $\textbf{e}_f = \textbf{Aw}_f + \textbf{Bw}_f - (1 + c)\textbf{Cw}_f$

This requires linear homomorphism only. Recursive operations require binding only (not hiding). Statistical hiding is possible using large lattice bases.

This system transforms an R1CS into an R3CS over the same constraint matrices with a public token witness.

=== RKCS digression

Is it possible to put some set of constraints `k` on the system? At first look the dimensions of the structures necessitate the use of a scalar and a vector (`s` and `E`).

An R1CS could be segmented into `k` groups of `m / k` constraints each. Each segment would be argued independently with final relation shown behind a lattice construction. Each zk argument of relation requires linear data in the witness length, so more outer (non-program?) variables = bigger proofs/longer proving time.

Proofs and proving refer to construction of arguments as information. Constructions are arguments unless stated otherwise.

A similar question I've been pondering: we can make arguments that a program was executed in the past.

Can we make arguments that a program was executed in the future?

Provided time is not totally ordered.

===
