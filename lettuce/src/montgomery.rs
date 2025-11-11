/// Montgomery representation of a finite field with modulus < 2^32
/// We hardcode types such that R is 2^32
pub(crate) struct Montgomery32 {
    q: u32,
    q_inv: u32, // -q^(-1) % 2^32
    r_sq: u32,  // R^2 % q
}

impl Montgomery32 {
    /// Uses newton-raphson modular inversion for q^-1
    ///
    /// See "Inverses modulo prime powers (including powers of 2)"
    /// https://en.wikipedia.org/wiki/Modular_multiplicative_inverse
    pub const fn new(q: u32) -> Self {
        let mut q_inv = q;
        // 5 = log2(32)
        let mut i = 0;
        while i < 5 {
            q_inv = q_inv.wrapping_mul(2u32.wrapping_sub(q.wrapping_mul(q_inv)));
            i += 1;
        }
        q_inv = q_inv.wrapping_neg();
        assert!(u32::MAX == q_inv.wrapping_mul(q));
        let r_q = (1u64 << 32) % q as u64;
        let r_sq = ((r_q * r_q) % q as u64) as u32;
        assert!(r_sq as u128 == (1u128 << 64) % q as u128);
        Self { q, q_inv, r_sq }
    }

    #[inline(always)]
    pub fn reduce(&self, t: u64) -> u32 {
        let m = (t as u32).wrapping_mul(self.q_inv);

        let mq = (m as u64) * (self.q as u64);

        let (sum, overflowed) = t.overflowing_add(mq);

        let mut reduced = sum >> 32;
        if overflowed {
            reduced += 1u64 << 32;
        }

        debug_assert!(
            reduced < 2 * (self.q as u64),
            "Montgomery reduced value is too large"
        );
        if let Some(out) = reduced.checked_sub(self.q as u64) {
            out as u32
        } else {
            reduced as u32
        }
    }

    #[inline(always)]
    pub fn mul(&self, a: u32, b: u32) -> u32 {
        debug_assert!(
            a < self.q && b < self.q,
            "Montgomery mul inputs must be < q"
        );
        self.reduce(a as u64 * b as u64)
    }

    #[inline(always)]
    pub fn neg(&self, a: u32) -> u32 {
        debug_assert!(a < self.q, "Montgomery neg inputs must be < q");
        if a == 0 {
            return 0;
        }
        self.q - a
    }

    #[inline(always)]
    pub fn add(&self, a: u32, b: u32) -> u32 {
        debug_assert!(
            a < self.q && b < self.q,
            "Montgomery add inputs must be < q"
        );
        if b == 0 {
            return a;
        }
        let sum = a as u64 + b as u64;
        if let Some(out) = sum.checked_sub(self.q as u64) {
            out as u32
        } else {
            sum as u32
        }
    }

    // Convert a value to montgomery form. Input must be < q.
    #[inline(always)]
    pub fn to_mont(&self, a: u32) -> u32 {
        debug_assert!(a < self.q, "Montgomery input is not in the prime field");
        self.mul(a, self.r_sq)
    }

    // Convert from Montgomery form
    #[inline(always)]
    pub fn from_mont(&self, a: u32) -> u32 {
        self.reduce(a as u64)
    }
}

#[test]
fn mont_roundtrip() {
    let q = 455 * 2u32.pow(20) * 9 + 1;
    let mont = Montgomery32::new(q);
    for i in 0..1000 {
        let a_roundtrip = mont.from_mont(mont.to_mont(i));
        assert_eq!(i, a_roundtrip);
    }
    for _ in 0..1000 {
        let a = rand::random::<u32>() % q;
        let a_roundtrip = mont.from_mont(mont.to_mont(a));
        assert_eq!(a, a_roundtrip);
    }
}

#[test]
fn mont_muls() {
    let q = 2u32.pow(31) - 1;
    let mont = Montgomery32::new(q);
    for _ in 0..1000 {
        let a = rand::random::<u32>() % q;
        let b = rand::random::<u32>() % q;
        let c = ((a as u64 * b as u64) % q as u64) as u32;
        let a = mont.to_mont(a);
        let b = mont.to_mont(b);
        let c_mont = mont.mul(a, b);
        assert_eq!(c, mont.from_mont(c_mont));
    }
}

#[test]
fn mont_adds() {
    let q = 2u32.pow(31) - 1;
    let mont = Montgomery32::new(q);
    for _ in 0..1000 {
        let a = rand::random::<u32>() % q;
        let b = rand::random::<u32>() % q;
        let c = ((a as u64 + b as u64) % q as u64) as u32;
        let a = mont.to_mont(a);
        let b = mont.to_mont(b);
        let c_mont = mont.add(a, b);
        assert_eq!(c, mont.from_mont(c_mont));
    }
}

#[test]
fn mont_subs() {
    let q = 2u32.pow(31) - 1;
    let mont = Montgomery32::new(q);
    for _ in 0..1000 {
        let a = rand::random::<u32>() % q;
        let b = rand::random::<u32>() % q;
        let c = ((a as u64 + (2 * q as u64 - b as u64)) % q as u64) as u32;
        let a = mont.to_mont(a);
        let b = mont.to_mont(b);
        let c_mont = mont.add(a, mont.neg(b));
        assert_eq!(c, mont.from_mont(c_mont));
    }
}

#[test]
fn mont_basic_properties() {
    let q = 2u32.pow(31) - 1;
    let mont = Montgomery32::new(q);

    // Test 0
    let zero = mont.to_mont(0);
    assert_eq!(zero, 0, "Montgomery form of 0 should be 0");
    assert_eq!(mont.from_mont(zero), 0);

    // Test 1
    let one = mont.to_mont(1);
    println!("Montgomery form of 1: {}", one);
    assert_eq!(mont.from_mont(one), 1);

    // Test identity: a * 1 = a
    for _ in 0..10 {
        let a = rand::random::<u32>() % q;
        let a_mont = mont.to_mont(a);
        let result = mont.mul(a_mont, one);
        assert_eq!(mont.from_mont(result), a, "a * 1 should equal a");
    }

    // Test: convert to mont and back
    for _ in 0..100 {
        let a = rand::random::<u32>() % q;
        let a_mont = mont.to_mont(a);
        let a_back = mont.from_mont(a_mont);
        assert_eq!(a, a_back, "Round trip failed for {}", a);
    }
}

#[test]
fn mont_r_squared_check() {
    let q = 2u32.pow(31) - 1;
    let mont = Montgomery32::new(q);

    // Verify r_squared is computed correctly
    let r_mod_q = ((1u64 << 32) % q as u64) as u32;
    let expected_r_sq = ((r_mod_q as u64 * r_mod_q as u64) % q as u64) as u32;

    println!("q: {}", q);
    println!("R mod q: {}", r_mod_q);
    println!("R^2 mod q (expected): {}", expected_r_sq);
    println!("R^2 mod q (computed): {}", mont.r_sq);

    assert_eq!(mont.r_sq, expected_r_sq);
}

#[test]
fn mont_q_inv_check() {
    let q = 2u32.pow(31) - 1;
    let mont = Montgomery32::new(q);

    let product = q.wrapping_mul(mont.q_inv);
    println!("q: {}", q);
    println!("q_inv: {}", mont.q_inv);
    println!("q * q_inv: {:#x} (should be 0xFFFFFFFF)", product);

    assert_eq!(product, u32::MAX);
}
