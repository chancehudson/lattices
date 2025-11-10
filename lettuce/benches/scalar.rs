use std::hint::black_box;

use criterion::*;

use lettuce::*;

fn criterion_benchmark(c: &mut Criterion) {
    type E = MilliScalarMont;
    // MilliScalar::sample_uniform(rng)
    c.bench_function("MilliScalar mul", |b| {
        b.iter_batched(
            || {
                let rng = &mut rand::rng();
                (0..1000)
                    .map(|_| (E::sample_uniform(rng), E::sample_uniform(rng)))
                    .collect::<Vec<_>>()
            },
            |data| {
                for (lhs, rhs) in data {
                    black_box(black_box(lhs) * black_box(rhs));
                }
            },
            BatchSize::SmallInput,
        )
    });
    c.bench_function("MilliScalar add", |b| {
        b.iter_batched(
            || {
                let rng = &mut rand::rng();
                (0..1000)
                    .map(|_| (E::sample_uniform(rng), E::sample_uniform(rng)))
                    .collect::<Vec<_>>()
            },
            |data| {
                for (lhs, rhs) in data {
                    black_box(black_box(lhs) + black_box(rhs));
                }
            },
            BatchSize::SmallInput,
        )
    });
    c.bench_function("MilliScalar sub", |b| {
        b.iter_batched(
            || {
                let rng = &mut rand::rng();
                (0..1000)
                    .map(|_| (E::sample_uniform(rng), E::sample_uniform(rng)))
                    .collect::<Vec<_>>()
            },
            |data| {
                for (lhs, rhs) in data {
                    black_box(black_box(lhs) - black_box(rhs));
                }
            },
            BatchSize::SmallInput,
        )
    });
    c.bench_function("MilliScalar inv", |b| {
        b.iter_batched(
            || {
                let rng = &mut rand::rng();
                (0..1000)
                    .map(|_| (E::sample_uniform(rng), E::sample_uniform(rng)))
                    .collect::<Vec<_>>()
            },
            |data| {
                for (lhs, rhs) in data {
                    black_box(black_box(lhs).inverse());
                    black_box(black_box(rhs).inverse());
                }
            },
            BatchSize::SmallInput,
        )
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
