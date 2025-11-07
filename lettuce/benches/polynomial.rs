use std::hint::black_box;

use criterion::*;

use lettuce::*;

fn criterion_benchmark(c: &mut Criterion) {
    type E = MilliScalar;
    const N: usize = 64;
    type P = Polynomial<N, E>;
    c.bench_function("Polynomial mul", |b| {
        b.iter_batched(
            || {
                let rng = &mut rand::rng();
                (0..1000)
                    .map(|_| (P::sample_uniform(rng), P::sample_uniform(rng)))
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
    c.bench_function("Polynomial add", |b| {
        b.iter_batched(
            || {
                let rng = &mut rand::rng();
                (0..1000)
                    .map(|_| (P::sample_uniform(rng), P::sample_uniform(rng)))
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
    c.bench_function("Polynomial sub", |b| {
        b.iter_batched(
            || {
                let rng = &mut rand::rng();
                (0..1000)
                    .map(|_| (P::sample_uniform(rng), P::sample_uniform(rng)))
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
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
