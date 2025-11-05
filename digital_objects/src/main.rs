use anyhow::Result;
use plotters::prelude::*;

use lettuce::*;

type E = CoolScalar;
const OUT_FILE_NAME: &str = "objects/test0.gif";

fn main() -> Result<()> {
    let rng = &mut rand::rng();
    let lattice = BDLOP::<64, LOLScalar>::lattice_for(1, rng);
    let val = Polynomial::sample_uniform(rng);
    let (commitment, secret) = BDLOP::commit(val.into(), lattice, rng);
    commitment.try_open(&secret).unwrap();

    // a root of unity exists, we need to find it
    let cycle_len = 64;
    let pow = (E::Q - 1) / (2 * cycle_len);
    let mut unity = E::one();
    let mut generator = E::one();
    for g in (20..(2u128.pow(31) - 2)).rev() {
        if g < 100 {
            panic!("unable to find generator and root of unity");
        }
        let w = E::from(g as u128).modpow(pow);
        println!("{w}");
        if w.modpow(cycle_len) == E::negone() && w.modpow(2 * cycle_len) == E::one() {
            // our generator is g and w is our root of unity
            println!(
                "Q: {} degree: {cycle_len} generator: {} root: {}",
                E::CARDINALITY,
                g,
                w
            );
            generator = (g as u128).into();
            unity = w;
            break;
        }
    }

    let root = BitMapBackend::gif(OUT_FILE_NAME, (700, 700), 30)?.into_drawing_area();

    let proj = Matrix::<E>::random(3, 3, rng);
    let midpoint = (E::Q / 2) as f64;

    let raw_points = (0..=(2 * cycle_len)).map(|i| {
        let i = i % (2 * cycle_len);
        let x = unity.modpow(i);
        let v = Vector::from(vec![x, x, x]);
        let o = &proj * &v;
        let percent = (i as f64) / (2 * cycle_len) as f64;
        (
            // percent * E::Q as f64,
            unity.modpow(i).displacement() as f64 + midpoint,
            (E::from(2u128) * unity.modpow(i)).displacement() as f64 + midpoint,
            (E::from(4u128) * unity.modpow(i)).displacement() as f64 + midpoint,
            // v.clone() * &vec![E::from(1u128), E::from(2u128), E::from(3u128)].into(),
            // v.clone() * &vec![E::from(3u128), E::from(1u128), E::from(2u128)].into(),
            // v.clone() * &vec![E::from(2u128), E::from(3u128), E::from(1u128)].into(),
            // unity.modpow(i).displacement() as f64 + midpoint,
            // percent * E::Q as f64,
            // midpoint, // generator.modpow(i).displacement() as f64,
            // ((E::Q / (2 * cycle_len)) * i) as f64,
        )
    });
    let lattice_points = (0..=(2 * cycle_len)).map(|i| {
        let i = i % (2 * cycle_len);
        let x = unity.modpow(i);
        let v = Vector::from(vec![x, x, x]);
        let o = &proj * &v;
        (
            o[0].displacement() as f64 + midpoint,
            o[1].displacement() as f64 + midpoint,
            o[2].displacement() as f64 + midpoint,
        )
    });
    for x in (0..(1000 * cycle_len)) {
        root.fill(&WHITE)?;

        let mut chart = ChartBuilder::on(&root)
            .caption(
                "primitive root of unity projected by scaled identity lattice (blue) vs random lattice (red) in Z/(3*2^30+1) aka CoolScalar",
                ("sans-serif", 16),
            )
            .build_cartesian_3d(0.0..2. * midpoint, 0.0..2. * midpoint, 0.0..2. * midpoint)?;

        chart.with_projection(|mut p| {
            // p.pitch = ((x % 360) as f64 / 36.0).abs();
            p.yaw = 2. * ((x as f64) / 40.0).sin() / 4.0;
            p.pitch = 2. * ((x as f64) / 80.).sin() / 4.0;
            p.scale = 0.7;
            // p.pitch = (x as f64 / 5000.);
            // p.yaw = (((x - 120) % 180) as f64 / 36.0).abs();
            // p.pitch = (x as f64 / 50.0).abs();
            // p.yaw = (1.0 + x as f64 / 50.0).abs();
            // p.scale = 0.7;
            // p.yaw = 1.57 - (1.57 - pitch as f64 / 10.0).abs();
            p.into_matrix() // build the projection matrix
        });
        chart.draw_series(LineSeries::new(raw_points.clone(), &BLUE))?;

        chart.draw_series(PointSeries::of_element(
            raw_points.clone(),
            3,
            &BLUE,
            &|c, s, st| EmptyElement::at(c) + Circle::new((0, 0), s, st.filled()),
        ))?;

        chart.draw_series(LineSeries::new(lattice_points.clone(), &RED))?;
        chart.draw_series(PointSeries::of_element(
            lattice_points.clone(),
            3,
            &RED,
            &|c, s, st| EmptyElement::at(c) + Circle::new((0, 0), s, st.filled()),
        ))?;

        chart
            .configure_axes()
            .light_grid_style(BLACK.mix(0.15))
            .max_light_lines(3)
            .draw()?;

        root.present()?;
    }

    // To avoid the IO failure being ignored silently, we manually call the present function
    root.present().expect("Unable to write result to file, please make sure 'plotters-doc-data' dir exists under current dir");
    println!("Result has been saved to {}", OUT_FILE_NAME);

    Ok(())
}
