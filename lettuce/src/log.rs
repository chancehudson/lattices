/// Standardize and identify library logs.
///
/// This exposes a macro api identical to the public api of `log`.
/// Use log::info!(), log::error!(), etc as normal
#[cfg(test)]
pub(crate) use test::SETUP;

#[cfg(test)]
mod test {
    use std::sync::LazyLock;

    pub static SETUP: LazyLock<()> = LazyLock::new(|| {
        // env_logger::builder()
        //     .filter_level(log_ext::LevelFilter::Debug)
        //     .init();
    });
}

pub(crate) const PREFIX: &'static str = "🐰🥬 ";

macro_rules! error {
    ($($input:expr),*) => {
        #[cfg(test)]
        std::sync::LazyLock::<()>::force(&log::SETUP);
        log_ext::error!("{}{}", crate::log::PREFIX, format!($($input,)*));
    };
}

macro_rules! info {
    ($($input:expr),*) => {
        #[cfg(test)]
        std::sync::LazyLock::<()>::force(&log::SETUP);
        log_ext::info!("{}{}", crate::log::PREFIX, format!($($input,)*));
    };
}

macro_rules! debug {
    ($($input:expr),*) => {
        #[cfg(test)]
        std::sync::LazyLock::<()>::force(&log::SETUP);
        log_ext::debug!("{}{}", crate::log::PREFIX, format!($($input,)*));
    };
}

pub(crate) use debug;
pub(crate) use error;
pub(crate) use info;
