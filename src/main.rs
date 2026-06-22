//! ohc_ingest driver — processes exactly one layer per run.
//!
//! Usage:
//!   ohc_ingest [config.toml] --tag NAME --layer T-B --years Y|A:B --months M|A:B|M,M,...
//!
//! `--tag`, `--layer`, `--years`, `--months` are REQUIRED (one run = one layer over one time
//! slice). Each may also be given via env (`OHC_TAG`, `OHC_LAYER`, `OHC_YEARS`, `OHC_MONTHS`);
//! CLI wins. `--tag` is the run identifier and labels the output store + metadata.
//! Static constants + paths come from `config.toml`, or from the defaults + path env vars
//! (`OHC_DIR_MEAN`, `OHC_DIR_ENSEMBLE`, `OHC_DIR_OUT`, `OHC_ETOPO`, `OHC_BASINMASK`) when no
//! config is given.
//!
//! Examples:
//!   ohc_ingest --layer 15-20 --years 2016 --months 8
//!   ohc_ingest config.toml --layer 300-700 --years 2004:2025 --months 1:12
//!
//! To process many layers, run one invocation per layer (e.g. a scheduler job array).

use anyhow::{bail, Context, Result};
use std::env;
use std::path::PathBuf;
use std::time::Instant;

use ohc_ingest::config::{parse_layer, parse_months, parse_years, LayerSpec, RunConfig, Slice};
use ohc_ingest::{basinmask, grid as gridmod, ingest, masks, ncread, zarrwrite, GridDef};

struct Cli {
    config_path: Option<String>,
    tag: Option<String>,
    layer: Option<LayerSpec>,
    years: Option<[i32; 2]>,
    months: Option<Vec<u32>>,
}

fn parse_cli() -> Result<Cli> {
    let args: Vec<String> = env::args().skip(1).collect();
    let mut cli = Cli { config_path: None, tag: None, layer: None, years: None, months: None };
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--tag" => {
                i += 1;
                cli.tag = Some(args.get(i).context("--tag needs a value")?.clone());
            }
            "--layer" => {
                i += 1;
                cli.layer = Some(parse_layer(args.get(i).context("--layer needs a value")?)?);
            }
            "--years" => {
                i += 1;
                cli.years = Some(parse_years(args.get(i).context("--years needs a value")?)?);
            }
            "--months" => {
                i += 1;
                cli.months = Some(parse_months(args.get(i).context("--months needs a value")?)?);
            }
            s if s.starts_with("--") => bail!("unknown flag {s}"),
            s => cli.config_path = Some(s.to_string()),
        }
        i += 1;
    }
    Ok(cli)
}

fn load_run_config(config_path: &Option<String>) -> Result<RunConfig> {
    if let Some(p) = config_path {
        let text = std::fs::read_to_string(p).with_context(|| format!("reading config {p}"))?;
        Ok(toml::from_str(&text)?)
    } else {
        let mut c = RunConfig::defaults();
        if let Some(v) = env::var_os("OHC_DIR_MEAN") { c.dir_mean = PathBuf::from(v); }
        if let Some(v) = env::var_os("OHC_DIR_ENSEMBLE") { c.dir_ensemble = PathBuf::from(v); }
        if let Some(v) = env::var_os("OHC_DIR_OUT") { c.dir_out = PathBuf::from(v); }
        if let Some(v) = env::var_os("OHC_ETOPO") { c.etopo_path = PathBuf::from(v); }
        if let Some(v) = env::var_os("OHC_BASINMASK") { c.basinmask_path = PathBuf::from(v); }
        Ok(c)
    }
}

fn resolve_slice(cli: &Cli) -> Result<Slice> {
    let layer = match cli.layer {
        Some(l) => l,
        None => match env::var("OHC_LAYER") {
            Ok(v) => parse_layer(&v)?,
            Err(_) => bail!("--layer is required (one layer per run), e.g. --layer 15-20"),
        },
    };
    let years = match cli.years {
        Some(y) => y,
        None => match env::var("OHC_YEARS") {
            Ok(v) => parse_years(&v)?,
            Err(_) => bail!("--years is required, e.g. --years 2016 or --years 2004:2025"),
        },
    };
    let months = match &cli.months {
        Some(m) => m.clone(),
        None => match env::var("OHC_MONTHS") {
            Ok(v) => parse_months(&v)?,
            Err(_) => bail!("--months is required, e.g. --months 8 or --months 1:12"),
        },
    };
    Ok(Slice { layer, years, months })
}

fn main() -> Result<()> {
    let cli = parse_cli()?;
    let mut cfg = load_run_config(&cli.config_path)?;
    cfg.run_tag = match &cli.tag {
        Some(t) => t.clone(),
        None => match env::var("OHC_TAG") {
            Ok(v) => v,
            Err(_) => bail!("--tag is required (run identifier), e.g. --tag OP20260110"),
        },
    };
    let slice = resolve_slice(&cli)?;
    let grid = GridDef::mapping();

    eprintln!(
        "Run {}: layer {} dbar, years {}..={}, months {:?}, {} timestep(s)",
        cfg.run_tag,
        slice.layer.tag(),
        slice.years[0],
        slice.years[1],
        slice.months,
        slice.time_axis().len(),
    );

    match &cli.config_path {
        Some(p) => eprintln!("config: {p}"),
        None => eprintln!("config: none — using built-in defaults + env vars"),
    }
    eprintln!("  dir_mean     = {}", cfg.dir_mean.display());
    eprintln!("  dir_ensemble = {}", cfg.dir_ensemble.display());
    eprintln!("  dir_out      = {}", cfg.dir_out.display());
    eprintln!("  etopo        = {}", cfg.etopo_path.display());
    eprintln!("  basinmask    = {}", cfg.basinmask_path.display());

    eprintln!("Loading ancillary grids…");
    let etopo = ncread::read_etopo(&cfg.etopo_path, &grid)
        .with_context(|| format!("reading etopo {}", cfg.etopo_path.display()))?;
    let basin_id = basinmask::read_basin_id(&cfg.basinmask_path, &grid)
        .with_context(|| format!("reading basin mask {}", cfg.basinmask_path.display()))?;
    let cell_area = gridmod::cell_area(&grid);

    let t0 = Instant::now();
    eprintln!(">>> layer {} dbar", slice.layer.tag());
    let data = ingest::ingest_layer(&cfg, &slice, &grid)
        .with_context(|| format!("ingesting layer {}", slice.layer.tag()))?;
    let (never, incomplete) = masks::compute_validity(&data.ohc_mean);
    let flags = masks::build_flags(
        &grid, &slice.layer, &etopo, Some(&basin_id),
        cfg.latitude_range_to_keep, &cfg.basins_to_remove, &never, &incomplete,
    )?;
    let store = zarrwrite::write_layer_store(
        &cfg, &slice, &grid, &data, &flags, &etopo, &basin_id, &cell_area,
    )?;
    eprintln!("wrote {} ({:.1}s)", store.display(), t0.elapsed().as_secs_f64());
    Ok(())
}
