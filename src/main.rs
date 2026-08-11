//! ohc_ingest driver — processes exactly one layer per run.
//!
//! Usage:
//!   ohc_ingest [config.toml] --tag NAME --provenance-link URL --layer T-B --years Y|A:B \
//!       --months M|A:B|M,M,... [--no-ensemble]
//!
//! `--tag`, `--provenance-link`, `--layer`, `--years`, `--months` are REQUIRED (one run = one layer
//! over one time slice). Each may also be given via env (`OHC_TAG`, `OHC_PROVENANCE_LINK`,
//! `OHC_LAYER`, `OHC_YEARS`, `OHC_MONTHS`); CLI wins. `--tag` is the run identifier: it names the
//! output store (`ohc_<tag>_plev<layer>.zarr`) and is written to the store's `provenance_tag`
//! attr (whitespace-stripped, never lowercased — must match the provenance record char-for-char).
//! `--provenance-link` points at that provenance record and is written to the `provenance_link` attr.
//! `--no-ensemble` (or `OHC_NO_ENSEMBLE`) ingests the mean only — skips the LocalCondSim files
//! and omits `ohc_ensemble` from the store (for mean-only products, or incomplete CondSim sets).
//! Static constants + paths come from `config.toml`, or from the defaults + path env vars
//! (`OHC_DIR_MEAN`, `OHC_DIR_ENSEMBLE`, `OHC_DIR_OUT`, `OHC_ETOPO`, `OHC_BASINMASK`) when no
//! config is given. `--dir_mean`, `--dir_ensemble`, `--dir_out` override those directories on the
//! command line (CLI wins over both config and env) — for munging paths per shell submission.
//! The config may omit those three dirs entirely (they default to `.`) and rely on the flags.
//!
//! Examples:
//!   ohc_ingest --layer 15_20 --years 2016 --months 8
//!   ohc_ingest config.toml --layer 300_700 --years 2004:2025 --months 1:12
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
    provenance_link: Option<String>,
    layer: Option<LayerSpec>,
    years: Option<[i32; 2]>,
    months: Option<Vec<u32>>,
    no_ensemble: bool,
    dir_mean: Option<PathBuf>,
    dir_ensemble: Option<PathBuf>,
    dir_out: Option<PathBuf>,
}

fn parse_cli() -> Result<Cli> {
    let args: Vec<String> = env::args().skip(1).collect();
    let mut cli = Cli {
        config_path: None, tag: None, provenance_link: None, layer: None, years: None, months: None,
        no_ensemble: env::var_os("OHC_NO_ENSEMBLE").is_some(),
        dir_mean: None, dir_ensemble: None, dir_out: None,
    };
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--no-ensemble" => {
                cli.no_ensemble = true;
            }
            "--dir_mean" => {
                i += 1;
                cli.dir_mean = Some(PathBuf::from(args.get(i).context("--dir_mean needs a value")?.clone()));
            }
            "--dir_ensemble" => {
                i += 1;
                cli.dir_ensemble = Some(PathBuf::from(args.get(i).context("--dir_ensemble needs a value")?.clone()));
            }
            "--dir_out" => {
                i += 1;
                cli.dir_out = Some(PathBuf::from(args.get(i).context("--dir_out needs a value")?.clone()));
            }
            "--tag" => {
                i += 1;
                cli.tag = Some(args.get(i).context("--tag needs a value")?.clone());
            }
            "--provenance-link" => {
                i += 1;
                cli.provenance_link =
                    Some(args.get(i).context("--provenance-link needs a value")?.clone());
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
    // CLI path overrides win over config.toml / env (handy for munging dirs per shell submission).
    if let Some(p) = &cli.dir_mean { cfg.dir_mean = p.clone(); }
    if let Some(p) = &cli.dir_ensemble { cfg.dir_ensemble = p.clone(); }
    if let Some(p) = &cli.dir_out { cfg.dir_out = p.clone(); }
    cfg.run_tag = match &cli.tag {
        Some(t) => t.clone(),
        None => match env::var("OHC_TAG") {
            Ok(v) => v,
            Err(_) => bail!("--tag is required (run identifier), e.g. --tag OP20260110"),
        },
    };
    // Strip whitespace; never lowercase or otherwise munge — the tag must match the provenance
    // record char-for-char (and it is the store's directory-name token).
    cfg.run_tag = cfg.run_tag.split_whitespace().collect();
    cfg.provenance_link = match &cli.provenance_link {
        Some(l) => l.clone(),
        None => match env::var("OHC_PROVENANCE_LINK") {
            Ok(v) => v,
            Err(_) => bail!("--provenance-link is required (pointer to the provenance record)"),
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
    eprintln!(">>> layer {} dbar{}", slice.layer.tag(),
        if cli.no_ensemble { " (mean-only, --no-ensemble)" } else { "" });
    let data = ingest::ingest_layer(&cfg, &slice, &grid, cli.no_ensemble)
        .with_context(|| format!("ingesting layer {}", slice.layer.tag()))?;
    let (never, incomplete) = masks::compute_validity(&data.ohc_mean);
    // Union the member NaN footprints (matches the original's mean∪members mask); None if mean-only.
    let ens_incomplete = data.ohc_ensemble.as_ref().map(masks::compute_ensemble_incomplete);
    let flags = masks::build_flags(
        &grid, &slice.layer, &etopo, Some(&basin_id),
        cfg.latitude_range_to_keep, &cfg.basins_to_remove, cfg.bathy_floor_m, &never, &incomplete,
        ens_incomplete.as_ref(),
    )?;
    let store = zarrwrite::write_layer_store(
        &cfg, &slice, &grid, &data, &flags, &etopo, &basin_id, &cell_area,
    )?;
    eprintln!("wrote {} ({:.1}s)", store.display(), t0.elapsed().as_secs_f64());
    Ok(())
}
