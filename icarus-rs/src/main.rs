use std::collections::HashMap;
use std::io::{BufRead, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

use anyhow::{Context, Result};
use clap::{Args, Parser, Subcommand};
use rayon::prelude::*;

use icarus::align::{align_pair, AlignParams};
use icarus::db::{read_db, DbWriter};
use icarus::dp::DpWork;
use icarus::output;
use icarus::prep::{prepare, PrepParams, Prepared};
use icarus::structure::read_structure;

#[derive(Parser)]
#[command(
    name = "icarus",
    version,
    about = "ICARUS: fast flexible protein structural alignment based on Protein Units"
)]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Args, Clone)]
struct AlignOpts {
    /// Maximum number of rigid bodies (Protein Units) in a flexible solution
    #[arg(long, default_value_t = 6)]
    max_bodies: usize,
    /// Score penalty per additional rigid body (TM-score units)
    #[arg(long, default_value_t = 0.0)]
    hinge_penalty: f64,
    /// Minimum Protein Unit size (residues)
    #[arg(long, default_value_t = 15)]
    min_pu_size: usize,
    /// Maximum number of finest-level Protein Units considered
    #[arg(long, default_value_t = 10)]
    max_pus: usize,
    /// Only peel the first structure (default: peel both, keep the best)
    #[arg(long)]
    one_direction: bool,
    /// Maximum number of seed superpositions per direction
    #[arg(long, default_value_t = 1500)]
    max_seeds: usize,
    /// Candidate placements refined per Protein Unit
    #[arg(long, default_value_t = 3)]
    per_pu: usize,
    /// Fragment dRMS tolerance for seeds (Å)
    #[arg(long, default_value_t = 1.0)]
    frag_tol: f32,
    /// Stride over query fragments for seeds
    #[arg(long, default_value_t = 2)]
    seed_stride: usize,
    /// Skip the refit stage (re-placing rigid bodies on the free target)
    #[arg(long)]
    no_refit: bool,
    /// Faster preset for large-scale runs: peel one structure only, 800 seeds,
    /// 2 candidate placements per PU (~2.5x faster, ~0.02 lower TM on RIPC)
    #[arg(long)]
    fast: bool,
    /// Drop residues whose C-alpha B-factor (pLDDT for AlphaFold models) is
    /// below this value before any processing
    #[arg(long)]
    min_plddt: Option<f32>,
}

impl AlignOpts {
    fn params(&self) -> (PrepParams, AlignParams) {
        (
            PrepParams {
                min_pu_size: self.min_pu_size,
                max_leaves: self.max_pus.clamp(1, 16),
            },
            AlignParams {
                max_segments: self.max_bodies.max(1),
                hinge_penalty: self.hinge_penalty,
                q_stride: self.seed_stride.max(1),
                max_seeds: if self.fast {
                    self.max_seeds.clamp(1, 800)
                } else {
                    self.max_seeds.max(1)
                },
                per_node: if self.fast {
                    self.per_pu.clamp(1, 2)
                } else {
                    self.per_pu.max(1)
                },
                frag_tol: self.frag_tol,
                both_directions: !(self.one_direction || self.fast),
                refit: !self.no_refit,
            },
        )
    }
}

#[derive(Subcommand)]
enum Cmd {
    /// Align two structures and print a report
    Align {
        structure1: PathBuf,
        structure2: PathBuf,
        #[arg(long)]
        chain1: Option<String>,
        #[arg(long)]
        chain2: Option<String>,
        /// Write the flexibly superposed (peeled) structure to this PDB file
        #[arg(long)]
        out_pdb: Option<PathBuf>,
        /// In --out-pdb, keep sequence order instead of target (chimera) order
        #[arg(long)]
        sequence_order: bool,
        /// Print a single TSV record instead of the report
        #[arg(long)]
        tsv: bool,
        #[command(flatten)]
        opts: AlignOpts,
    },
    /// Align a list of pairs (TSV: id1 id2 per line); structures are read and
    /// preprocessed once, pairs are aligned in parallel
    Pairs {
        pairs: PathBuf,
        /// Directory holding the structures (<id><ext>)
        #[arg(long, default_value = ".")]
        dir: PathBuf,
        /// File name suffix appended to ids
        #[arg(long, default_value = "")]
        ext: String,
        /// Output TSV (default: stdout)
        #[arg(short, long)]
        output: Option<PathBuf>,
        /// Also write each superposed model into this directory
        #[arg(long)]
        models: Option<PathBuf>,
        /// Threads (0 = all cores)
        #[arg(short, long, default_value_t = 0)]
        threads: usize,
        #[command(flatten)]
        opts: AlignOpts,
    },
    /// Score two superposed structures as gdt2.pl does (no superposition):
    /// sequential DP alignment, TM-score normalised by the shortest chain
    Gdt {
        structure1: PathBuf,
        structure2: PathBuf,
        /// Normalisation length (default: shortest chain)
        #[arg(long)]
        len: Option<usize>,
        /// Print aligned residue pairs
        #[arg(long)]
        pairs: bool,
    },
    /// Preprocess structures (parse, DSSP, Protein Peeling) into a database
    Createdb {
        /// Directory (searched recursively) or text file listing structure paths
        input: PathBuf,
        /// Output database (.icdb)
        output: PathBuf,
        /// Threads (0 = all cores)
        #[arg(short, long, default_value_t = 0)]
        threads: usize,
        #[command(flatten)]
        opts: AlignOpts,
    },
    /// Flexible alignment of database structures: all-vs-all, or the candidate
    /// pairs of a list (e.g. a Foldseek .m8 prefilter result)
    Search {
        query_db: PathBuf,
        target_db: PathBuf,
        /// Candidate pairs (first two columns: query and target names)
        #[arg(long)]
        pairs: Option<PathBuf>,
        /// Output TSV (default: stdout)
        #[arg(short, long)]
        output: Option<PathBuf>,
        /// Only report pairs with a flexible TM-score >= this value
        #[arg(long, default_value_t = 0.0)]
        min_tm: f64,
        /// Only report pairs with a connectivity-aware TM-score (tm_conn,
        /// normalised by the longer chain) >= this value
        #[arg(long, default_value_t = 0.0)]
        min_conn: f64,
        /// Only report pairs with a rigid TM-score normalised by the longer
        /// chain (tm_rigid_max, the homology score) >= this value
        #[arg(long, default_value_t = 0.0)]
        min_rigid: f64,
        /// Add a column with the superposition of every rigid body
        #[arg(long)]
        transforms: bool,
        /// Threads (0 = all cores)
        #[arg(short, long, default_value_t = 0)]
        threads: usize,
        #[command(flatten)]
        opts: AlignOpts,
    },
    /// Print the Protein Unit hierarchy of a structure
    Peel {
        structure: PathBuf,
        #[arg(long)]
        chain: Option<String>,
        #[arg(long, default_value_t = 15)]
        min_pu_size: usize,
        #[arg(long, default_value_t = 10)]
        max_pus: usize,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.cmd {
        Cmd::Align {
            structure1,
            structure2,
            chain1,
            chain2,
            out_pdb,
            sequence_order,
            tsv,
            opts,
        } => {
            let (pp, ap) = opts.params();
            let keep = out_pdb.is_some();
            let t0 = Instant::now();
            let a = prepare(load(&structure1, chain1.as_deref(), keep, &opts)?, &pp);
            let b = prepare(load(&structure2, chain2.as_deref(), keep, &opts)?, &pp);
            let t1 = Instant::now();
            let mut work = DpWork::default();
            let r = align_pair(&a, &b, &ap, &mut work);
            let t2 = Instant::now();
            if tsv {
                println!("{}", output::TSV_HEADER);
                println!("{}", output::tsv_line(&r, &a, &b));
            } else {
                print!("{}", output::report(&r, &a, &b));
                println!(
                    "  Time: preprocessing {:.1} ms, alignment {:.1} ms",
                    (t1 - t0).as_secs_f64() * 1e3,
                    (t2 - t1).as_secs_f64() * 1e3
                );
            }
            if let Some(path) = out_pdb {
                let mv = if r.reversed { &b } else { &a };
                let mut w = BufWriter::new(std::fs::File::create(&path)?);
                output::write_moved_pdb(&mut w, &r.flex, mv, !sequence_order)?;
            }
        }
        Cmd::Pairs {
            pairs,
            dir,
            ext,
            output: out,
            models,
            threads,
            opts,
        } => {
            if threads > 0 {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build_global()
                    .ok();
            }
            let (pp, ap) = opts.params();
            let list = read_pairs(&pairs)?;
            let mut ids: Vec<String> = list
                .iter()
                .flat_map(|(a, b)| [a.clone(), b.clone()])
                .collect();
            ids.sort();
            ids.dedup();
            let keep = models.is_some();
            let t0 = Instant::now();
            let prepared: HashMap<String, Arc<Prepared>> = ids
                .par_iter()
                .filter_map(|id| {
                    let path = dir.join(format!("{id}{ext}"));
                    match load(&path, None, keep, &opts) {
                        Ok(s) => Some((id.clone(), Arc::new(prepare(s, &pp)))),
                        Err(e) => {
                            eprintln!("warning: {id}: {e:#}");
                            None
                        }
                    }
                })
                .collect();
            let t1 = Instant::now();
            if let Some(m) = &models {
                std::fs::create_dir_all(m)?;
            }
            let lines: Vec<String> = list
                .par_iter()
                .map_init(DpWork::default, |work, (ia, ib)| {
                    let (Some(a), Some(b)) = (prepared.get(ia), prepared.get(ib)) else {
                        return None;
                    };
                    let ts = Instant::now();
                    let r = align_pair(a, b, &ap, work);
                    let ms = ts.elapsed().as_secs_f64() * 1e3;
                    if let Some(m) = &models {
                        let mv = if r.reversed { b } else { a };
                        let path = m.join(format!("{ia}__{ib}.pdb"));
                        if let Ok(f) = std::fs::File::create(&path) {
                            let mut w = BufWriter::new(f);
                            let _ = output::write_moved_pdb(&mut w, &r.flex, mv, true);
                        }
                    }
                    Some(format!("{}\t{:.2}", output::tsv_line(&r, a, b), ms))
                })
                .flatten()
                .collect();
            let t2 = Instant::now();
            let mut w: Box<dyn Write> = match &out {
                Some(p) => Box::new(BufWriter::new(std::fs::File::create(p)?)),
                None => Box::new(BufWriter::new(std::io::stdout())),
            };
            writeln!(w, "{}\tms", output::TSV_HEADER)?;
            for l in &lines {
                writeln!(w, "{l}")?;
            }
            eprintln!(
                "{} structures prepared in {:.2} s, {} pairs aligned in {:.2} s ({:.1} pairs/s)",
                prepared.len(),
                (t1 - t0).as_secs_f64(),
                lines.len(),
                (t2 - t1).as_secs_f64(),
                lines.len() as f64 / (t2 - t1).as_secs_f64().max(1e-9)
            );
        }
        Cmd::Createdb {
            input,
            output,
            threads,
            opts,
        } => {
            if threads > 0 {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build_global()
                    .ok();
            }
            let (pp, _) = opts.params();
            let files = collect_files(&input)?;
            let t0 = Instant::now();
            // chunks bound memory on large collections (e.g. all of Swiss-Prot)
            let mut db = DbWriter::create(&output)?;
            let mut nres = 0usize;
            for (ci, chunk) in files.chunks(20_000).enumerate() {
                let items: Vec<Prepared> = chunk
                    .par_iter()
                    .filter_map(|p| match load(p, None, false, &opts) {
                        Ok(s) => Some(prepare(s, &pp)),
                        Err(e) => {
                            eprintln!("warning: {}: {e:#}", p.display());
                            None
                        }
                    })
                    .collect();
                for p in &items {
                    nres += p.len();
                    db.push(p)?;
                }
                if files.len() > 20_000 {
                    eprintln!(
                        "  {}/{} files ({:.0} s)",
                        (ci * 20_000 + chunk.len()),
                        files.len(),
                        t0.elapsed().as_secs_f64()
                    );
                }
            }
            let n = db.finish()?;
            eprintln!(
                "{} structures ({} residues) preprocessed in {:.2} s -> {}",
                n,
                nres,
                t0.elapsed().as_secs_f64(),
                output.display()
            );
        }
        Cmd::Search {
            query_db,
            target_db,
            pairs,
            output: out,
            min_tm,
            min_conn,
            min_rigid,
            transforms,
            threads,
            opts,
        } => {
            if threads > 0 {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build_global()
                    .ok();
            }
            let (_, ap) = opts.params();
            let t0 = Instant::now();
            let same = query_db == target_db;
            let qs = read_db(&query_db)?;
            let ts = if same {
                Vec::new()
            } else {
                read_db(&target_db)?
            };
            let tset: &Vec<Prepared> = if same { &qs } else { &ts };
            let qidx: HashMap<&str, usize> = qs
                .iter()
                .enumerate()
                .map(|(i, p)| (p.s.name.as_str(), i))
                .collect();
            let tidx: HashMap<&str, usize> = tset
                .iter()
                .enumerate()
                .map(|(i, p)| (p.s.name.as_str(), i))
                .collect();
            let lookup = |m: &HashMap<&str, usize>, n: &str| -> Option<usize> {
                let n = strip_name(n);
                m.get(n.as_str())
                    .copied()
                    .or_else(|| n.rsplit_once('_').and_then(|(a, _)| m.get(a).copied()))
            };
            let list: Vec<(usize, usize)> = match &pairs {
                Some(p) => {
                    let mut v: Vec<(usize, usize)> = read_pairs(p)?
                        .iter()
                        .filter_map(|(a, b)| Some((lookup(&qidx, a)?, lookup(&tidx, b)?)))
                        .filter(|(a, b)| !(same && a == b))
                        .collect();
                    if same {
                        v.iter_mut().for_each(|x| *x = (x.0.min(x.1), x.0.max(x.1)));
                    }
                    v.sort_unstable();
                    v.dedup();
                    v
                }
                None if same => (0..qs.len())
                    .flat_map(|i| (i + 1..qs.len()).map(move |j| (i, j)))
                    .collect(),
                None => (0..qs.len())
                    .flat_map(|i| (0..tset.len()).map(move |j| (i, j)))
                    .collect(),
            };
            let t1 = Instant::now();
            eprintln!(
                "loaded {} + {} structures in {:.2} s; {} pairs to align",
                qs.len(),
                tset.len(),
                (t1 - t0).as_secs_f64(),
                list.len()
            );
            let mut w: Box<dyn Write> = match &out {
                Some(p) => Box::new(BufWriter::new(std::fs::File::create(p)?)),
                None => Box::new(BufWriter::new(std::io::stdout())),
            };
            writeln!(
                w,
                "{}\tms{}",
                output::TSV_HEADER,
                if transforms { "\ttransforms" } else { "" }
            )?;
            let mut done = 0usize;
            let mut kept = 0usize;
            for chunk in list.chunks(20_000) {
                let lines: Vec<String> = chunk
                    .par_iter()
                    .map_init(DpWork::default, |work, &(i, j)| {
                        let (a, b) = (&qs[i], &tset[j]);
                        let ts = Instant::now();
                        let r = align_pair(a, b, &ap, work);
                        if r.tm_flex() < min_tm
                            || r.tm_conn() < min_conn
                            || (min_rigid > 0.0 && output::tm_rigid_max(&r, a, b) < min_rigid)
                        {
                            return None;
                        }
                        Some(format!(
                            "{}\t{:.2}",
                            output::tsv_line(&r, a, b),
                            ts.elapsed().as_secs_f64() * 1e3
                        ))
                    })
                    .flatten()
                    .collect();
                for l in &lines {
                    writeln!(w, "{l}")?;
                }
                done += chunk.len();
                kept += lines.len();
                let el = t1.elapsed().as_secs_f64();
                eprintln!(
                    "  {done}/{} pairs ({:.1} pairs/s), {kept} reported",
                    list.len(),
                    done as f64 / el.max(1e-9)
                );
            }
            w.flush()?;
        }
        Cmd::Gdt {
            structure1,
            structure2,
            len,
            pairs,
        } => {
            let a = read_structure(&structure1, None, false)?;
            let b = read_structure(&structure2, None, false)?;
            let g = output::gdt(&a.ca, &b.ca, len);
            println!("{:.4}\t{:.4}\t{}", g.tm, g.tm_search, g.n_aligned);
            if pairs {
                for (i, j, d) in g.pairs {
                    println!("{}\t{}\t{:.3}", a.resid[i as usize], b.resid[j as usize], d);
                }
            }
        }
        Cmd::Peel {
            structure,
            chain,
            min_pu_size,
            max_pus,
        } => {
            let s = read_structure(&structure, chain.as_deref(), false)?;
            let p = prepare(
                s,
                &PrepParams {
                    min_pu_size,
                    max_leaves: max_pus,
                },
            );
            println!("{} residues; secondary structure:", p.len());
            println!(
                "{}",
                p.ss.iter()
                    .map(|&c| if c == ' ' { '-' } else { c })
                    .collect::<String>()
            );
            for (lvl, nodes) in p.tree.levels.iter().enumerate() {
                let spans: Vec<String> = nodes
                    .iter()
                    .map(|&n| {
                        format!(
                            "{}-{}",
                            p.s.resid[p.tree.nodes[n].start], p.s.resid[p.tree.nodes[n].end]
                        )
                    })
                    .collect();
                println!("level {lvl}: {} PUs  {}", nodes.len(), spans.join(" "));
            }
        }
    }
    Ok(())
}

fn load(
    path: &Path,
    chain: Option<&str>,
    keep: bool,
    opts: &AlignOpts,
) -> Result<icarus::structure::Structure> {
    let mut s = read_structure(path, chain, keep)?;
    if let Some(t) = opts.min_plddt {
        s.mask_low_confidence(t);
        anyhow::ensure!(
            s.len() >= 10,
            "fewer than 10 residues left after pLDDT masking"
        );
    }
    Ok(s)
}

fn is_structure_file(p: &Path) -> bool {
    let s = p.to_string_lossy().to_ascii_lowercase();
    let s = s.strip_suffix(".gz").unwrap_or(&s);
    if [".pdb", ".cif", ".ent", ".mmcif"]
        .iter()
        .any(|e| s.ends_with(e))
    {
        return true;
    }
    // other names (e.g. SCOP/ASTRAL domain files "d1abca_", "d1apy.1") are
    // accepted unless they carry an obviously non-structural extension
    let deny = [
        ".txt", ".tsv", ".csv", ".json", ".md", ".log", ".m8", ".icdb", ".fasta", ".fa", ".py",
        ".sh", ".tar", ".zip", ".dbtype", ".index", ".lookup", ".source",
    ];
    !deny.iter().any(|e| s.ends_with(e))
}

fn collect_files(input: &Path) -> Result<Vec<PathBuf>> {
    let mut out = Vec::new();
    if input.is_dir() {
        let mut stack = vec![input.to_path_buf()];
        while let Some(d) = stack.pop() {
            for e in std::fs::read_dir(&d)? {
                let p = e?.path();
                if p.is_dir() {
                    stack.push(p);
                } else if is_structure_file(&p) {
                    out.push(p);
                }
            }
        }
    } else {
        for line in std::io::BufReader::new(std::fs::File::open(input)?).lines() {
            let l = line?;
            let l = l.trim();
            if !l.is_empty() && !l.starts_with('#') {
                out.push(PathBuf::from(l));
            }
        }
    }
    out.sort();
    // one file per structure name (e.g. AFDB ships both .cif.gz and .pdb.gz)
    let mut seen = std::collections::HashSet::new();
    out.retain(|p| {
        seen.insert(strip_name(
            &p.file_name().unwrap_or_default().to_string_lossy(),
        ))
    });
    Ok(out)
}

/// Structure name as stored in databases: file name without structure extensions.
fn strip_name(s: &str) -> String {
    let mut n = s.to_string();
    for ext in [".gz", ".pdb", ".cif", ".ent", ".mmcif"] {
        if let Some(x) = n.strip_suffix(ext) {
            n = x.to_string();
        }
    }
    n
}

fn read_pairs(path: &Path) -> Result<Vec<(String, String)>> {
    let f = std::fs::File::open(path).with_context(|| format!("cannot open {}", path.display()))?;
    let mut out = Vec::new();
    for line in std::io::BufReader::new(f).lines() {
        let line = line?;
        let t: Vec<&str> = line.split_whitespace().collect();
        if t.len() >= 2 && !t[0].starts_with('#') {
            out.push((t[0].to_string(), t[1].to_string()));
        }
    }
    Ok(out)
}
