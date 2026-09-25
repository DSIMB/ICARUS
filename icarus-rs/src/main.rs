use std::collections::HashMap;
use std::io::{BufRead, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

use anyhow::{Context, Result};
use clap::{Args, Parser, Subcommand};
use rayon::prelude::*;

use icarus::align::{align_pair, AlignParams};
use icarus::dp::DpWork;
use icarus::output;
use icarus::prep::{prepare, PrepParams, Prepared};
use icarus::structure::read_structure;

#[derive(Parser)]
#[command(name = "icarus", version, about = "ICARUS: fast flexible protein structural alignment based on Protein Units")]
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
}

impl AlignOpts {
    fn params(&self) -> (PrepParams, AlignParams) {
        (
            PrepParams { min_pu_size: self.min_pu_size, max_leaves: self.max_pus.clamp(1, 16) },
            AlignParams {
                max_segments: self.max_bodies.max(1),
                hinge_penalty: self.hinge_penalty,
                q_stride: self.seed_stride.max(1),
                max_seeds: self.max_seeds.max(1),
                per_node: self.per_pu.max(1),
                frag_tol: self.frag_tol,
                both_directions: !self.one_direction,
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
        Cmd::Align { structure1, structure2, chain1, chain2, out_pdb, sequence_order, tsv, opts } => {
            let (pp, ap) = opts.params();
            let keep = out_pdb.is_some();
            let t0 = Instant::now();
            let a = prepare(read_structure(&structure1, chain1.as_deref(), keep)?, &pp);
            let b = prepare(read_structure(&structure2, chain2.as_deref(), keep)?, &pp);
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
        Cmd::Pairs { pairs, dir, ext, output: out, models, threads, opts } => {
            if threads > 0 {
                rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().ok();
            }
            let (pp, ap) = opts.params();
            let list = read_pairs(&pairs)?;
            let mut ids: Vec<String> = list.iter().flat_map(|(a, b)| [a.clone(), b.clone()]).collect();
            ids.sort();
            ids.dedup();
            let keep = models.is_some();
            let t0 = Instant::now();
            let prepared: HashMap<String, Arc<Prepared>> = ids
                .par_iter()
                .filter_map(|id| {
                    let path = dir.join(format!("{id}{ext}"));
                    match read_structure(&path, None, keep) {
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
        Cmd::Gdt { structure1, structure2, len, pairs } => {
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
        Cmd::Peel { structure, chain, min_pu_size, max_pus } => {
            let s = read_structure(&structure, chain.as_deref(), false)?;
            let p = prepare(s, &PrepParams { min_pu_size, max_leaves: max_pus });
            println!("{} residues; secondary structure:", p.len());
            println!("{}", p.ss.iter().map(|&c| if c == ' ' { '-' } else { c }).collect::<String>());
            for (lvl, nodes) in p.tree.levels.iter().enumerate() {
                let spans: Vec<String> = nodes
                    .iter()
                    .map(|&n| format!("{}-{}", p.s.resid[p.tree.nodes[n].start], p.s.resid[p.tree.nodes[n].end]))
                    .collect();
                println!("level {lvl}: {} PUs  {}", nodes.len(), spans.join(" "));
            }
        }
    }
    Ok(())
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
