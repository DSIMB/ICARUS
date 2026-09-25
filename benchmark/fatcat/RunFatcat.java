import java.io.*;
import java.util.*;
import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.io.PDBFileReader;
import org.biojava.nbio.structure.align.AFPTwister;
import org.biojava.nbio.structure.align.fatcat.FatCatFlexible;
import org.biojava.nbio.structure.align.fatcat.FatCatRigid;
import org.biojava.nbio.structure.align.model.AFPChain;

/**
 * Batch jFATCAT driver for the ICARUS benchmark.
 * Usage: RunFatcat flexible|rigid pairs.tsv struct_dir out_dir
 * For each pair "query target", superposes the query onto the (fixed) target and
 * writes out_dir/query__target.pdb with the (twisted) query C-alpha atoms, plus
 * a timing table.
 */
public class RunFatcat {
    public static void main(String[] args) throws Exception {
        boolean flex = args[0].equals("flexible");
        File dir = new File(args[2]);
        File out = new File(args[3]);
        out.mkdirs();
        PDBFileReader reader = new PDBFileReader();
        PrintWriter times = new PrintWriter(new FileWriter(new File(out, "times.tsv")));
        BufferedReader br = new BufferedReader(new FileReader(args[1]));
        String line;
        while ((line = br.readLine()) != null) {
            String[] f = line.trim().split("\\s+");
            if (f.length < 2) continue;
            String q = f[0], t = f[1];
            try {
                Structure sq = reader.getStructure(new File(dir, q + ".pdb"));
                Structure st = reader.getStructure(new File(dir, t + ".pdb"));
                Atom[] caT = StructureTools.getRepresentativeAtomArray(st.getChainByIndex(0));
                Atom[] caQ = StructureTools.getRepresentativeAtomArray(sq.getChainByIndex(0));
                long t0 = System.nanoTime();
                AFPChain afp = flex ? new FatCatFlexible().align(caT, caQ) : new FatCatRigid().align(caT, caQ);
                Group[] twisted = AFPTwister.twistOptimized(afp, caT, StructureTools.cloneAtomArray(caQ));
                double secs = (System.nanoTime() - t0) / 1e9;
                PrintWriter pw = new PrintWriter(new FileWriter(new File(out, q + "__" + t + ".pdb")));
                int n = 1;
                for (Group g : twisted) {
                    Atom a = g.getAtom("CA");
                    if (a == null) continue;
                    ResidueNumber rn = g.getResidueNumber();
                    char ic = rn.getInsCode() == null ? ' ' : rn.getInsCode();
                    pw.printf(Locale.US, "ATOM  %5d  CA  %3s A%4d%c   %8.3f%8.3f%8.3f  1.00  0.00           C%n",
                            n++, g.getPDBName(), rn.getSeqNum(), ic, a.getX(), a.getY(), a.getZ());
                }
                pw.println("END");
                pw.close();
                times.printf(Locale.US, "%s\t%s\t%.3f\t%d%n", q, t, secs, afp.getBlockNum());
                times.flush();
            } catch (Exception e) {
                System.err.println("FAILED " + q + " " + t + ": " + e);
            }
        }
        times.close();
    }
}
