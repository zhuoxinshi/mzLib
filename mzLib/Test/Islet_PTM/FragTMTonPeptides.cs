using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using NUnit.Framework;

namespace Test.Islet_PTM
{
    // =============================================================================
    // C# port of PsmRatio_PeptideQuant.Rmd
    // "Ratio-based PSM-to-peptide quantification (TMT-Integrator style, MM PSMs)".
    //
    // Re-implements the TMT-Integrator ratio-based quantification (Chang et al.,
    // Nat Commun 2026) on MetaMorpheus PSM tables. All PSMs from the three TMT18
    // plexes (plates R2-R3, R4-R5, R6-R7) are combined into one integrated report
    // in ratio space, exactly as FragPipe does.
    //
    // Per-plex pipeline (each plex = 9 fraction *_PSMs.psmtsv files):
    //   1. FILTER      drop decoy/contaminant, require QValue <= 0.01, require the
    //                  TMT18 label, drop the lowest 5% of PSMs by summed reporter
    //                  intensity, keep the best (max-summed) PSM per (peptide,file).
    //   2. RATIO       log2 each reporter channel and subtract the "virtual
    //                  reference" = the MEAN of that PSM's own reporter intensities.
    //   3. GROUP       group PSMs by peptide (Full Sequence).
    //   4. OUTLIERS    per (peptide, channel) drop IQR-fence outliers when >= 4 vals.
    //   5. AGGREGATE   median of the surviving PSM ratios per channel.
    //   6. NORM II     optional per-sample median centering (the "MD" table).
    //   7. ABUNDANCE   anchor each peptide to Ref_i (weighted sum of top-3 PSM MS1
    //                  precursor intensities, median across plexes); log2 abundance
    //                  = ratio + log2(Ref_i).
    //
    // R -> C# cheat-sheet for the data.table calls in the Rmd:
    //   fread(select=..)          -> ParsePsmFile: read a TSV, keep only needed cols
    //   rbindlist(lapply(fs,..))  -> concatenate rows over the 9 fraction files
    //   chm[chm<=0] <- NA         -> intensities <= 0 become double.NaN (missing)
    //   rowSums(x, na.rm=TRUE)    -> Sum ignoring NaN
    //   quantile(x, p, type=7)    -> Quantile7 (R's default interpolation)
    //   which.max by group        -> best-PSM pick per (FullSeq, File)
    //   sweep(log2(chm),1,..,'-') -> per-row subtract of log2(reference)
    //   melt + dcast + agg_cell   -> group by peptide, AggCell per channel
    //   merge(all=TRUE) / Reduce  -> full outer join of the 3 plex tables on FullSeq
    //   apply(.,1,median,na.rm)   -> Ref_i = row-wise median of the 3 Ref columns
    //   sweep(R,2,colmedian,'-')  -> MD: subtract each column's median
    //
    // Not ported: the ggplot2 boxplot/PCA QC panels (sections 5-6) are plotting-only
    // and produce no data file. This port reproduces the three data outputs:
    //   TMTI_abundance_peptide_None.tsv, _MD.tsv, and TMTI_ratio_peptide_None.tsv.
    // =============================================================================
    internal class FragTMTonPeptides
    {
        // ---- 1. Paths and parameters (mirrors the "params" chunk) -----------------
        private const string InDir =
            @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\UPLC\MM\LFgptmdPrunedDb_search-cali-search\Task1-SearchTask\Individual File Results";
        private const string OutDir =
            @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\DataAnalysis";

        // plate id -> the token that appears in its file names ("R2_R3" -> "TMT_R2-R3")
        private static readonly (string Plate, string FileToken)[] Plates =
        {
            ("R2_R3", "TMT_R2-R3"), ("R4_R5", "TMT_R4-R5"), ("R6_R7", "TMT_R6-R7")
        };

        // 18 TMT18 reporter channels, in TSV column order.
        private static readonly string[] Channels =
        {
            "126","127N","127C","128N","128C","129N","129C","130N","130C",
            "131N","131C","132N","132C","133N","133C","134N","134C","135N"
        };

        private const double QVAL_MAX = 0.01;   // 1% PSM-level FDR
        private const double LOW_FRAC = 0.05;   // drop lowest 5% by summed reporter intensity
        private const double IQR_MULT = 1.5;    // IQR outlier fence
        private const int    MIN_GRP  = 4;      // min ratios in a channel to apply IQR removal

        // <plate>_<channel>, with the two reference channels renamed to Pool1/Pool2.
        private static string SampleCol(string plate, string ch) =>
            ch == "134C" ? plate + "_Pool1" : ch == "135N" ? plate + "_Pool2" : plate + "_" + ch;

        // One filtered PSM row (only the columns the algorithm uses).
        private sealed class Psm
        {
            public string FullSeq, File, BaseSeq, Accession;
            public double PrecInt;
            public double[] Ch;       // 18 reporter intensities, NaN = missing
            public double Summ;       // sum of Ch ignoring NaN
            public int NnonNa;        // count of non-NaN channels
        }

        // Per-plex aggregated result: peptide -> 18 channel ratios, plus Ref and nPSM.
        private sealed class PlexResult
        {
            public readonly Dictionary<string, double[]> Ratio = new();  // 18 medians (NaN ok)
            public readonly Dictionary<string, double> Ref = new();
            public readonly Dictionary<string, int> NPsm = new();
            public readonly Dictionary<string, (string BaseSeq, string Accession, string Orf)> Meta = new();
        }

        [Test]
        public static void FragTMTonPeptidesTest()
        {
            if (!Directory.Exists(InDir))
                Assert.Ignore($"Input directory not found (local-data test): {InDir}");

            // ---- 2. Process each plex independently ------------------------------
            var results = Plates.ToDictionary(p => p.Plate, p => ProcessPlate(p.Plate, p.FileToken));

            // ---- 3. Integrate plexes: full outer join on FullSeq -----------------
            // Union of every peptide seen in any plex (data.table merge sorts by key,
            // so we sort ordinally to keep output row order comparable to R).
            var peptides = results.Values.SelectMany(r => r.Ratio.Keys)
                                  .Distinct().OrderBy(s => s, StringComparer.Ordinal).ToList();

            // meta: first occurrence across plexes in plate order (unique(.., by=FullSeq)).
            var meta = new Dictionary<string, (string BaseSeq, string Accession, string Orf)>();
            foreach (var (plate, _) in Plates)
                foreach (var kv in results[plate].Meta)
                    if (!meta.ContainsKey(kv.Key)) meta[kv.Key] = kv.Value;

            // The 54 integrated sample columns, in plate-major / channel order.
            var sampCols = new List<string>();
            foreach (var (plate, _) in Plates)
                foreach (var ch in Channels) sampCols.Add(SampleCol(plate, ch));

            int nPep = peptides.Count, nCol = sampCols.Count;
            var R      = new double[nPep, nCol];    // integrated log2 ratios (None)
            var refI   = new double[nPep];          // Ref_i (median of the 3 Ref columns)
            var numPsm = new int[nPep];             // total PSMs across plexes

            for (int i = 0; i < nPep; i++)
            {
                string pep = peptides[i];
                var refVals = new List<double>(3);
                int col = 0, psmSum = 0;
                foreach (var (plate, _) in Plates)
                {
                    var res = results[plate];
                    bool present = res.Ratio.TryGetValue(pep, out var ratios);
                    for (int c = 0; c < Channels.Length; c++, col++)
                        R[i, col] = present ? ratios[c] : double.NaN;
                    if (present)
                    {
                        refVals.Add(res.Ref[pep]);
                        psmSum += res.NPsm[pep];
                    }
                }
                refI[i]   = Median(refVals);          // median across plexes, NaN-free input
                numPsm[i] = psmSum;
            }

            // ---- MD normalization: subtract each sample column's median ----------
            var colMedian = new double[nCol];
            for (int c = 0; c < nCol; c++)
            {
                var vals = new List<double>(nPep);
                for (int i = 0; i < nPep; i++) if (!double.IsNaN(R[i, c])) vals.Add(R[i, c]);
                colMedian[c] = vals.Count > 0 ? Median(vals) : 0.0;
            }

            // ---- Ratio -> abundance: A = R (+/- MD) + log2(Ref_i) -----------------
            var aNone = new double[nPep, nCol];
            var aMd   = new double[nPep, nCol];
            for (int i = 0; i < nPep; i++)
            {
                double log2ref = Math.Log2(refI[i]);
                for (int c = 0; c < nCol; c++)
                {
                    double r = R[i, c];
                    aNone[i, c] = r + log2ref;
                    aMd[i, c]   = (r - colMedian[c]) + log2ref;
                }
            }
            Console.WriteLine($"Integrated: {nPep} peptides x {nCol} channels");

            // ---- 4. Write the three data outputs ---------------------------------
            Directory.CreateDirectory(OutDir);
            WriteTable(Path.Combine(OutDir, "TMTI_abundance_peptide_None.tsv"), peptides, meta, numPsm, sampCols, aNone);
            WriteTable(Path.Combine(OutDir, "TMTI_abundance_peptide_MD.tsv"),   peptides, meta, numPsm, sampCols, aMd);
            WriteTable(Path.Combine(OutDir, "TMTI_ratio_peptide_None.tsv"),     peptides, meta, numPsm, sampCols, R);
            Console.WriteLine("Wrote 2 abundance tables (None, MD) + 1 ratio-only table");
        }

        // -------------------------------------------------------------------------
        // Per-plex processing: filter -> ratio -> group -> IQR + median per channel.
        // -------------------------------------------------------------------------
        private static PlexResult ProcessPlate(string plate, string fileToken)
        {
            var files = Directory.GetFiles(InDir, "*_PSMs.psmtsv")
                                 .Where(f => Path.GetFileName(f).Contains(fileToken))
                                 .OrderBy(f => f, StringComparer.Ordinal).ToArray();
            if (files.Length != 9)
                throw new Exception($"Expected 9 PSM files for {plate}, found {files.Length}");

            // rbindlist: read all 9 fractions, then apply Step-1 filters row by row.
            var psms = new List<Psm>();
            foreach (var f in files) ParsePsmFile(f, psms);

            // Drop the lowest 5% of PSMs by summed reporter intensity (strict >).
            var summs = psms.Select(p => p.Summ).OrderBy(x => x).ToArray();
            double thr = Quantile7(summs, LOW_FRAC);
            psms = psms.Where(p => p.Summ > thr).ToList();

            // Best-PSM: keep the max-summed PSM per (FullSeq, File); first wins on ties.
            var best = new Dictionary<(string, string), Psm>();
            foreach (var p in psms)
            {
                var key = (p.FullSeq, p.File);
                if (!best.TryGetValue(key, out var cur) || p.Summ > cur.Summ) best[key] = p;
            }
            psms = best.Values.ToList();
            Console.WriteLine($"  {plate}: {psms.Count} PSMs after filtering");

            // Group PSMs by peptide (Full Sequence).
            var byPeptide = new Dictionary<string, List<Psm>>();
            foreach (var p in psms)
            {
                if (!byPeptide.TryGetValue(p.FullSeq, out var lst)) byPeptide[p.FullSeq] = lst = new List<Psm>();
                lst.Add(p);
            }

            var res = new PlexResult();
            foreach (var kv in byPeptide)
            {
                string pep = kv.Key;
                var group = kv.Value;

                // Steps 2-5: per channel, gather each PSM's virtual-reference log2 ratio
                //   ratio = log2(channel) - log2(mean of that PSM's reporters),
                // then IQR-trim and take the median across the peptide's PSMs.
                var med = new double[Channels.Length];
                for (int c = 0; c < Channels.Length; c++)
                {
                    var vals = new List<double>(group.Count);
                    foreach (var p in group)
                    {
                        double v = p.Ch[c];
                        if (double.IsNaN(v)) continue;                       // missing channel
                        double reference = p.Summ / p.NnonNa;               // mean = virtual ref
                        vals.Add(Math.Log2(v) - Math.Log2(reference));
                    }
                    med[c] = AggCell(vals);
                }
                res.Ratio[pep] = med;
                res.NPsm[pep] = group.Count;

                // Step-7 prep: Ref = weighted sum of the top-3 PSMs by precursor intensity,
                // each PSM weighted by (1 / nNonNa) * precursor intensity.
                res.Ref[pep] = group.OrderByDescending(p => p.PrecInt).Take(3)
                                    .Sum(p => (1.0 / p.NnonNa) * p.PrecInt);

                // meta from the first PSM; ORF = accession up to the first space.
                var first = group[0];
                int sp = first.Accession.IndexOf(' ');
                res.Meta[pep] = (first.BaseSeq, first.Accession,
                                 sp < 0 ? first.Accession : first.Accession.Substring(0, sp));
            }
            return res;
        }

        // -------------------------------------------------------------------------
        // Parse one *_PSMs.psmtsv, applying the row-level Step-1 filters:
        //   Decoy != Y, Contaminant != Y, QValue <= 0.01, has the TMT18 label,
        //   and at least one positive reporter channel.
        // -------------------------------------------------------------------------
        private static void ParsePsmFile(string path, List<Psm> outPsms)
        {
            using var sr = new StreamReader(path);
            var header = sr.ReadLine()?.Split('\t') ?? throw new Exception($"Empty file: {path}");
            int Col(string n) { int i = Array.IndexOf(header, n); if (i < 0) throw new Exception($"Missing '{n}' in {path}"); return i; }

            int iFile = Col("File Name"), iFull = Col("Full Sequence"), iBase = Col("Base Sequence");
            int iAcc = Col("Accession"), iCont = Col("Contaminant"), iDec = Col("Decoy");
            int iQ = Col("QValue"), iPrec = Col("Precursor Intensity");
            var iCh = Channels.Select(Col).ToArray();
            var inv = CultureInfo.InvariantCulture;

            string line;
            while ((line = sr.ReadLine()) != null)
            {
                var p = line.Split('\t');
                if (p.Length < header.Length) continue;
                if (p[iDec] == "Y" || p[iCont] == "Y") continue;
                if (!double.TryParse(p[iQ], NumberStyles.Float, inv, out double q) || q > QVAL_MAX) continue;
                if (!p[iFull].Contains("Multiplex Label:TMT18")) continue;

                var ch = new double[Channels.Length];
                double summ = 0; int nn = 0;
                for (int c = 0; c < Channels.Length; c++)
                {
                    if (double.TryParse(p[iCh[c]], NumberStyles.Float, inv, out double v) && v > 0)
                    { ch[c] = v; summ += v; nn++; }
                    else ch[c] = double.NaN;               // <= 0 / unparseable -> missing
                }
                if (summ <= 0) continue;                   // keep summ > 0

                double.TryParse(p[iPrec], NumberStyles.Float, inv, out double prec);
                outPsms.Add(new Psm
                {
                    FullSeq = p[iFull], File = p[iFile], BaseSeq = p[iBase], Accession = p[iAcc],
                    PrecInt = prec, Ch = ch, Summ = summ, NnonNa = nn
                });
            }
        }

        // agg_cell: drop IQR-fence outliers (only when >= MIN_GRP values), then median.
        private static double AggCell(List<double> x)
        {
            if (x.Count == 0) return double.NaN;
            if (x.Count >= MIN_GRP)
            {
                var sorted = x.OrderBy(v => v).ToArray();
                double q1 = Quantile7(sorted, 0.25), q3 = Quantile7(sorted, 0.75), iqr = q3 - q1;
                double lo = q1 - IQR_MULT * iqr, hi = q3 + IQR_MULT * iqr;
                x = x.Where(v => v >= lo && v <= hi).ToList();
            }
            return Median(x);
        }

        // R's type-7 quantile (the quantile()/median() default): linear interpolation
        // at position (n-1)*p on the sorted values. Input must be sorted ascending.
        private static double Quantile7(double[] sorted, double p)
        {
            int n = sorted.Length;
            if (n == 0) return double.NaN;
            if (n == 1) return sorted[0];
            double h = (n - 1) * p;
            int lo = (int)Math.Floor(h);
            if (lo + 1 >= n) return sorted[n - 1];
            return sorted[lo] + (h - lo) * (sorted[lo + 1] - sorted[lo]);
        }

        private static double Median(List<double> x)
        {
            if (x.Count == 0) return double.NaN;
            var a = x.ToArray(); Array.Sort(a);
            int n = a.Length;
            return n % 2 == 1 ? a[n / 2] : 0.5 * (a[n / 2 - 1] + a[n / 2]);
        }

        // -------------------------------------------------------------------------
        // Write one output table: meta columns + the 54 sample columns (round 4 dp,
        // NaN -> empty, matching data.table fwrite defaults).
        // -------------------------------------------------------------------------
        private static void WriteTable(string path, List<string> peptides,
            Dictionary<string, (string BaseSeq, string Accession, string Orf)> meta,
            int[] numPsm, List<string> sampCols, double[,] mat)
        {
            var inv = CultureInfo.InvariantCulture;
            using var w = new StreamWriter(path);
            w.WriteLine(string.Join('\t', new[] { "Index", "BaseSequence", "ProteinAccession", "ORF", "NumPSMs" }.Concat(sampCols)));
            for (int i = 0; i < peptides.Count; i++)
            {
                var m = meta[peptides[i]];
                var sb = new StringBuilder();
                sb.Append(peptides[i]).Append('\t').Append(m.BaseSeq).Append('\t')
                  .Append(m.Accession).Append('\t').Append(m.Orf).Append('\t').Append(numPsm[i]);
                for (int c = 0; c < sampCols.Count; c++)
                {
                    double v = mat[i, c];
                    sb.Append('\t');
                    if (!double.IsNaN(v)) sb.Append(Math.Round(v, 4, MidpointRounding.ToEven).ToString("0.####", inv));
                }
                w.WriteLine(sb.ToString());
            }
        }
    }
}
