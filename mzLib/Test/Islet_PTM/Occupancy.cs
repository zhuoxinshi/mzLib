using NUnit.Framework;
using Readers;
using Readers.Generated;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace Test.Islet_PTM
{
    public record PtmRecord
    {
        public string Accession { get; set; }
        public int ModSite { get; set; }
        public string Mod { get; set; }
    }

    public class PtmSiteGroup
    {
        public PtmRecord PTM { get; set; }
        public IEnumerable<PsmFromTsv> UnmodPsms { get; set; }
        public IEnumerable<PsmFromTsv> ModPsms { get; set; }

        [Test]
        public static void Occupancy()
        {
            var notInteresting = new List<string> { "TMT18", "Fixed", "Artifact", "Variable", "Metal" };

            string allPSMs_path = @"E:\Islets\Brian_data\PTM\MS3_all_LFgptmdFilterPrunedDb\Task1-SearchTask\AllPSMs.psmtsv";
            var allPSMsTMT_file = new PsmFromTsvFile(allPSMs_path);
            var allPsms = allPSMsTMT_file.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T" && p.AmbiguityLevel== "1");
            var allPeptides_path = @"E:\Islets\Brian_data\PTM\MS3_all_LFgptmdFilterPrunedDb\Task1-SearchTask\AllPeptides.psmtsv";
            var allPeptidesTMT_file = new PsmFromTsvFile(allPeptides_path);
            var allPeptides = allPeptidesTMT_file.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T" && p.AmbiguityLevel == "1");
            var allPeptidesWithMod = allPeptides.Where(p => SpectrumMatchFromTsv.ParseModifications(p.FullSequence).Values.Any(v => !notInteresting.Any(key => v.Contains(key))));
            var uniquePsmsWithMod = allPsms.GroupBy(p => p.FullSequence).Select(g => g.First()).Where(p => SpectrumMatchFromTsv.ParseModifications(p.FullSequence).Values.Any(v => !notInteresting.Any(key => v.Contains(key))));

            var allPtms = new List<PtmRecord>();
            foreach (var peptide in uniquePsmsWithMod)
            {
                var allMods = SpectrumMatchFromTsv.ParseModifications(peptide.FullSequence).Where(kvp => !notInteresting.Any(key => kvp.Value.Contains(key)));
                int startAA = int.Parse(peptide.StartAndEndResiduesInProtein.Split()[0].Split("[")[1]);
                foreach (var mod in allMods)
                {
                    int modSite = mod.Key + startAA;
                    if (mod.Key != 0) modSite = modSite - 1;
                    if (mod.Key == 0 && startAA != 1) continue;
                    var modName = mod.Value.Split(':')[1].Split(' ')[0];
                    var parsedModName =  DeepLC.ParseModNameForDeepLC(modName);
                    allPtms.Add(new PtmRecord { Accession = peptide.ProteinAccession, ModSite = modSite, Mod = parsedModName });
                }
            }
            var distinctPtms = allPtms.Distinct().ToList();
            var allPtmGroups = allPtms.GroupBy(ptm => new {a = ptm.Accession, m = ptm.ModSite}).ToList();
        }

        // =====================================================================================
        // C# port of psmsToOccupancies_islet.ipynb, up to (and including) the CODE CELL 49
        // chunk that builds `output_df` and prints its shape.
        //
        // What the pipeline does: from a MetaMorpheus AllPSMs table with 18 TMT channels,
        // compute per-site PTM occupancy = (summed reporter intensity of PSMs carrying the mod
        // at a protein position) / (summed reporter intensity of ALL PSMs covering that
        // position), channel by channel.
        //
        // pandas/numpy -> C# cheat-sheet for the notebook idioms:
        //   pd.read_csv(sep='\t')                 -> hand-parsed TSV (channels "126".."135N")
        //   df[mask]                              -> LINQ Where on the parsed rows
        //   iloc[:, tmt1:tmt18+1]                 -> the 18 reporter columns, indices tmt1..tmt1+17
        //   (channels==0).any(axis=1)            -> drop a PSM if ANY channel is exactly 0
        //   col*target/col.sum()                 -> sample-loading normalization per channel
        //   Series.map(remove_mods)              -> strip "insignificant" mod brackets
        //   df.join / dict lookup (CELL 32)      -> remap Start/End by (Accession, Base Sequence)
        //   str.contains("|")                    -> drop position/sequence-ambiguous rows
        //   df.apply(row -> lookup, axis=1)      -> the two accumulation loops (CELL 42 / 43)
        //   site.mod_intensity / total_intensity -> SiteModification.CalculateOccupancy
        //   np.array([...]).T + dict-comprehension columns -> BuildOutput
        //
        // Not ported: plotting/stats-only cells (matplotlib histograms, limma_py eBayes,
        // volcano plot) after CELL 49; those produce no part of `output_df`.
        // =====================================================================================
        private const string FilePath = @"E:\Islets\Brian_data\PTM\Occupancy\AllPSMs.psmtsv";
        private const string RunSuffix = "AllPSMs";
        private const int NCh = 18;                             // number of TMT channels
        private const string Normalization = "median";         // "median" or "mean" (CELL 20)
        // If true, also writes the output table (CELL 51, one line past the requested chunk),
        // with a "_csharp" suffix so it never clobbers the notebook's own output.
        private const bool WriteOutputForInspection = true;

        private static readonly string[] InsignificantMods =
            { "Multiplex Label", "Common Variable", "Common Fixed", "Metal", "Artifact" };

        private static readonly Regex ModRe  = new(@"\[([^:\]]+):\s*([^\]]+)\]", RegexOptions.Compiled);
        private static readonly Regex SpanRe = new(@"^\[(\d+) to (\d+)\]", RegexOptions.Compiled);

        private sealed class Row
        {
            public string FullSeq, BaseSeq, Accession, Gene, StartEnd;
            public double[] Ch;                                 // 18 reporter intensities
        }

        private sealed class SiteModification
        {
            public string ModificationName, ParentAccession, ParentGene;
            public int ParentPosition;
            public double[] ModIntensity = new double[NCh];
            public double[] TotalIntensity = new double[NCh];
            public int NumberOfPsms;

            // occupancy = mod / total per channel; zeros if the position has no signal at all.
            public double[] CalculateOccupancy()
            {
                if (!TotalIntensity.Any(t => t > 0)) return new double[NCh];
                var occ = new double[NCh];
                for (int i = 0; i < NCh; i++) occ[i] = ModIntensity[i] / TotalIntensity[i];
                return occ;
            }
        }

        [Test]
        public static void PythonTranslation()
        {
            if (!File.Exists(FilePath)) { Assert.Ignore($"Input file not found (local-data test): {FilePath}"); return; }

            // ---- CELL 8/9: read the PSM table; locate the 18 reporter columns ----------------
            using var sr = new StreamReader(FilePath);
            var header = sr.ReadLine()?.Split('\t') ?? throw new Exception("Empty file");
            int Col(string n) { int i = Array.IndexOf(header, n); if (i < 0) throw new Exception($"Missing '{n}'"); return i; }
            int iFull = Col("Full Sequence"), iBase = Col("Base Sequence"), iAcc = Col("Accession"),
                iGene = Col("Gene Name"), iAmb = Col("Ambiguity Level"), iDct = Col("Decoy/Contaminant/Target"),
                iQ = Col("QValue"), iSE = Col("Start and End Residues In Full Sequence");
            int tmt1 = Col("126"), tmt18 = Col("135N");            // reporter block is tmt1..tmt18 (18 wide)
            var channelNames = Enumerable.Range(0, NCh).Select(i => header[tmt1 + i]).ToArray();
            var inv = CultureInfo.InvariantCulture;

            // Single pass: (a) build the raw (Accession, Base Sequence) -> Start/End map used by
            // CELL 32 (last occurrence in file wins, exactly like the Python dict), and
            // (b) apply the CELL 11 + CELL 16 row filters, keeping the reporter intensities.
            var seqToPos = new Dictionary<(string, string), string>();
            var rows = new List<Row>();
            string line;
            while ((line = sr.ReadLine()) != null)
            {
                var p = line.Split('\t');
                if (p.Length < header.Length) continue;
                seqToPos[(p[iAcc], p[iBase])] = p[iSE];            // CELL 32 map, from ALL raw rows

                // CELL 11: 1% FDR, target only, unambiguous localization (Ambiguity Level == "1").
                if (!double.TryParse(p[iQ], NumberStyles.Float, inv, out double q) || !(q < 0.01)) continue;
                if (p[iDct] != "T" || p[iAmb] != "1") continue;

                // CELL 16: drop the PSM if ANY of the 18 channels is exactly 0.
                var ch = new double[NCh]; bool hasZero = false;
                for (int i = 0; i < NCh; i++)
                {
                    double.TryParse(p[tmt1 + i], NumberStyles.Float, inv, out double v);
                    ch[i] = v; if (v == 0) { hasZero = true; break; }
                }
                if (hasZero) continue;
                rows.Add(new Row { FullSeq = p[iFull], BaseSeq = p[iBase], Accession = p[iAcc], Gene = p[iGene], StartEnd = p[iSE], Ch = ch });
            }
            Console.WriteLine($"Shape after removing rows with any zeros in TMT channels: ({rows.Count}, {header.Length})");

            // ---- CELL 20: sample-loading normalization -------------------------------------
            // Scale each channel so its column total meets a common target (median/mean of totals).
            var colTotals = new double[NCh];
            foreach (var r in rows) for (int i = 0; i < NCh; i++) colTotals[i] += r.Ch[i];
            double target = Normalization == "mean" ? colTotals.Average() : Median(colTotals);
            var scale = colTotals.Select(t => target / t).ToArray();
            foreach (var r in rows) for (int i = 0; i < NCh; i++) r.Ch[i] *= scale[i];

            // ---- CELL 27: strip insignificant modifications from the full sequence ---------
            foreach (var r in rows) r.FullSeq = RemoveMods(r.FullSeq);

            // ---- CELL 32: standardize Start/End via the (Accession, Base Sequence) map -----
            rows = rows.Where(r => seqToPos.ContainsKey((r.Accession, r.BaseSeq))).ToList();
            foreach (var r in rows) r.StartEnd = seqToPos[(r.Accession, r.BaseSeq)];

            // ---- CELL 38: drop position- or sequence-ambiguous rows ("|" present) ----------
            rows = rows.Where(r => !r.StartEnd.Contains('|') && !r.FullSeq.Contains('|')).ToList();

            // ---- CELL 42: total reporter intensity at every (protein position, accession) --
            // This is the occupancy denominator. If a peptide starts at the protein N-terminus
            // (span.start == 1) we also fold in position 0 (the protein N-terminal residue).
            var positionTotals = new Dictionary<(int, string), double[]>();
            foreach (var r in rows)
            {
                var (start, end) = ParseSpan(r.StartEnd);
                int lo = start == 1 ? 0 : start;
                for (int pos = lo; pos <= end; pos++)
                {
                    var key = (pos, r.Accession);
                    if (!positionTotals.TryGetValue(key, out var acc)) positionTotals[key] = acc = new double[NCh];
                    for (int i = 0; i < NCh; i++) acc[i] += r.Ch[i];
                }
            }
            Console.WriteLine($"Number of positions with intensity data: {positionTotals.Count}");

            // ---- CELL 43: accumulate modified-intensity per site modification --------------
            var siteMods = new Dictionary<(string, int, string), SiteModification>();
            foreach (var r in rows)
            {
                var posmods = ParseModifications(r.FullSeq);
                var (start, _) = ParseSpan(r.StartEnd);
                foreach (var (pepPos, mods) in posmods)
                    foreach (var mod0 in mods)
                    {
                        // A peptide-N-terminal mod (pepPos==0) only counts when the peptide is at
                        // the protein N-terminus; otherwise it isn't a protein N-terminal site.
                        if (pepPos == 0 && start != 1) continue;
                        string mod = UniversalPhospho(mod0) ?? mod0;   // fold phospho-S variants together
                        int posInProtein = start + pepPos - 1;
                        var key = (mod, posInProtein, r.Accession);
                        if (!siteMods.TryGetValue(key, out var site))
                        {
                            var total = positionTotals[(posInProtein, r.Accession)];  // present by construction
                            siteMods[key] = site = new SiteModification
                            {
                                ModificationName = mod, ParentPosition = posInProtein,
                                ParentAccession = r.Accession, ParentGene = r.Gene,
                                TotalIntensity = (double[])total.Clone()
                            };
                        }
                        for (int i = 0; i < NCh; i++) site.ModIntensity[i] += r.Ch[i];
                        site.NumberOfPsms++;
                    }
            }
            Console.WriteLine($"Number of site modifications: {siteMods.Count}");

            // ---- CELL 49: keep changing sites (0 < summed occupancy < 18) and build output --
            var outputSites = siteMods.Values.Where(s => { double sum = s.CalculateOccupancy().Sum(); return sum > 0 && sum < 18; }).ToList();
            var outputDf = BuildOutput(outputSites, channelNames);
            Console.WriteLine($"Output DataFrame shape: ({outputDf.Count}, {4 + 4 * NCh})");

            if (WriteOutputForInspection)                              // CELL 51 (verification only)
            {
                string outPath = Path.Combine(Path.GetDirectoryName(FilePath)!, $"site_modification_occupancies_{RunSuffix}_csharp.tsv");
                WriteOutput(outPath, outputDf, channelNames);
                Console.WriteLine($"Wrote {outPath}");
            }
        }

        // One assembled output record (mirrors a row of `output_df`).
        private sealed class OutRow
        {
            public string Modification, Accession, Gene;
            public int Position;
            public double[] Mod, Total, Weighted, Occupancy;      // 18 each
        }

        // Builds the equivalent of the CELL 49 `output_df` (4 meta columns + 4x18 value columns).
        private static List<OutRow> BuildOutput(List<SiteModification> sites, string[] channelNames)
        {
            var res = new List<OutRow>(sites.Count);
            foreach (var s in sites)
            {
                var occ = s.CalculateOccupancy();
                double meanTotal = s.TotalIntensity.Average();    // total_intensities.mean(axis=0) per site
                var weighted = occ.Select(o => o * meanTotal).ToArray();
                res.Add(new OutRow
                {
                    Modification = s.ModificationName, Position = s.ParentPosition,
                    Accession = s.ParentAccession, Gene = s.ParentGene,
                    Mod = s.ModIntensity, Total = s.TotalIntensity, Weighted = weighted, Occupancy = occ
                });
            }
            return res;
        }

        // Writes the output table (sorted by Accession, Position, Modification as in CELL 51).
        private static void WriteOutput(string path, List<OutRow> df, string[] ch)
        {
            var inv = CultureInfo.InvariantCulture;
            using var w = new StreamWriter(path);
            var head = new List<string> { "Modification", "Position", "Accession", "Gene Name" };
            head.AddRange(ch.Select(c => c + " Mod Intensity"));
            head.AddRange(ch.Select(c => c + " Total Intensity"));
            head.AddRange(ch.Select(c => c + " Mean Sample Intensity-Weighted Site Abundance"));
            head.AddRange(ch.Select(c => c + " Occupancy"));
            w.WriteLine(string.Join('\t', head));

            foreach (var r in df.OrderBy(r => r.Accession, StringComparer.Ordinal).ThenBy(r => r.Position).ThenBy(r => r.Modification, StringComparer.Ordinal))
            {
                var sb = new StringBuilder();
                sb.Append(r.Modification).Append('\t').Append(r.Position).Append('\t').Append(r.Accession).Append('\t').Append(r.Gene);
                foreach (var block in new[] { r.Mod, r.Total, r.Weighted, r.Occupancy })
                    foreach (var v in block) sb.Append('\t').Append(v.ToString("R", inv));
                w.WriteLine(sb.ToString());
            }
        }

        // ---- helpers (the notebook's utility functions) ------------------------------------

        // remove_mods: strip every "[<insignificant class> ...]" bracket from the sequence.
        private static string RemoveMods(string seq)
        {
            foreach (var mod in InsignificantMods)
                seq = Regex.Replace(seq, @"\[" + Regex.Escape(mod) + @".*?\]", "");
            return seq;
        }

        // parse_modifications: walk the full sequence, mapping (residues seen so far) -> mod names.
        // base position 0 == a peptide N-terminal mod; k == a mod on the k-th residue.
        private static Dictionary<int, List<string>> ParseModifications(string full)
        {
            var posmods = new Dictionary<int, List<string>>();
            int basePos = 0, i = 0;
            while (i < full.Length)
            {
                if (full[i] == '[')
                {
                    var m = ModRe.Match(full, i);
                    if (m.Success && m.Index == i)                 // a real mod bracket anchored at i
                    {
                        string cls = m.Groups[1].Value.Trim(), name = m.Groups[2].Value.Trim();
                        if (!InsignificantMods.Contains(cls) && !InsignificantMods.Contains(name))
                        {
                            if (!posmods.TryGetValue(basePos, out var lst)) posmods[basePos] = lst = new List<string>();
                            lst.Add(name);
                        }
                        i = m.Index + m.Length;
                        continue;
                    }
                }
                basePos++; i++;                                    // ordinary residue
            }
            return posmods;
        }

        // stringToSpan: "[start to end]" -> (start, end).
        private static (int start, int end) ParseSpan(string s)
        {
            var m = SpanRe.Match(s);
            return (int.Parse(m.Groups[1].Value), int.Parse(m.Groups[2].Value));
        }

        // universal_phosphoserine: treat "Phosphorylation on S"/"Phosphoserine on S" as one mod.
        private static string UniversalPhospho(string m) =>
            m == "Phosphorylation on S" || m == "Phosphoserine on S" ? "Phosphoserine on S" : null;

        private static double Median(IReadOnlyList<double> vals)
        {
            var a = vals.OrderBy(x => x).ToArray();
            int n = a.Length;
            return n == 0 ? 0 : n % 2 == 1 ? a[n / 2] : 0.5 * (a[n / 2 - 1] + a[n / 2]);
        }
    }

    //var allSitePsmGroups = new List<PtmSiteGroup>();
    //        foreach (var ptm in allPtms.Distinct())
    //        {
    //            var sitePsms = new List<PsmFromTsv>();
    //var allPsmsFromTheProtein = allPsms.Where(p => p.ProteinAccession == ptm.Accession);
    //            foreach(var psm in allPsmsFromTheProtein)
    //            {
    //                int startAA = int.Parse(psm.StartAndEndResiduesInProtein.Split()[0].Split("[")[1]);
    //int endAA = int.Parse(psm.StartAndEndResiduesInProtein.Split()[2].Split("]")[0]);
    //                if ((ptm.ModSite >= startAA && ptm.ModSite <= endAA) || (ptm.ModSite == 0 && startAA == 1))
    //                {
    //                    var allModsFromThisPsm = SpectrumMatchFromTsv.ParseModifications(psm.FullSequence).Where(kvp => !notInteresting.Any(key => kvp.Value.Contains(key)));
    //                    if (allModsFromThisPsm.Any(kvp => kvp.Key + startAA == ptm.ModSite && kvp.Value == ptm.Mod) || (ptm.ModSite == 0 && !allModsFromThisPsm.Any()))
    //                    sitePsms.Add(psm);
    //                }
    //            }
    //            //allSitePsmGroups.Add(new PtmSiteGroup { PTM = ptm, AllPsms = sitePsms });
    //        }
}
