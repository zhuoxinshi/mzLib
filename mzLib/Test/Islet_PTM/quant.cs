//using NUnit.Framework;
//using System;
//using System.Collections.Generic;
//using System.Globalization;
//using System.IO;
//using System.Linq;
//using System.Text;
//using System.Threading.Tasks;

//namespace Test.Islet_PTM
//{
//    // =============================================================================
//    // TMT peptide-level quant from PSMs (FragPipe TMT-Integrator style), C# port.
//    //
//    // Steps:
//    //   1. Per PSM: raw 16-channel intensities (134C, 135N are reference channels
//    //      and ignored everywhere). Reference = median of those 16 channels;
//    //      compute log2(channel / reference) for each of the 16.
//    //   2. Aggregate PSM log2-ratios to peptide (FullSequence) by the median,
//    //      within each TMT batch (R2 / R4 / R6 file).
//    //   3. Build a peptide x 48-sample table, with samples ordered exactly as
//    //      they appear in design.xlsx.
//    //   4. Align: median-center each of the 48 sample columns in log2 space.
//    //   5. Convert back to intensity:
//    //         intensity[sample] = 2^aligned_log2[sample] * scalar
//    //      where scalar = sum of reporter intensities across every PSM,
//    //      every 16 sample channels, every batch where that peptide was seen.
//    //      One scalar per peptide => the 48 columns stay on a comparable scale.
//    //
//    // Usage:
//    //   dotnet run --project . -- (or compile the file standalone with `dotnet script`/csc)
//    //   Paths are hard-coded below — edit if the layout changes.
//    // =============================================================================


//    internal static class TMTPeptideQuant
//    {
//        // 18 reporter channels in TSV column order
//        private static readonly string[] AllChannels =
//        {
//        "126","127N","127C","128N","128C","129N","129C","130N","130C",
//        "131N","131C","132N","132C","133N","133C","134N","134C","135N"
//    };

//        // Sample channels = first 16 (drop refs 134C, 135N)
//        private static readonly string[] SampleChannels = AllChannels.Take(16).ToArray();

//        private const string InDir =
//            @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\UPLC\LFgptmdPrunedDb-AAsub_search-cali-search";

//        private const string DesignXlsx =
//            @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\DataAnalysis\design.xlsx";

//        private const string OutTsv =
//            @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\DataAnalysis\Peptide_Intensities_48samples.tsv";

//        private static readonly (string Batch, string Path)[] Files =
//        {
//        ("R2", Path.Combine(InDir, "AllPSMs_withRI_AAsub_R2.tsv")),
//        ("R4", Path.Combine(InDir, "AllPSMs_withRI_AAsub_R4.tsv")),
//        ("R6", Path.Combine(InDir, "AllPSMs_withRI_AAsub_R6.tsv"))
//    };

//        // Batch -> starting index in the 48-sample layout (R2: 0-15, R4: 16-31, R6: 32-47)
//        private static readonly Dictionary<string, int> BatchStart = new()
//        {
//            ["R2"] = 0,
//            ["R4"] = 16,
//            ["R6"] = 32
//        };

//        private sealed class PeptideData
//        {
//            public string BaseSequence = "";
//            public string Accession = "";
//            public int PsmCount;
//            public double TotalReporterSum;                       // scalar accumulator
//            public readonly List<double>[] RatiosPerSample = new List<double>[48];

//            public PeptideData()
//            {
//                for (int i = 0; i < 48; i++) RatiosPerSample[i] = new List<double>();
//            }
//        }

//        public static void Main()
//        {
//            Console.WriteLine("Reading design.xlsx ...");
//            var sampleNames = ReadDesignSamples(DesignXlsx);     // 48 names, design order
//            if (sampleNames.Count != 48)
//                throw new Exception($"Expected 48 samples in design, got {sampleNames.Count}");

//            var pepData = new Dictionary<string, PeptideData>(capacity: 200_000);

//            foreach (var (batch, path) in Files)
//            {
//                Console.WriteLine($"Processing {batch}: {Path.GetFileName(path)}");
//                ProcessPsmFile(batch, path, pepData);
//            }
//            Console.WriteLine($"Total peptides (by FullSequence): {pepData.Count:N0}");

//            // Build peptide x 48 log2-ratio matrix (median of PSM log2-ratios per cell)
//            var peptideKeys = pepData.Keys.ToArray();
//            int P = peptideKeys.Length;
//            var matrix = new double[P, 48];
//            for (int p = 0; p < P; p++)
//            {
//                var pd = pepData[peptideKeys[p]];
//                for (int s = 0; s < 48; s++)
//                {
//                    matrix[p, s] = pd.RatiosPerSample[s].Count > 0
//                        ? Median(pd.RatiosPerSample[s])
//                        : double.NaN;
//                }
//            }

//            // Median-center each of the 48 sample columns (in log2 space)
//            Console.WriteLine("Aligning 48 sample distributions (median centering in log2)...");
//            var preMedians = new double[48];
//            var postMedians = new double[48];
//            for (int s = 0; s < 48; s++)
//            {
//                var col = new List<double>(P);
//                for (int p = 0; p < P; p++)
//                    if (!double.IsNaN(matrix[p, s])) col.Add(matrix[p, s]);
//                double med = col.Count > 0 ? Median(col) : 0.0;
//                preMedians[s] = med;
//                for (int p = 0; p < P; p++)
//                    if (!double.IsNaN(matrix[p, s])) matrix[p, s] -= med;
//                // Post for sanity
//                var col2 = new List<double>(col.Count);
//                for (int p = 0; p < P; p++)
//                    if (!double.IsNaN(matrix[p, s])) col2.Add(matrix[p, s]);
//                postMedians[s] = col2.Count > 0 ? Median(col2) : 0.0;
//            }
//            Console.WriteLine($"  pre-centering medians  range: [{preMedians.Min():F4}, {preMedians.Max():F4}]");
//            Console.WriteLine($"  post-centering medians range: [{postMedians.Min():E2}, {postMedians.Max():E2}]");

//            // Write output TSV with columns in design order
//            Console.WriteLine($"Writing {OutTsv}");
//            using var w = new StreamWriter(OutTsv);
//            var header = new List<string> { "FullSequence", "BaseSequence", "Accession", "nPSM_total", "scalar" };
//            header.AddRange(sampleNames);
//            w.WriteLine(string.Join('\t', header));

//            var inv = CultureInfo.InvariantCulture;
//            for (int p = 0; p < P; p++)
//            {
//                var pd = pepData[peptideKeys[p]];
//                var cells = new List<string>(5 + 48)
//            {
//                peptideKeys[p],
//                pd.BaseSequence,
//                pd.Accession,
//                pd.PsmCount.ToString(inv),
//                pd.TotalReporterSum.ToString("R", inv)
//            };
//                for (int s = 0; s < 48; s++)
//                {
//                    double v = matrix[p, s];
//                    if (double.IsNaN(v))
//                        cells.Add("");                            // NA where peptide unobserved
//                    else
//                        cells.Add((Math.Pow(2.0, v) * pd.TotalReporterSum).ToString("R", inv));
//                }
//                w.WriteLine(string.Join('\t', cells));
//            }
//            Console.WriteLine($"Done. Wrote {P:N0} peptides x 48 samples.");
//        }

//        [Test]
//        public static void PeptideQuant()
//        {

//        }
//        // ------------------------------------------------------------------------
//        // Per-PSM processing: parse TSV, compute log2-ratios, accumulate into pepData
//        // ------------------------------------------------------------------------
//        private static void ProcessPsmFile(string batch, string path,
//                                           Dictionary<string, PeptideData> pepData)
//        {
//            using var sr = new StreamReader(path);
//            var headerLine = sr.ReadLine()
//                ?? throw new Exception($"Empty file: {path}");
//            var headerCols = headerLine.Split('\t');

//            int idxFull = Array.IndexOf(headerCols, "FullSequence");
//            int idxBase = Array.IndexOf(headerCols, "BaseSequence");
//            int idxAcc = Array.IndexOf(headerCols, "Accession");
//            if (idxFull < 0 || idxBase < 0 || idxAcc < 0)
//                throw new Exception($"Missing expected columns in {path}");

//            var channelIdx = SampleChannels.Select(c => Array.IndexOf(headerCols, c)).ToArray();
//            if (channelIdx.Any(i => i < 0))
//                throw new Exception($"Missing one or more sample-channel columns in {path}");

//            int batchStart = BatchStart[batch];
//            var inv = CultureInfo.InvariantCulture;
//            var vals = new double[16];

//            string? line;
//            int psms = 0, skipped = 0;
//            while ((line = sr.ReadLine()) != null)
//            {
//                var parts = line.Split('\t');
//                if (parts.Length < headerCols.Length) { skipped++; continue; }

//                // Parse 16 sample-channel intensities; <=0 / NA treated as missing.
//                int nValid = 0;
//                double rowSum = 0.0;
//                for (int i = 0; i < 16; i++)
//                {
//                    if (double.TryParse(parts[channelIdx[i]], NumberStyles.Float, inv, out double v)
//                        && v > 0 && !double.IsInfinity(v))
//                    {
//                        vals[i] = v;
//                        rowSum += v;
//                        nValid++;
//                    }
//                    else vals[i] = double.NaN;
//                }
//                if (nValid == 0) { skipped++; continue; }

//                double refMed = MedianNonNaN(vals);
//                if (!(refMed > 0)) { skipped++; continue; }

//                string fullSeq = parts[idxFull];
//                if (!pepData.TryGetValue(fullSeq, out var pd))
//                {
//                    pd = new PeptideData
//                    {
//                        BaseSequence = parts[idxBase],
//                        Accession = parts[idxAcc]
//                    };
//                    pepData[fullSeq] = pd;
//                }
//                pd.PsmCount++;
//                pd.TotalReporterSum += rowSum;

//                for (int i = 0; i < 16; i++)
//                {
//                    if (double.IsNaN(vals[i])) continue;
//                    double lr = Math.Log2(vals[i] / refMed);
//                    pd.RatiosPerSample[batchStart + i].Add(lr);
//                }
//                psms++;
//            }
//        }

//        // ------------------------------------------------------------------------
//        // Median helpers
//        // ------------------------------------------------------------------------
//        private static double Median(List<double> values)
//        {
//            var arr = values.ToArray();
//            Array.Sort(arr);
//            int n = arr.Length;
//            if (n == 0) return double.NaN;
//            return (n % 2 == 1) ? arr[n / 2] : 0.5 * (arr[n / 2 - 1] + arr[n / 2]);
//        }

//        private static double MedianNonNaN(double[] values)
//        {
//            var list = new List<double>(values.Length);
//            foreach (var v in values) if (!double.IsNaN(v)) list.Add(v);
//            return Median(list);
//        }

//    }
