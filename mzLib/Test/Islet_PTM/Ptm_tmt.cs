using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Readers;
using NUnit.Framework;
using System.Windows.Markup;
using Omics;
using System.Data.Entity.Core.Common.CommandTrees.ExpressionBuilder;
using System.IO;
using System.Text.RegularExpressions;
using Easy.Common.Extensions;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.SupportedModels.RetentionTimeModels;
using TopDownProteomics;
using System.Globalization;

namespace Test
{
    public class Ptm_tmt
    {
        [Test]
        public static void HighResTMT()
        {
            var notInteresting = new List<string> { "TMT18", "Fixed", "Artifact", "Variable", "Metal" };

            var allPeptidesHigh_path = @"E:\Fly_TMT\MM\calied_gptmd_secondPass\Task2-SearchTask\AllPeptides.psmtsv";
            var allPeptidesHigh_file = new PsmFromTsvFile(allPeptidesHigh_path);
            var allPeptidesHigh = allPeptidesHigh_file.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T");
            var allPeptidesNoModHigh = allPeptidesHigh.Where(p => !SpectrumMatchFromTsv.ParseModifications(p.FullSequence).Values.Any(v => v.Contains("Fixed") || v.Contains("Variable") || v.Contains("TMT")) && !p.FullSequence.Contains("|"));
            var allPeptidesWithModsHigh = allPeptidesHigh.Where(p => SpectrumMatchFromTsv.ParseModifications(p.FullSequence).Values.Any(v => !notInteresting.Any(key => v.Contains(key)))).ToList();

            var allPeptidesLow_path = @"E:\Islets\Brian_data\Real_islets\Frxn\All_LFgptmdFilterPrunedDb\Task1-SearchTask\AllPeptides.psmtsv";
            var allPeptidesLow_file = new PsmFromTsvFile(allPeptidesLow_path);
            var allPeptidesLow = allPeptidesLow_file.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T");
            //var allPeptidesWithModsLow = allPeptidesLow.Where(p => SpectrumMatchFromTsv.ParseModifications(p.FullSequence).Values.Any(v => !notInteresting.Any(key => v.Contains(key))) || p.Description.Contains("chain")).ToList();

            var rtPairs = new List<(double high, double low)>();
            var highResOverlap = new List<PsmFromTsv>();
            foreach (var peptide in allPeptidesWithModsHigh)
            {
                var lowResPeptide = allPeptidesLow.FirstOrDefault(p => p.FullSequence == peptide.FullSequence);
                if (lowResPeptide != null)
                {
                    rtPairs.Add((peptide.RetentionTime, lowResPeptide.RetentionTime));
                    highResOverlap.Add(lowResPeptide);
                }
            }

            var lowResOutPath = @"E:\Islets\Brian_data\Real_islets\Frxn\All_LFgptmdFilterPrunedDb\ModifiedPeptidesTMT_HighResOverlap_BestPsm.tsv";
            TmtPair.WriteResults(highResOverlap, lowResOutPath);
        }

        [Test]
        public static void TMT_BestPsm()
        {
            var notInteresting = new List<string> { "TMT18", "Fixed", "Artifact", "Variable", "Metal" };

            string allPeptidesTMT_path = @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\UPLC\MM\LFgptmdPrunedDb_search-cali-search\Task1-SearchTask\AllPeptides.psmtsv";
            var allPeptidesTMT_file = new PsmFromTsvFile(allPeptidesTMT_path);
            var allPeptidesTMT = allPeptidesTMT_file.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T");
            var allPeptidesTMTWithMods = allPeptidesTMT.Where(p => SpectrumMatchFromTsv.ParseModifications(p.FullSequence).Values.Any(v => !notInteresting.Any(key => v.Contains(key))) || p.Description.Contains("chain")).ToList();

            var model = new Prosit2020iRTTMT();
            var rtPredictionOutPath = @"E:\Islets\Brian_data\Real_islets\Frxn\All_gptmdPrunedDb-second\TMT_BestPsm_RtPrediction.tsv";
            var rtColumns = new List<string> { "FullSequence", "Mod", "ObservedRT", "PredictedRT" };
            using (StreamWriter writer = new StreamWriter(rtPredictionOutPath))
            {
                writer.WriteLine(string.Join("\t", rtColumns));
                foreach (var peptide in allPeptidesTMT)
                {
                    var predictionInput = new List<RetentionTimePredictionInput> { new RetentionTimePredictionInput(peptide.FullSequence.Replace("X", "N-terminus")) };
                    if (predictionInput.First().SequenceWarning != null)
                    {
                        continue;
                    }
                    var predictedRt = model.Predict(predictionInput).First().PredictedRetentionTime;
                    var mod = "Unmodified";
                    if (!SpectrumMatchFromTsv.ParseModifications(peptide.FullSequence).Values.All(v => v.Contains("Fixed") || v.Contains("Variable")))
                    {
                        mod = "GptmdMod";
                        if (peptide.FullSequence.Contains("Biological") || peptide.FullSequence.Contains("Uniprot"))
                        {
                            mod = "BioMod";
                            if (peptide.FullSequence.Contains("Phospho"))
                            {
                                mod = "Phospho";
                            }
                        }
                    }
                    var outString = new List<string> { peptide.FullSequence, mod, peptide.RetentionTime.ToString(), predictedRt.ToString() };
                    writer.WriteLine(string.Join("\t", outString));
                }
            }
            //string outPath = @"E:\Islets\Brian_data\Real_islets\Frxn\All_gptmdPrunedDb-second\ModifiedPeptidesTMT.tsv";
            //TmtPair.WriteResults(allPeptidesTMTWithMods, outPath);
        }

        [Test]
        public static void TMT_SumPsms()
        {
            var notInteresting = new List<string> { "TMT18", "Fixed", "Artifact", "Variable", "Metal" };

            string allPeptidesTMT_path = @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\UPLC\MM\LFgptmdPrunedDb_search-cali-search\Task1-SearchTask\AllPSMs.psmtsv";
            var allPeptidesTMT_file = new PsmFromTsvFile(allPeptidesTMT_path);
            var allPsms = allPeptidesTMT_file.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T");
            var psmGroup = allPsms.Where(p => p.FileName.Contains("R6"));

            string outPath = @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\UPLC\MM\LFgptmdPrunedDb_search-cali-search\calied_Peptides_PsmSum_R6-R7.tsv";
            var columns = new List<string> { "BaseSequence", "FullSequence", "Mods", "Protein Accession", "Protein Name", "GeneName", "Description", "Start_and_End_Residues" };
            var labels = new List<string> { "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C", "132N", "132C", "133N", "133C", "134N", "134C", "135N" };
            columns.AddRange(labels);
            using (StreamWriter writer = new StreamWriter(outPath))
            {
                writer.WriteLine(string.Join("\t", columns));
                var peptides = psmGroup.GroupBy(p => p.FullSequence);
                foreach (var peptide in peptides)
                {
                    var mods = String.Join(", ", SpectrumMatchFromTsv.ParseModifications(peptide.First().FullSequence).Values.Where(m => !notInteresting.Any(x => m.Contains(x))));
                    var outString = new List<string> { peptide.First().BaseSequence, peptide.First().FullSequence, mods, peptide.First().ProteinAccession, peptide.First().ProteinName, peptide.First().GeneName, peptide.First().Description, peptide.First().StartAndEndResiduesInProtein };

                    var reporterIonIntensities = new double[peptide.First().Intensities.Count()];
                    foreach (var psm in peptide)
                    {
                        for (int i = 0; i < peptide.First().Intensities.Count(); i++)
                        {
                            reporterIonIntensities[i] += psm.Intensities[i];
                        }
                    }
                    outString.AddRange(reporterIonIntensities.Select(i => i.ToString()));
                    writer.WriteLine(string.Join("\t", outString));
                }
            }
        }

        [Test]
        public static void TMT_SumPsms_multiPlates()
        {
            var dir = @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\UPLC\LFgptmdPrunedDb-noFilter_search\Task1-SearchTask";
            var individualDir = @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\UPLC\MM\LFgptmdPrunedDb_search-cali-search\Task1-SearchTask\Individual File Results";
            var allIndividualFiles = Directory.GetFiles(individualDir, "*PSMs.psmtsv", SearchOption.AllDirectories);
            var allPsms = @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\UPLC\LFgptmdPrunedDb-noFilter_search\Task1-SearchTask\AllPSMs.psmtsv";
            var allPsms_file = new PsmFromTsvFile(allPsms);
            var filteredPsms = allPsms_file.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T");
            var keys = new List<string> { "R2-R3", "R4-R5", "R6-R7"};

            var columns = new List<string> { "BaseSequence", "FullSequence", "Protein Accession", "Protein Name", "GeneName", "Description", "Start_and_End_Residues" };
            var labels = new List<string> { "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C", "132N", "132C", "133N", "133C", "134N", "134C", "135N" };

            foreach (var key in keys)
            {
                var psmsForKey = filteredPsms.Where(p => p.FileName.Contains(key));
                string outPath = Path.Combine(dir, $"Peptide_quant_sumPSM_{key}.tsv");
                using (StreamWriter writer = new StreamWriter(outPath))
                {
                    writer.WriteLine(string.Join("\t", columns));
                    var peptides = psmsForKey.GroupBy(p => p.FullSequence);
                    foreach (var peptide in peptides)
                    {
                        var outString = new List<string> { peptide.First().BaseSequence, peptide.First().FullSequence, peptide.First().ProteinAccession, peptide.First().ProteinName, peptide.First().GeneName, peptide.First().Description, peptide.First().StartAndEndResiduesInProtein };

                        var reporterIonIntensities = new double[peptide.First().Intensities.Count()];
                        foreach (var psm in peptide)
                        {
                            for (int i = 0; i < peptide.First().Intensities.Count(); i++)
                            {
                                reporterIonIntensities[i] += psm.Intensities[i];
                            }
                        }
                        //var bestPsm = peptide.OrderByDescending(p => p.Score).First();
                        //var reporterIonIntensities = bestPsm.Intensities;
                        
                        outString.AddRange(reporterIonIntensities.Select(i => i.ToString()));
                        writer.WriteLine(string.Join("\t", outString));
                    }
                }
            }
        }

        [Test]
        public static void TMT_ConvertIntensitiesToRatios()
        {
            var dir = @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\UPLC\MM\LFgptmdPrunedDb_search-cali-search\Task1-SearchTask";
            var allPsms = @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\UPLC\MM\LFgptmdPrunedDb_search-cali-search\Task1-SearchTask\AllPSMs.psmtsv";
            var outPath = Path.Combine(dir, "AllPSMs_ratios.psmtsv");

            using (var reader = new StreamReader(allPsms))
            using (var writer = new StreamWriter(outPath))
            {
                // Header is written back verbatim and used to locate the columns we need.
                var headerLine = reader.ReadLine();
                writer.WriteLine(headerLine);
                var header = headerLine.Split('\t');
                var headerIndex = new Dictionary<string, int>();
                for (int i = 0; i < header.Length; i++)
                {
                    headerIndex[header[i].Trim()] = i;
                }

                var channelIndices = SpectrumMatchFromTsvHeader.TmtChannelNames.Where(name => headerIndex.ContainsKey(name)).Select(name => headerIndex[name]).ToList();

                int qValueCol = headerIndex[SpectrumMatchFromTsvHeader.QValue];
                int targetCol = headerIndex[SpectrumMatchFromTsvHeader.DecoyContaminantTarget];

                string line;
                while ((line = reader.ReadLine()) != null)
                {
                    if (string.IsNullOrWhiteSpace(line)) continue;
                    var spl = line.Split('\t');

                    // Same filter as the original filteredPsms: 1% FDR targets only.
                    if (!double.TryParse(spl[qValueCol].Trim(), NumberStyles.Any, CultureInfo.InvariantCulture, out double qValue) || qValue > 0.01 || spl[targetCol].Trim() != "T")
                    {
                        continue;
                    }

                    // Ratio of each channel = channel intensity / mean of the row (across all channels).
                    var intensities = channelIndices.Select(c => double.TryParse(spl[c].Trim(), NumberStyles.Any, CultureInfo.InvariantCulture, out double v) ? v : 0.0).ToList();
                    double mean = intensities.Average();

                    if (mean != 0)
                    {
                        for (int i = 0; i < channelIndices.Count; i++)
                        {
                            spl[channelIndices[i]] = (intensities[i] / mean).ToString(CultureInfo.InvariantCulture);
                        }
                    }
                    writer.WriteLine(string.Join("\t", spl));
                }
            }
        }


        private static double Median(IEnumerable<double> values)
        {
            var sorted = values.OrderBy(x => x).ToArray();
            int n = sorted.Length;
            if (n == 0) return 0;
            return n % 2 == 1 ? sorted[n / 2] : 0.5 * (sorted[n / 2 - 1] + sorted[n / 2]);
        }

        // Median of each channel across a set of ratio vectors; null per channel when there are no vectors.
        private static double?[] MedianPerChannel(IEnumerable<double[]> vectors, int nCh)
        {
            var list = vectors.ToList();
            var res = new double?[nCh];
            if (list.Count == 0) return res;                    // all null -> blank cells in the output
            for (int i = 0; i < nCh; i++)
                res[i] = Median(list.Select(v => v[i]));
            return res;
        }

        // Sum of each channel across a set of intensity vectors; zeros when there are no vectors.
        private static double[] SumPerChannel(IEnumerable<double[]> vectors, int nCh)
        {
            var res = new double[nCh];
            foreach (var v in vectors)
                for (int i = 0; i < nCh; i++)
                    res[i] += v[i];
            return res;
        }

        [Test]
        public static void TMT_RatioAggregation()
        {
            var dir = @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\UPLC\MM\LFgptmdPrunedDb_search-cali-search\Task1-SearchTask";
            var allPsms_path = @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\UPLC\MM\LFgptmdPrunedDb_search-cali-search\Task1-SearchTask\AllPSMs_R2-R3.psmtsv";
            var allPsms_file = new PsmFromTsvFile(allPsms_path);

            var labels = new List<string> { "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C", "132N", "132C", "133N", "133C", "134N", "134C", "135N" };
            var startEndParser = new Regex(@"\[(\d+)\s+to\s+(\d+)\]");
            int nCh = labels.Count;

            // All 1% FDR target PSMs with reporter ions (the set used for channel normalization).
            var targetPsms = allPsms_file.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T" && p.Intensities != null && p.Intensities.Length == nCh).ToList();

            // Sample-loading normalization (before any further filtering): scale each channel so that
            // its summed intensity across all 1% FDR target PSMs is equal across the 18 channels.
            var channelSums = new double[nCh];
            foreach (var p in targetPsms)
                for (int i = 0; i < nCh; i++)
                    channelSums[i] += p.Intensities[i];

            double targetSum = channelSums.Average();   // common total to scale every channel to
            var scale = new double[nCh];
            for (int i = 0; i < nCh; i++)
                scale[i] = channelSums[i] > 0 ? targetSum / channelSums[i] : 1.0;

            foreach (var p in targetPsms)
                for (int i = 0; i < nCh; i++)
                    p.Intensities[i] *= scale[i];

            // Keep unambiguous full sequences (single site mapping).
            var filteredPsms = targetPsms
                .Where(p => !string.IsNullOrEmpty(p.FullSequence) && !p.FullSequence.Contains("|")
                            //&& !string.IsNullOrEmpty(p.ProteinAccession) && !p.ProteinAccession.Contains("|")
                            //&& !string.IsNullOrEmpty(p.StartAndEndResiduesInProtein) && !p.StartAndEndResiduesInProtein.Contains("|")
                            )
                .ToList();

            // Convert each PSM's reporter intensities to within-row ratios (channel / mean of the row),
            // then record its protein span and the (protein position, mod) bio-mods it carries.
            var psmInfos = new List<(PsmFromTsv psm, int start, int end, string acc, string gene, HashSet<(int pos, string mod)> carried)>();
            foreach (var psm in filteredPsms)
            {
                double mean = psm.Intensities.Average();
                if (mean == 0) continue;                        // all-zero row: no usable ratio
                for (int i = 0; i < psm.Intensities.Length; i++)
                    psm.Intensities[i] = psm.Intensities[i] / mean;

                var m = startEndParser.Match(psm.StartAndEndResiduesInProtein);
                if (!m.Success) continue;
                int start = int.Parse(m.Groups[1].Value);
                int end = int.Parse(m.Groups[2].Value);

                // Bio-mods this PSM carries, mapped to protein positions (same formula as Occupancy).
                var carried = new HashSet<(int pos, string mod)>();
                foreach (var mod in SpectrumMatchFromTsv.ParseModifications(psm.FullSequence).Where(kvp => kvp.Value.Contains("Common Biological")))
                {
                    if (mod.Key == 0 && start != 1) continue;   // N-terminal mod only counts at protein N-term
                    int pos = mod.Key + start;
                    if (mod.Key != 0) pos = pos - 1;
                    carried.Add((pos, mod.Value.Split(':')[1].Trim()));
                }
                psmInfos.Add((psm, start, end, psm.ProteinAccession, psm.GeneName, carried));
            }

            // Every distinct modification site = (accession, protein position, mod).
            var sites = psmInfos.SelectMany(x => x.carried.Select(c => (acc: x.acc, pos: c.pos, mod: c.mod))).Distinct().ToList();

            // Index PSMs by accession so the span lookup only scans the relevant protein.
            var psmsByAcc = psmInfos.GroupBy(x => x.acc).ToDictionary(g => g.Key, g => g.ToList());

            string outPath = Path.Combine(dir, "Normed_site_ratio_median_R2-R3_accessionAmb_siteAmb.tsv");
            var columns = new List<string> { "Modification", "Accession", "Position", "GeneName", "NumModPsms", "NumUnmodPsms" };
            columns.AddRange(labels.Select(l => "Mod_" + l));
            columns.AddRange(labels.Select(l => "Unmod_" + l));

            using (var writer = new StreamWriter(outPath))
            {
                writer.WriteLine(string.Join("\t", columns));
                foreach (var site in sites.OrderBy(s => s.acc).ThenBy(s => s.pos).ThenBy(s => s.mod))
                {
                    var accPsms = psmsByAcc[site.acc];

                    // PSMs carrying the mod at this site vs. PSMs spanning the site without it.
                    var modPsms = accPsms.Where(x => x.carried.Contains((site.pos, site.mod))).ToList();
                    var unmodPsms = accPsms.Where(x => x.start <= site.pos && site.pos <= x.end
                                                      && !x.carried.Contains((site.pos, site.mod))).ToList();

                    var gene = modPsms[0].gene;
                    var modMedian = MedianPerChannel(modPsms.Select(x => x.psm.Intensities), labels.Count);
                    var unmodMedian = MedianPerChannel(unmodPsms.Select(x => x.psm.Intensities), labels.Count);

                    var outString = new List<string> { site.mod, site.acc, site.pos.ToString(), gene, modPsms.Count.ToString(), unmodPsms.Count.ToString() };
                    outString.AddRange(modMedian.Select(v => v.HasValue ? v.Value.ToString(CultureInfo.InvariantCulture) : ""));
                    outString.AddRange(unmodMedian.Select(v => v.HasValue ? v.Value.ToString(CultureInfo.InvariantCulture) : ""));
                    writer.WriteLine(string.Join("\t", outString));
                }
            }
        }

        [Test]
        public static void TMT_Occupancy() { 
            var dir = @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\UPLC\MM\LFgptmdPrunedDb_search-cali-search\Task1-SearchTask";
            var allPsms_path = @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\UPLC\MM\LFgptmdPrunedDb_search-cali-search\Task1-SearchTask\AllPSMs_R6-R7.psmtsv";
            var allPsms_file = new PsmFromTsvFile(allPsms_path);

            var labels = new List<string> { "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C", "132N", "132C", "133N", "133C", "134N", "134C", "135N" };
            var startEndParser = new Regex(@"\[(\d+)\s+to\s+(\d+)\]");
            int nCh = labels.Count;

            // All 1% FDR target PSMs with reporter ions (the set used for channel normalization).
            var targetPsms = allPsms_file.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T" && p.Intensities != null && p.Intensities.Length == nCh).ToList();

            // Sample-loading normalization (before any further filtering): scale each channel so that
            // its summed intensity across all 1% FDR target PSMs is equal across the 18 channels.
            var channelSums = new double[nCh];
            foreach (var p in targetPsms)
                for (int i = 0; i < nCh; i++)
                    channelSums[i] += p.Intensities[i];

            double targetSum = channelSums.Average();   // common total to scale every channel to
            var scale = new double[nCh];
            for (int i = 0; i < nCh; i++)
                scale[i] = channelSums[i] > 0 ? targetSum / channelSums[i] : 1.0;

            foreach (var p in targetPsms)
                for (int i = 0; i < nCh; i++)
                    p.Intensities[i] *= scale[i];

            // Keep unambiguous full sequences (single site mapping).
            var filteredPsms = targetPsms
                .Where(p => !string.IsNullOrEmpty(p.FullSequence) && !p.FullSequence.Contains("|")
                            //&& !string.IsNullOrEmpty(p.ProteinAccession) && !p.ProteinAccession.Contains("|")
                            //&& !string.IsNullOrEmpty(p.StartAndEndResiduesInProtein) && !p.StartAndEndResiduesInProtein.Contains("|")
                            )
                .ToList();

            // Record each PSM's protein span and the (protein position, mod) bio-mods it carries.
            // Reporter intensities are kept as (normalized) intensities — NOT converted to ratios —
            // so they can be summed per site below.
            var psmInfos = new List<(PsmFromTsv psm, int start, int end, string acc, string gene, HashSet<(int pos, string mod)> carried)>();
            foreach (var psm in filteredPsms)
            {
                var m = startEndParser.Match(psm.StartAndEndResiduesInProtein);
                if (!m.Success) continue;
                int start = int.Parse(m.Groups[1].Value);
                int end = int.Parse(m.Groups[2].Value);

                // Bio-mods this PSM carries, mapped to protein positions (same formula as Occupancy).
                var carried = new HashSet<(int pos, string mod)>();
                foreach (var mod in SpectrumMatchFromTsv.ParseModifications(psm.FullSequence).Where(kvp => kvp.Value.Contains("Common Biological")))
                {
                    if (mod.Key == 0 && start != 1) continue;   // N-terminal mod only counts at protein N-term
                    int pos = mod.Key + start;
                    if (mod.Key != 0) pos = pos - 1;
                    carried.Add((pos, mod.Value.Split(':')[1].Trim()));
                }
                psmInfos.Add((psm, start, end, psm.ProteinAccession, psm.GeneName, carried));
            }

            // Every distinct modification site = (accession, protein position, mod).
            var sites = psmInfos.SelectMany(x => x.carried.Select(c => (acc: x.acc, pos: c.pos, mod: c.mod))).Distinct().ToList();

            // Index PSMs by accession so the span lookup only scans the relevant protein.
            var psmsByAcc = psmInfos.GroupBy(x => x.acc).ToDictionary(g => g.Key, g => g.ToList());

            string outPath = Path.Combine(dir, "site_occupancy_R6-R7.tsv");
            var columns = new List<string> { "Modification", "Accession", "Position", "GeneName", "NumModPsms", "NumUnmodPsms" };
            columns.AddRange(labels.Select(l => "Mod_" + l));
            columns.AddRange(labels.Select(l => "Unmod_" + l));
            columns.AddRange(labels.Select(l => "Occupancy_" + l));

            using (var writer = new StreamWriter(outPath))
            {
                writer.WriteLine(string.Join("\t", columns));
                foreach (var site in sites.OrderBy(s => s.acc).ThenBy(s => s.pos).ThenBy(s => s.mod))
                {
                    var accPsms = psmsByAcc[site.acc];

                    // PSMs carrying the mod at this site vs. PSMs spanning the site without it.
                    var modPsms = accPsms.Where(x => x.carried.Contains((site.pos, site.mod))).ToList();
                    var unmodPsms = accPsms.Where(x => x.start <= site.pos && site.pos <= x.end
                                                      && !x.carried.Contains((site.pos, site.mod))).ToList();

                    var gene = modPsms[0].gene;

                    // Per-site quant = summed reporter intensity of its PSMs, per channel.
                    var modSum = SumPerChannel(modPsms.Select(x => x.psm.Intensities), nCh);
                    var unmodSum = SumPerChannel(unmodPsms.Select(x => x.psm.Intensities), nCh);

                    // Occupancy per channel = modified / (modified + unmodified).
                    var occupancy = new double[nCh];
                    for (int i = 0; i < nCh; i++)
                    {
                        double total = modSum[i] + unmodSum[i];
                        occupancy[i] = total > 0 ? modSum[i] / total : 0.0;
                    }

                    var outString = new List<string> { site.mod, site.acc, site.pos.ToString(), gene, modPsms.Count.ToString(), unmodPsms.Count.ToString() };
                    outString.AddRange(modSum.Select(v => v.ToString(CultureInfo.InvariantCulture)));
                    outString.AddRange(unmodSum.Select(v => v.ToString(CultureInfo.InvariantCulture)));
                    outString.AddRange(occupancy.Select(v => v.ToString(CultureInfo.InvariantCulture)));
                    writer.WriteLine(string.Join("\t", outString));
                }
            }
        }

        [Test]
        public static void TransferIds()
        {
            string fullSequence = "AAGGAGAQVGGSISSGSSASSVTVTR";
            string converted = TmtPair.AddTmtLabelsToSequence(fullSequence);
            string fullSeq2 = "EVIITPNSAWGGEGSLGC[Common Fixed:Carbamidomethyl on C]GIGYGYLHR";
            string converted2 = TmtPair.AddTmtLabelsToSequence(fullSeq2);
            string fullSeq3 = "FGIVTSSAGTGTTEDTEAK[]";
            string converted3 = TmtPair.AddTmtLabelsToSequence(fullSeq3);
            string fullSeq4 = "HAEGTFTSDVSS[Common Biological:Phosphorylation on S]YLEGQAAK";
            string converted4 = TmtPair.AddTmtLabelsToSequence(fullSeq4);
            string removed1 = TmtPair.RemoveTmtLabels(converted);
            string removed2 = TmtPair.RemoveTmtLabels(converted2);
            string removed3 = TmtPair.RemoveTmtLabels(converted3);
            string removed4 = TmtPair.RemoveTmtLabels(converted4);

            var fullSeq = "[Multiplex Label:TMT18 on X]GQAGPEGAAP[Common Biological:Hydroxylation on P]APEEDK[Multiplex Label:TMT18 on K]";
            var parsed = SpectrumMatchFromTsv.ParseModifications(fullSeq);
            var notInteresting = new List<string> { "TMT18", "Fixed", "Artifact", "Variable" };
            var contains = parsed.Values.Any(v => !notInteresting.Any(key => v.Contains(key)));
        }


        public static string ModType(string fullSequence)
        {
            var mod = "Unmodified";
            if (!SpectrumMatchFromTsv.ParseModifications(fullSequence).Values.All(v => v.Contains("Fixed") || v.Contains("Variable") || v.Contains("TMT")))
            {
                mod = "GptmdMod";
                if (fullSequence.Contains("Biological") || fullSequence.Contains("Uniprot"))
                {
                    mod = "BioMod";
                    if (fullSequence.IndexOf("Phospho", StringComparison.OrdinalIgnoreCase) >= 0)
                    {
                        mod = "Phospho";
                    }
                }
            }
            return mod;
        }

        public static IEnumerable<PsmFromTsv> GetUnmodifiedPeptides(IEnumerable<PsmFromTsv> psms)
        {
            return psms.Where(p => SpectrumMatchFromTsv.ParseModifications(p.FullSequence).Values.All(v => v.Contains("Fixed") || v.Contains("Variable") || v.Contains("TMT") || v.Contains("sub")) && !p.FullSequence.Contains("|"));
        }

        public static IEnumerable<PsmFromTsv> GetModifiedPeptides(IEnumerable<PsmFromTsv> psms)
        {
            var notInteresting = new List<string> { "TMT18", "Fixed", "Artifact", "Variable", "Metal" };
            return psms.Where(p => SpectrumMatchFromTsv.ParseModifications(p.FullSequence).Values.Any(v => !notInteresting.Any(key => v.Contains(key))));
        }

        [Test]
        public static void TIC()
        {
            var directory = @"E:\Aneuploidy\Mistranslation_project\011626\040326_IS_test\";
            var files = Directory.GetFiles(directory, "*04-18*");
            var files2 = Directory.GetFiles(directory, "*04-19*");
            var fileList3 = files.Concat(files2).ToList();
            var fileLists2 = new List<string> { @"E:\Aneuploidy\Mistranslation_project\011626\040326_IS_test\04-16-26_TMT_mix.raw" ,
            @"E:\Aneuploidy\Mistranslation_project\011626\040326_IS_test\04-17-26_YL_TMT.raw"};
            var ticAreas = new Dictionary<string, string>();
            foreach (var file in fileList3)
            {
                var dataFile = MsDataFileReader.GetDataFile(file);
                var ms1TicArea = dataFile.GetMS1Scans().Where(s => s.RetentionTime >= 15 && s.RetentionTime <= 110).Sum(s => s.TotalIonCurrent);
                ticAreas[file] = $"{ms1TicArea:E2}";
            }
        }
    }

    public class TmtPair
    {
        public PsmFromTsv BestPsm_LF { get; set; }
        public List<PsmFromTsv> AllPsmsTMT { get; set; }
        public double[] ReporterIonIntensities { get; set; }

        public TmtPair(PsmFromTsv bestPsm_LF, IEnumerable<PsmFromTsv> psms = null)
        {
            BestPsm_LF = bestPsm_LF;
            AllPsmsTMT = psms == null ? new List<PsmFromTsv>() : psms.ToList();
        }

        public void FilterPsms(int fractionTol, double rtTol, double similarityTol)
        {
            string fraction_LF = Regex.Match(BestPsm_LF.FileNameWithoutExtension, @"Frxn\d+").Value;
            for (int i = AllPsmsTMT.Count - 1; i >= 0; i--)
            {
                var psm = AllPsmsTMT[i];
                var fraction_tmt = Regex.Match(psm.FileNameWithoutExtension, @"Frxn\d+").Value;
                if (Math.Abs(int.Parse(fraction_LF.Substring(4)) - int.Parse(fraction_tmt.Substring(4))) > fractionTol)
                {
                    AllPsmsTMT.RemoveAt(i); continue;
                }

                //RT filter
                var predictedRt = 0;
                if (Math.Abs(predictedRt - psm.RetentionTime) > rtTol)
                {
                    AllPsmsTMT.RemoveAt(i); continue;
                }

                //predicted spectrum
                var spectralSimilarity = 0;
                if (spectralSimilarity < similarityTol)
                {
                    AllPsmsTMT.RemoveAt(i); continue;
                }
            }
        }

        public static void WriteResults(IEnumerable<PsmFromTsv> allPsms, string outPath)
        {
            var modsToexclude = new List<string> { "TMT18", "Fixed", "Artifact", "Variable", "Metal" };
            var columns = new List<string> { "BaseSequence", "FullSequence", "Mods", "Protein Accession", "Protein Name", "GeneName", "Description" };
            var labels = new List<string> { "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C", "132N", "132C", "133N", "133C", "134N", "134C", "135N" };
            columns.AddRange(labels);
            using (StreamWriter writer = new StreamWriter(outPath))
            {
                writer.WriteLine(string.Join("\t", columns));
                var peptides = allPsms.GroupBy(p => p.FullSequence);
                foreach (var peptide in peptides)
                {
                    var mods = String.Join(", ", SpectrumMatchFromTsv.ParseModifications(peptide.First().FullSequence).Values.Where(m => !modsToexclude.Any(x => m.Contains(x))));
                    var outString = new List<string> { peptide.First().BaseSequence, peptide.First().FullSequence, mods, peptide.First().ProteinAccession, peptide.First().ProteinName, peptide.First().GeneName, peptide.First().Description };

                    var reporterIonIntensities = new double[peptide.First().Intensities.Count()];
                    foreach (var psm in peptide)
                    {
                        for (int i = 0; i < peptide.First().Intensities.Count(); i++)
                        {
                            reporterIonIntensities[i] += psm.Intensities[i];
                        }
                    }
                    outString.AddRange(reporterIonIntensities.Select(i => i.ToString()));
                    writer.WriteLine(string.Join("\t", outString));
                }
            }
        }

        public static string AddTmtLabelsToSequence(string fullSequence)
        {
            var sb = new StringBuilder();

            if (fullSequence.StartsWith("["))
            {
                // Already modified N-term
                int endBracket = fullSequence.IndexOf(']');
                sb.Append(fullSequence.Substring(0, endBracket + 1));
                fullSequence = fullSequence.Substring(endBracket + 1);
            }
            else
            {
                // Add TMT to N-term
                sb.Append("[Multiplex Label:TMT18 on X]");
            }

            // Process the rest
            for (int i = 0; i < fullSequence.Length; i++)
            {
                char c = fullSequence[i];
                sb.Append(c);

                if (c == 'K')
                {
                    // Check if next char is a modification
                    if (i + 1 < fullSequence.Length && fullSequence[i + 1] == '[')
                    {
                        // Already modified, do nothing
                        continue;
                    }
                    else
                    {
                        sb.Append("[Multiplex Label:TMT18 on K]");
                    }
                }
            }

            return sb.ToString();
        }

        public static string RemoveTmtLabels(string fullSequence)
        {
            // Remove N-term TMT label if present
            fullSequence = Regex.Replace(
                fullSequence,
                @"^\[Multiplex Label:TMT18 on X\]",
                "",
                RegexOptions.Compiled);

            // Remove TMT label on K (after K)
            fullSequence = Regex.Replace(
                fullSequence,
                @"K\[Multiplex Label:TMT18 on K\]",
                "K",
                RegexOptions.Compiled);

            return fullSequence;
        }

        
    }
}

