using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Readers;
using NUnit.Framework;
using System.Data;
using System.Net.WebSockets;
using System.IO;
using OxyPlot;
using System.Collections;
using System.Globalization;
using System.Windows.Documents;
using System.Windows.Input;
using System.Xml.Linq;
using PredictionClients.Koina.SupportedModels.RetentionTimeModels;
using PredictionClients.Koina.SupportedModels.FragmentIntensityModels;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.Util;
using MassSpectrometry;
using MassSpectrometry.MzSpectra;
using MzLibUtil;
using Omics;
using Omics.SpectrumMatch;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using Chemistry;
using Readers.SpectralLibrary;
using System.Data.Entity.Core.Metadata.Edm;
using Chromatography.RetentionTimePrediction;
using Chromatography.RetentionTimePrediction.Chronologer;
using Omics.Modifications;
using OxyPlot.Reporting;
using IncompatibleModHandlingMode = Chromatography.RetentionTimePrediction.IncompatibleModHandlingMode;

namespace Test
{

    internal class AAsub
    {
        [Test]
        public static void AAsubQuant()
        {
            var dir = @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\UPLC\PAW\UPLC_all\msn_files";
            var files = System.IO.Directory.GetFiles(dir, "*PAW_tmt.txt");
            var psmPath = @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\UPLC\LFgptmdPrunedDb-AAsub_search-cali-search\Task1-SearchTask\AllPSMs.psmtsv";
            var psmFile = new PsmFromTsvFile(psmPath);
            var allPsms = psmFile.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T" && p.AmbiguityLevel == "1");

            var columns = new List<string> { "BaseSequence", "FullSequence", "Mods", "Protein Accession", "Protein Name", "GeneName", "Description", "Start_and_End_Residues" };
            var labels = new List<string> { "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C", "132N", "132C", "133N", "133C", "134N", "134C", "135N" };
            columns.AddRange(labels);
            var notInteresting = new List<string> { "TMT18", "Fixed", "Artifact", "Variable", "Metal" };
            var outDir = @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\UPLC\LFgptmdPrunedDb-AAsub_search-cali-search";

            var exps = new List<string> { "R2", "R4", "R6" };
            foreach (var exp in exps)
            {
                var expPsms = allPsms.Where(p => p.FileName.Contains(exp));
                var psmsGroupedByFile = expPsms.GroupBy(p => p.FileName);

                foreach (var psmGroup in psmsGroupedByFile)
                {
                    var reporterFile = files.FirstOrDefault(f => Path.GetFileName(f).Contains(psmGroup.Key));
                    var lines = File.ReadAllLines(reporterFile).Skip(1);
                    var reporterDic = lines.Select(line => line.Split('\t')).ToDictionary(parts => int.Parse(parts[1]), parts => parts.Skip(parts.Length - 18)
                        .Select(s => double.Parse(s, CultureInfo.InvariantCulture)).ToArray());
                    foreach(var psm in psmGroup)
                    {
                        psm.Intensities = reporterDic[psm.Ms2ScanNumber];
                    }
                }

                var outPath = Path.Combine(outDir, $"{exp}_AAsub_quant.tsv");
                using (StreamWriter writer = new StreamWriter(outPath))
                {
                    writer.WriteLine(string.Join("\t", columns));
                    var peptides = expPsms.GroupBy(p => p.FullSequence);
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
        }

        [Test]
        public static void WriteOutReporterIons()
        {
            var dir = @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\UPLC\PAW\UPLC_all\msn_files";
            var files = System.IO.Directory.GetFiles(dir, "*PAW_tmt.txt");
            var psmPath = @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\UPLC\LFgptmdPrunedDb-AAsub_search-cali-search\Task1-SearchTask\AllPSMs.psmtsv";
            var psmFile = new PsmFromTsvFile(psmPath);
            var allPsms = psmFile.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T" && p.AmbiguityLevel == "1");

            var columns = new List<string> { "FileName", "ScanNumber", "BaseSequence", "FullSequence", "Accession" };
            var labels = new List<string> { "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C", "132N", "132C", "133N", "133C", "134N", "134C", "135N" };
            columns.AddRange(labels);

            var outPath = @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\UPLC\LFgptmdPrunedDb-AAsub_search-cali-search\FilteredPsmsWithReporterIonIntensities.tsv";
            var psmsGroupedByFile = allPsms.GroupBy(p => p.FileName);
            foreach (var psmGroup in psmsGroupedByFile)
            {
                var reporterFile = files.FirstOrDefault(f => Path.GetFileName(f).Contains(psmGroup.Key));
                var lines = File.ReadAllLines(reporterFile).Skip(1);
                var reporterDic = lines.Select(line => line.Split('\t')).ToDictionary(parts => int.Parse(parts[1]), parts => parts.Skip(parts.Length - 18)
                    .Select(s => double.Parse(s, CultureInfo.InvariantCulture)).ToArray());
                foreach (var psm in psmGroup)
                {
                    psm.Intensities = reporterDic[psm.Ms2ScanNumber];
                }
            }

            var exps = new List<string> { "R2", "R4", "R6" };
            var outDir = @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\UPLC\LFgptmdPrunedDb-AAsub_search-cali-search";
            foreach (var exp in exps)
            {
                var expOutPath = Path.Combine(outDir, $"AllPSMs_withRI_AAsub_{exp}.tsv");
                var expPsms = allPsms.Where(p => p.FileName.Contains(exp));
                using (StreamWriter writer = new StreamWriter(expOutPath))
                {
                    writer.WriteLine(string.Join("\t", columns));
                    foreach (var psm in expPsms)
                    {
                        var outString = new List<string> { psm.FileName, psm.Ms2ScanNumber.ToString(), psm.BaseSeq, psm.FullSequence, psm.Accession };
                        outString.AddRange(psm.Intensities.Select(i => i.ToString()));
                        writer.WriteLine(string.Join("\t", outString));
                    }
                }
            }
        }

        /// <summary>
        /// Reads AllPeptides.psmtsv (MetaMorpheus), predicts retention time with Prosit2019iRT,
        /// and writes two output TSV files next to the input:
        ///   1) AllPeptides_RT_predictions_all.tsv   — every peptide that yields a prediction
        ///   2) AllPeptides_RT_predictions_filtered.tsv — same, but excluding peptides whose
        ///      Mass Diff (Da) has |shift| &lt; 2 Da (too low to be a real AA sub) or is within
        ///      2 Da of common PTM masses (e.g. oxidation 15.9949 Da).
        /// </summary>
        [Test]
        public static void PredictRtForAllPeptides()
        {
            var psmPath = @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\LF\raw_gptmd-2nd-AAsub_writePrunedDb\Task2-SearchTask\AllPeptides.psmtsv";
            var outDir = @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\LF\raw_gptmd-2nd-AAsub_writePrunedDb";
            var allOutPath = Path.Combine(outDir, "AllPeptides_RT_predictions_all.tsv");
            var filteredOutPath = Path.Combine(outDir, "AllPeptides_RT_predictions_filtered.tsv");

            var psmFile = new PsmFromTsvFile(psmPath);
            var peptides = psmFile.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T" && p.AmbiguityLevel == "1").ToList();
            var unmodPeptides = Ptm_tmt.GetUnmodifiedPeptides(peptides).ToList();

            // Build the sequence we feed to Prosit:
            //  - peptides carrying a "nucleotide substitution" mod get the substitution
            //    applied so the model sees the translated sequence;
            //  - otherwise the FullSequence is passed through unchanged.
            var sequencesForPrediction = unmodPeptides
                .Select(p => p.FullSequence.Contains("sub") ? DeepLC.ApplyAaSubstitutions(p.FullSequence) : p.FullSequence)
                .ToList();

            // Predict RT with Prosit2019iRT. RemoveIncompatibleMods strips mods the model
            // can't handle (anything other than Met-ox and Carbamidomethyl on C) and still
            // returns a prediction on the base sequence.
            var modelInputs = sequencesForPrediction.Select(s => new RetentionTimePredictionInput(s)).ToList();
            var model = new Prosit2019iRT();
            var predictions = model.Predict(modelInputs);

            using (var allWriter = new StreamWriter(allOutPath))
            using (var filtWriter = new StreamWriter(filteredOutPath))
            {
                var header = string.Join("\t", "BaseSeq", "FullSequence", "SubstitutedSequence", "ObservedRt", "PredictedRt", "Score", "AAsub", "AAsubMassShifts");
                allWriter.WriteLine(header);
                filtWriter.WriteLine(header);

                int allCount = 0, filtCount = 0;
                for (int i = 0; i < unmodPeptides.Count; i++)
                {
                    var p = unmodPeptides[i];
                    var pred = predictions[i];
                    if (pred.PredictedRetentionTime == null) continue; // no prediction available

                    // Reconstruct the AA-sub annotations (semicolon-joined) for traceability.
                    var aasub = SpectrumMatchFromTsv.ParseModifications(p.FullSequence)
                        .Values.Where(v => v.Contains("nucleotide substitution"));
                    var aasubString = aasub.Any() ? string.Join(";", aasub) : "";

                    // Per-sub annotation in the form "X -> Y: ±N.NNNN Da", same order as the
                    // AAsub annotations; multiple subs joined with ';'.
                    var annotations = DeepLC.GetAaSubstitutionAnnotations(p.FullSequence);
                    var shiftsString = annotations.Count > 0
                        ? string.Join(";", annotations)
                        : "";

                    var row = string.Join("\t",
                        p.BaseSeq,
                        p.FullSequence,
                        sequencesForPrediction[i],
                        p.RetentionTime.ToString(CultureInfo.InvariantCulture),
                        pred.PredictedRetentionTime.Value.ToString(CultureInfo.InvariantCulture),
                        p.Score.ToString(CultureInfo.InvariantCulture),
                        aasubString,
                        shiftsString);

                    allWriter.WriteLine(row);
                    allCount++;

                    // Filter using the per-substitution residue-mass delta (mass(to) - mass(from)).
                    // Drops peptides whose any sub has |delta| < 2 Da (e.g. I/L, Q/K, N/D)
                    // or sits within 0.05 Da of a common PTM mass (oxidation 15.9949, etc.).
                    if (!DeepLC.HasPtmConfoundedAaSubstitution(p.FullSequence))
                    {
                        filtWriter.WriteLine(row);
                        filtCount++;
                    }
                }

                TestContext.WriteLine($"Wrote {allCount} peptides to {allOutPath}");
                TestContext.WriteLine($"Wrote {filtCount} filtered peptides to {filteredOutPath}");
            }
        }

        /// <summary>
        /// For every identified peptide in AllPeptides.psmtsv, predict its fragment-ion spectrum
        /// with Prosit2020IntensityHCD (label-free HCD model) and compute the cosine spectral
        /// similarity against the matching MS2 scan in the raw file. Substitutions are applied
        /// to the FullSequence before prediction so the model scores the translated peptide.
        /// Output TSV columns: FileName, ScanNumber, BaseSeq, FullSequence, SubstitutedSequence,
        /// Charge, AAsub, AAsubMassShifts, CosineSimilarity.
        /// </summary>
        [Test]
        public static void SpectralSimilarityForAaSubPeptides()
        {
            var psmPath = @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\LF\raw_2nd-AAsub-gptmdDb_writePrunedDb_oldBranch\Task1-SearchTask\AllPeptides.psmtsv";
            var rawDir = @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\LF";
            var outPath = @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\LF\raw_2nd-AAsub-gptmdDb_writePrunedDb_oldBranch\AllPeptides_SpectralSimilarity2.tsv";

            // 1. Load PSMs. Keep all q-value/target/unambiguous peptides regardless of mod
            //    state so Fixed/Variable mod (Carbamidomethyl C, Oxidation M) AND BioMod
            //    peptides are scored — Prosit2020IntensityHCD with RemoveIncompatibleMods
            //    strips mods it can't predict (e.g. Phospho) and still returns a spectrum.
            var psmFile = new PsmFromTsvFile(psmPath);
            var peptides = psmFile.Results
                .Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T"
                            && p.AmbiguityLevel == "1" && !p.FullSequence.Contains("|"))
                .ToList();

            // 2. Apply AA substitutions so the model sees the translated sequence.
            //    Each PSM gets (substituted sequence, charge) — that's the prediction key.
            var psmPredictionKeys = peptides
                .Select(p => (
                    seq:    p.FullSequence.Contains("sub") ? DeepLC.ApplyAaSubstitutions(p.FullSequence) : p.FullSequence,
                    charge: p.ChargeState))
                .ToList();

            // 3. Predict only once per unique (sequence, charge).
            var uniqueKeys = psmPredictionKeys.Distinct().ToList();
            var modelInputs = uniqueKeys
                .Select(k => new FragmentIntensityPredictionInput(
                    FullSequence: k.seq,
                    PrecursorCharge: k.charge,
                    CollisionEnergy: 30,
                    InstrumentType: null,
                    FragmentationType: null))
                .ToList();

            // Fully qualified because the file aliases `IncompatibleModHandlingMode`
            // to the Chronologer enum (line 35 using-alias). Prosit needs the Koina one.
            var model = new Prosit2020IntensityHCD(
                modHandlingMode: PredictionClients.Koina.Util.IncompatibleModHandlingMode.RemoveIncompatibleMods,
                fragmentIonMappingMode: FragmentIonMappingMode.MapToValidatedFullSequence);
            model.Predict(modelInputs);

            // 4. Build predicted (m/z, intensity) arrays per (sequence, charge). We DON'T use
            //    GenerateLibrarySpectraFromPredictions because its bulk ToDictionary on
            //    Product.Annotation can crash with duplicate keys (e.g. "M0-64.00") and kill
            //    the whole batch. Replicate its logic per-peptide with dedup + try/catch.
            double minIntensityFilter = 1e-6;
            var keyToPredicted = new Dictionary<(string, int), (double[] mz, double[] intensity)>();
            for (int i = 0; i < uniqueKeys.Count; i++)
            {
                if (!model.ValidInputsMask[i]) continue;
                var prediction = model.Predictions[i];
                if (prediction == null) continue;

                try
                {
                    string fullSeq = model.FragmentIonMappingMode == FragmentIonMappingMode.MapToValidatedFullSequence
                        ? prediction.ValidatedFullSequence
                        : prediction.FullSequence;
                    var peptide = new PeptideWithSetModifications(fullSeq);
                    var theoretical = new List<Product>();
                    peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, theoretical);
                    // GroupBy().First() instead of ToDictionary to tolerate duplicate annotations.
                    var tpLookup = theoretical
                        .GroupBy(tp => tp.Annotation)
                        .ToDictionary(g => g.Key, g => g.First());

                    var mzList = new List<double>();
                    var intList = new List<double>();
                    for (int j = 0; j < prediction.FragmentAnnotations.Count; j++)
                    {
                        var ann = prediction.FragmentAnnotations[j];
                        var intens = prediction.FragmentIntensities[j];
                        if (intens == -1 || intens < minIntensityFilter || ann == null || !ann.Contains("+")) continue;
                        var parts = ann.Split("+");
                        if (!tpLookup.TryGetValue(parts[0], out var tp)) continue;
                        int charge = int.Parse(parts[1]);
                        mzList.Add(tp.ToMz(charge));
                        intList.Add(intens);
                    }

                    if (mzList.Count == 0) continue;
                    // Sort by m/z (SpectralSimilarity expects ascending).
                    var idx = Enumerable.Range(0, mzList.Count).OrderBy(k => mzList[k]).ToArray();
                    var mzSorted = idx.Select(k => mzList[k]).ToArray();
                    var intSorted = idx.Select(k => intList[k]).ToArray();
                    keyToPredicted[uniqueKeys[i]] = (mzSorted, intSorted);
                }
                catch
                {
                    // Skip peptides that can't be fragmented cleanly.
                    continue;
                }
            }

            // 5. Cache raw files by base filename.
            var rawFilesDict = Directory.GetFiles(rawDir, "*.raw")
                .ToDictionary(Path.GetFileNameWithoutExtension, MsDataFileReader.GetDataFile);

            // 6. Write results.
            var productTol = new PpmTolerance(20);
            using (var writer = new StreamWriter(outPath))
            {
                writer.WriteLine(string.Join("\t",
                    "FileName", "ScanNumber", "BaseSeq", "FullSequence", "SubstitutedSequence",
                    "Charge", "AAsub", "AAsubMassShifts", "PTM", "CosineSimilarity"));

                int wrote = 0, skippedNoSpec = 0, skippedNoRaw = 0;
                for (int i = 0; i < peptides.Count; i++)
                {
                    var p = peptides[i];
                    var key = psmPredictionKeys[i];

                    if (!keyToPredicted.TryGetValue(key, out var pred))
                    { skippedNoSpec++; continue; }

                    if (!rawFilesDict.TryGetValue(p.FileName, out var rawFile))
                    { skippedNoRaw++; continue; }
                    var rawScan = rawFile.GetOneBasedScan(p.Ms2ScanNumber);
                    if (rawScan == null) { skippedNoRaw++; continue; }

                    // Use the 2-MzSpectrum constructor — the 4-array overload has a bug
                    // where LocalTolerance is never initialized, leading to a NullRef in
                    // GetIntensityPairs. Wrap our predicted arrays in an MzSpectrum.
                    var predictedSpectrum = new MzSpectrum(pred.mz, pred.intensity, shouldCopy: false);
                    var similarity = new SpectralSimilarity(
                        rawScan.MassSpectrum,
                        predictedSpectrum,
                        SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum,
                        toleranceInPpm: 20,
                        allPeaks: false,
                        tol: productTol);
                    double? cos = similarity.CosineSimilarity();

                    var aasub = SpectrumMatchFromTsv.ParseModifications(p.FullSequence)
                        .Values.Where(v => v.Contains("nucleotide substitution"));
                    var aasubString = aasub.Any() ? string.Join(";", aasub) : "";
                    var annotations = DeepLC.GetAaSubstitutionAnnotations(p.FullSequence);
                    var shiftsString = annotations.Count > 0 ? string.Join(";", annotations) : "";
                    // PTM label from FullSequence so BioMod/GptmdMod peptides are flagged
                    // (after AA-sub annotations are accounted for separately).
                    var ptmLabel = Ptm_tmt.ModType(p.FullSequence);

                    writer.WriteLine(string.Join("\t",
                        p.FileName,
                        p.Ms2ScanNumber.ToString(CultureInfo.InvariantCulture),
                        p.BaseSeq,
                        p.FullSequence,
                        key.seq,
                        p.ChargeState.ToString(CultureInfo.InvariantCulture),
                        aasubString,
                        shiftsString,
                        ptmLabel,
                        cos?.ToString("F4", CultureInfo.InvariantCulture) ?? ""));
                    wrote++;
                }

                TestContext.WriteLine($"Wrote {wrote} rows to {outPath}");
                TestContext.WriteLine($"Skipped {skippedNoSpec} (no predicted spectrum) and {skippedNoRaw} (no raw scan).");
            }
        }

        /// <summary>
        /// Same flow as PredictRtForAllPeptides but uses the Chronologer retention-time
        /// predictor (Chromatography.RetentionTimePrediction.Chronologer) instead of
        /// Prosit2019iRT. AA substitutions are applied before prediction so Chronologer
        /// scores the translated peptide. Writes two TSVs next to the input:
        ///   1) AllPeptides_RT_predictions_chronologer_all.tsv
        ///   2) AllPeptides_RT_predictions_chronologer_filtered.tsv (drops PTM-confounded subs)
        /// </summary>
        [Test]
        public static void PredictRtForAllPeptidesChronologer()
        {
            var psmPath = @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\LF\raw_2nd-AAsub-gptmdDb_writePrunedDb_oldBranch\Task1-SearchTask\AllPeptides.psmtsv";
            var outDir  = @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\LF\raw_2nd-AAsub-gptmdDb_writePrunedDb_oldBranch";
            var allOutPath      = Path.Combine(outDir, "AllPeptides_RT_predictions_chronologer_all2.tsv");
            var filteredOutPath = Path.Combine(outDir, "AllPeptides_RT_predictions_chronologer_filtered2.tsv");

            var psmFile = new PsmFromTsvFile(psmPath);
            var peptides = psmFile.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T" && p.AmbiguityLevel == "1").ToList();
            //var unmodPeptides = Ptm_tmt.GetUnmodifiedPeptides(peptides).ToList();

            using var predictor = new ChronologerRetentionTimePredictor(IncompatibleModHandlingMode.RemoveIncompatibleMods);

            using var allWriter  = new StreamWriter(allOutPath);
            using var filtWriter = new StreamWriter(filteredOutPath);
            var header = string.Join("\t",
                "BaseSeq", "FullSequence", "SubstitutedSequence",
                "ObservedRt", "PredictedRt", "Score", "AAsub", "AAsubMassShifts", "PTM");
            allWriter.WriteLine(header);
            filtWriter.WriteLine(header);

            int allCount = 0, filtCount = 0, skipped = 0;
            foreach(var peptide in peptides)
            {
                double? predicted = null;
                var baseSeq = peptide.BaseSeq;
                var fullSeq = peptide.FullSequence;
                if (fullSeq.Contains("nucleotide substitution"))
                {
                    fullSeq = DeepLC.ApplyAaSubstitutions(peptide.FullSequence); // only test AAsub peptides
                    baseSeq = IBioPolymerWithSetMods.GetBaseSequenceFromFullSequence(fullSeq);
                }
                var allModsDic = SpectrumMatchFromTsv.ParseModifications(fullSeq);
                try
                {
                    // Pass null (not an empty dict) so PWSM falls back to
                    // Mods.AllKnownProteinModsDictionary. IBioPolymerWithSetMods.
                    // GetModificationDictionaryFromFullSequence strips the
                    // "Common Fixed:" / "Common Variable:" / "Common Biological:"
                    // prefix before lookup, so all standard MM annotations resolve.
                    var pep = new PeptideWithSetModifications(fullSeq, null);
                    predicted = predictor.PredictRetentionTime(pep, out var failure);
                }
                catch
                {
                    skipped++;
                    continue;
                }

                if (predicted == null) { skipped++; continue; }

                var aasub = SpectrumMatchFromTsv.ParseModifications(peptide.FullSequence).Values.Where(v => v.Contains("nucleotide substitution"));
                var aasubString = aasub.Any() ? string.Join(";", aasub) : "";
                var annotations = DeepLC.GetAaSubstitutionAnnotations(peptide.FullSequence);
                var shiftsString = annotations.Count > 0 ? string.Join(";", annotations) : "";
                var mods = Ptm_tmt.ModType(fullSeq);

                var row = string.Join("\t", peptide.BaseSeq, peptide.FullSequence, fullSeq, peptide.RetentionTime.ToString(CultureInfo.InvariantCulture),
                    predicted.Value.ToString(CultureInfo.InvariantCulture), peptide.Score.ToString(CultureInfo.InvariantCulture), aasubString, shiftsString, mods);

                allWriter.WriteLine(row);
                allCount++;

                // Filter: drop PTM-confounded subs (same logic as the Prosit test).
                if (!DeepLC.HasPtmConfoundedAaSubstitution(peptide.FullSequence))
                {
                    filtWriter.WriteLine(row);
                    filtCount++;
                }
            }
        }

        [Test]
        public static void TestMs2PIPinputAAsub()
        {
            var psmPath = @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\LF\raw_2nd-AAsub-gptmdDb_writePrunedDb_oldBranch\Task1-SearchTask\AllPeptides.psmtsv";
            var psmFile = new PsmFromTsvFile(psmPath);
            var peptides = psmFile.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T" && p.AmbiguityLevel == "1").ToList();
            //var unmodPeptides = Ptm_tmt.GetUnmodifiedPeptides(peptides).ToList();

            var outPath = @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\LF\raw_2nd-AAsub-gptmdDb_writePrunedDb_oldBranch\allPeptides_input_3_ms2pip.peprec";
            using (StreamWriter writer = new StreamWriter(outPath))
            {
                writer.WriteLine("spec_id modifications peptide charge");
                foreach (var peptide in peptides)
                {
                    var baseSeq = peptide.BaseSeq;
                    var fullSeq = peptide.FullSequence;
                    if (fullSeq.Contains("nucleotide substitution"))
                    {
                        fullSeq = DeepLC.ApplyAaSubstitutions(peptide.FullSequence); // only test AAsub peptides
                        baseSeq = IBioPolymerWithSetMods.GetBaseSequenceFromFullSequence(fullSeq);
                    }
                    var parsedMods = DeepLC.ParseModsForDeepLC(fullSeq);
                    var outString = new List<string> { $"scan{peptide.Ms2ScanNumber}", parsedMods, baseSeq, peptide.ChargeState.ToString() };
                    writer.WriteLine(string.Join(" ", outString));
                }
            }
        }

        [Test]
        public static void TestSpectralSimilarityAAsub()
        {
            var psmPath = @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\LF\raw_2nd-AAsub-gptmdDb_writePrunedDb_oldBranch\Task1-SearchTask\AllPeptides.psmtsv";
            var psmFile = new PsmFromTsvFile(psmPath);
            var peptides = psmFile.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T" && p.AmbiguityLevel == "1").ToList();
            var unmodPeptides = Ptm_tmt.GetUnmodifiedPeptides(peptides).ToList();

            var raw_directory = @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\LF";
            var raw_paths = Directory.GetFiles(raw_directory, "*.raw");//, SearchOption.AllDirectories
            var rawFilesDict = raw_paths.ToDictionary(r => r, r => MsDataFileReader.GetDataFile(r));

            var libraryPath = @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\LF\raw_2nd-AAsub-gptmdDb_writePrunedDb_oldBranch\allPeptides_2_ms2pip.msp";
            var spectralLibrary = new SpectralLibrary(new List<string> { libraryPath });
            var allLibrarySepctra = spectralLibrary.Results;
            var similarities = new List<double>();

            var similarityOutPath = @"E:\Aneuploidy\Mistranslation_project\011626\042426_TMT\LF\raw_2nd-AAsub-gptmdDb_writePrunedDb_oldBranch\SpectralSimilarity_allPeptides_3.tsv";
            using (StreamWriter writer = new StreamWriter(similarityOutPath))
            {
                writer.WriteLine(string.Join("\t", "ScanNumber", "BaseSeq", "FullSequence", "SubstitutedSequence", "AAsub", "PTM", "massshift", "CosineSimilarity"));

                var precursorTol = new PpmTolerance(20);
                var productTol = new PpmTolerance(20);
                foreach (var peptide in peptides)
                {
                    var rawFile = rawFilesDict.FirstOrDefault(kvp => kvp.Key.Contains(peptide.FileName)).Value;
                    var rawScan = rawFile.GetOneBasedScan(peptide.Ms2ScanNumber);

                    var baseSeq = peptide.BaseSeq;
                    var fullSeq = peptide.FullSequence;
                    if (fullSeq.Contains("nucleotide substitution"))
                    {
                        fullSeq = DeepLC.ApplyAaSubstitutions(peptide.FullSequence); // only test AAsub peptides
                        baseSeq = IBioPolymerWithSetMods.GetBaseSequenceFromFullSequence(fullSeq);
                    }

                    var librarySpectrum = spectralLibrary.Results.FirstOrDefault(s => IBioPolymerWithSetMods.GetBaseSequenceFromFullSequence(s.Sequence) == baseSeq && s.ChargeState == peptide.ChargeState && precursorTol.Within(s.PrecursorMz.ToMass(s.ChargeState), peptide.MonoisotopicMass));//Math.Round(s.PrecursorMz, 0) == Math.Round(peptide.PrecursorMass.ToMz(peptide.ChargeState), 0)
                    if (librarySpectrum != null)
                    {
                        var aasub = SpectrumMatchFromTsv.ParseModifications(peptide.FullSequence).Values.Where(v => v.Contains("nucleotide substitution"));
                        var aasubString = aasub.Any() ? string.Join(";", aasub) : "";
                        var annotations = DeepLC.GetAaSubstitutionAnnotations(peptide.FullSequence);
                        var shiftsString = annotations.Count > 0 ? string.Join(";", annotations) : "";
                        var ptm = Ptm_tmt.ModType(fullSeq);

                        var similarity = new SpectralSimilarity(rawScan.MassSpectrum, librarySpectrum, SpectralSimilarity.SpectrumNormalizationScheme.SquareRootSpectrumSum, 20, false, tol: productTol);
                        similarities.Add(similarity.CosineSimilarity().Value);
                        var outString = new List<string> { peptide.Ms2ScanNumber.ToString(), peptide.BaseSeq, peptide.FullSequence, fullSeq, aasubString, ptm, shiftsString, similarity.CosineSimilarity().Value.ToString() };
                        writer.WriteLine(string.Join("\t", outString));
                    }
                }
            }
        }

    }

}
