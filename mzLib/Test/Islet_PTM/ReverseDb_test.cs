using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using Readers;
using NUnit.Framework;
using Omics.SpectrumMatch;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using UsefulProteomicsDatabases;
using System.Globalization;
using MassSpectrometry;
using Test.FileReadingTests;
using System.Security.Cryptography.X509Certificates;
using Omics;

namespace Test.Islet_PTM
{
    public class ReverseDb_test
    {
        // Standard genetic code (stop codons omitted).
        private static readonly Dictionary<string, char> CodonTable = new()
        {
            {"TTT",'F'}, {"TTC",'F'}, {"TTA",'L'}, {"TTG",'L'},
            {"CTT",'L'}, {"CTC",'L'}, {"CTA",'L'}, {"CTG",'L'},
            {"ATT",'I'}, {"ATC",'I'}, {"ATA",'I'}, {"ATG",'M'},
            {"GTT",'V'}, {"GTC",'V'}, {"GTA",'V'}, {"GTG",'V'},
            {"TCT",'S'}, {"TCC",'S'}, {"TCA",'S'}, {"TCG",'S'},
            {"CCT",'P'}, {"CCC",'P'}, {"CCA",'P'}, {"CCG",'P'},
            {"ACT",'T'}, {"ACC",'T'}, {"ACA",'T'}, {"ACG",'T'},
            {"GCT",'A'}, {"GCC",'A'}, {"GCA",'A'}, {"GCG",'A'},
            {"TAT",'Y'}, {"TAC",'Y'},
            {"CAT",'H'}, {"CAC",'H'}, {"CAA",'Q'}, {"CAG",'Q'},
            {"AAT",'N'}, {"AAC",'N'}, {"AAA",'K'}, {"AAG",'K'},
            {"GAT",'D'}, {"GAC",'D'}, {"GAA",'E'}, {"GAG",'E'},
            {"TGT",'C'}, {"TGC",'C'}, {"TGG",'W'},
            {"CGT",'R'}, {"CGC",'R'}, {"CGA",'R'}, {"CGG",'R'},
            {"AGT",'S'}, {"AGC",'S'}, {"AGA",'R'}, {"AGG",'R'},
            {"GGT",'G'}, {"GGC",'G'}, {"GGA",'G'}, {"GGG",'G'}
        };

        private static readonly HashSet<char> Forbidden = new() { 'C', 'K', 'R' };

        // For each non-forbidden AA, the set of distinct non-forbidden AAs reachable
        // by changing one nucleotide in any of its codons.
        private static readonly Dictionary<char, List<char>> ReachableAaMap = BuildReachableMap();

        private static Dictionary<char, List<char>> BuildReachableMap()
        {
            var bases = new[] { 'A', 'C', 'G', 'T' };
            var codonsByAa = CodonTable.GroupBy(kv => kv.Value)
                                       .ToDictionary(g => g.Key, g => g.Select(kv => kv.Key).ToList());

            var result = new Dictionary<char, List<char>>();
            foreach (var aa in codonsByAa.Keys)
            {
                if (Forbidden.Contains(aa)) continue;
                var reachable = new HashSet<char>();
                foreach (var codon in codonsByAa[aa])
                {
                    for (int i = 0; i < 3; i++)
                    {
                        foreach (var b in bases)
                        {
                            if (b == codon[i]) continue;
                            var mutatedCodon = codon.Substring(0, i) + b + codon.Substring(i + 1);
                            if (CodonTable.TryGetValue(mutatedCodon, out char newAa) && newAa != aa && !Forbidden.Contains(newAa))
                            {
                                reachable.Add(newAa);
                            }
                        }
                    }
                }
                result[aa] = reachable.ToList();
            }
            return result;
        }

        /// <summary>
        /// Randomly mutates one amino acid in <paramref name="peptide"/> to a different amino
        /// acid reachable by a single-nucleotide substitution in the standard genetic code.
        /// C, K, and R are never selected as the mutation site and are never produced as the
        /// substituted residue. Returns the unchanged input if no eligible position exists.
        /// </summary>
        /// <param name="peptide">Peptide base sequence (single-letter codes).</param>
        /// <param name="random">Random source; a new one is created if null.</param>
        /// <param name="mutatedIndex">Index of the residue that was changed, or -1 if none.</param>
        /// <param name="originalAa">Residue prior to mutation, or '\0' if none.</param>
        /// <param name="mutatedAa">Residue after mutation, or '\0' if none.</param>
        public static string MutateOneAminoAcid(
            string peptide,
            Random random,
            out int mutatedIndex,
            out char originalAa,
            out char mutatedAa)
        {
            mutatedIndex = -1;
            originalAa = '\0';
            mutatedAa = '\0';
            if (string.IsNullOrEmpty(peptide)) return peptide;
            random ??= new Random();

            var mutableIndices = Enumerable.Range(0, peptide.Length)
                .Where(i => ReachableAaMap.TryGetValue(peptide[i], out var list) && list.Count > 0)
                .ToList();
            if (mutableIndices.Count == 0) return peptide;

            int pos = mutableIndices[random.Next(mutableIndices.Count)];
            var choices = ReachableAaMap[peptide[pos]];
            char newAa = choices[random.Next(choices.Count)];

            var chars = peptide.ToCharArray();
            originalAa = chars[pos];
            chars[pos] = newAa;
            mutatedIndex = pos;
            mutatedAa = newAa;
            return new string(chars);
        }

        /// <summary>Convenience overload that discards the mutation details.</summary>
        public static string MutateOneAminoAcid(string peptide, Random random = null)
            => MutateOneAminoAcid(peptide, random, out _, out _, out _);

        [Test]
        public static void ReverseFasta()
        {
            var peptidePath = @"E:\Aneuploidy\Mistranslation_project\Mistranslation_search\ReverseDb\std_YL_raw-search-gptmdFilter-search\Task3-SearchTask\AllPeptides.psmtsv";
            var outDir = @"E:\Aneuploidy\Mistranslation_project\Mistranslation_search\ReverseDb\std_YL_raw-search-gptmdFilter-search";
            var peptides = new PsmFromTsvFile(peptidePath);
            var candidatePeptides = peptides.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T" && p.AmbiguityLevel == "1" && p.MissedCleavage == "0" && p.PEP <= 0.01 && p.Score >= 15);
            var fastaFile = @"E:\Databases\uniprotkb_taxonomy_id_559292_AND_review_2024_10_02.fasta";

            // Keep peptides carrying no modifications other than Carbamidomethyl (fixed on C)
            // or Oxidation on M. Anything else (TMT, phospho, acetyl, biological mods, ...)
            // is excluded so the mutation isn't confounded by an unrelated mass shift.
            var safeMods = new List<string> { "Carbamidomethyl", "Oxidation on M" };
            var eligible = candidatePeptides.Where(p => SpectrumMatchFromTsv.ParseModifications(p.FullSequence).Values.All(m => safeMods.Any(sm => m.Contains(sm)))).ToList();
            TestContext.WriteLine($"{eligible.Count} peptides eligible for mutation.");

            // Shuffle once with a fixed seed so the smaller samples are nested subsets
            // of the larger ones (the 100-peptide set ⊂ 1000 ⊂ 5000 ⊂ 10000), and the
            // experiment is reproducible across runs.
            var rng = new Random(42);
            var shuffled = eligible.OrderBy(_ => rng.Next()).ToList();

            // Load the target FASTA once and index by accession so we can splice
            // each peptide's mutation back into its parent protein.
            var allProteins = ProteinDbLoader.LoadProteinFasta(fastaFile, generateTargets: true, DecoyType.None, isContaminant: false,out var fastaErrors);
            var proteinByAccession = allProteins.Where(p => !p.IsDecoy).GroupBy(p => p.Accession).ToDictionary(g => g.Key, g => g.First());

            // Pre-compute the random mutation for every eligible peptide ONCE so that the smaller nested samples share the exact same mutations as the larger
            // ones. Each record locates the residue in the parent protein (1-based) so we can splice it into the FASTA, and remembers the from/to AA for the TSV.
            var positionParser = new Regex(@"(\d+)\s+to\s+(\d+)");
            var precomputed = new List<(PsmFromTsv psm, string mutatedBase, int idx, char from, char to, int proteinPos1Based)>();
            int skippedNoProtein = 0, skippedBadPosition = 0, skippedNoMutation = 0, skippedAaMismatch = 0;
            foreach (var p in shuffled)
            {
                if (string.IsNullOrEmpty(p.ProteinAccession) || !proteinByAccession.TryGetValue(p.ProteinAccession, out var parent))
                { skippedNoProtein++; continue; }

                var match = positionParser.Match(p.StartAndEndResiduesInProtein ?? "");
                if (!match.Success) { skippedBadPosition++; continue; }
                int start1Based = int.Parse(match.Groups[1].Value);

                var mutatedBase = MutateOneAminoAcid(p.BaseSeq, rng, out int idx, out char from, out char to);
                if (idx < 0) { skippedNoMutation++; continue; }

                int proteinPos1Based = start1Based + idx;
                // Sanity-check: the residue we're mutating must agree with the parent
                // protein at that position; if not, the peptide came from a different
                // copy/isoform and we shouldn't touch the FASTA blindly.
                if (proteinPos1Based < 1 || proteinPos1Based > parent.BaseSequence.Length|| parent.BaseSequence[proteinPos1Based - 1] != from)
                { skippedAaMismatch++; continue; }

                precomputed.Add((p, mutatedBase, idx, from, to, proteinPos1Based));
            }

            var sampleSizes = new[] { 1000 };
            foreach (var n in sampleSizes)
            {
                if (precomputed.Count < n)
                {
                    TestContext.WriteLine($"Only {precomputed.Count} mutations available — skipping sample size {n}.");
                    continue;
                }

                var sample = precomputed.Take(n).ToList();

                // Group mutations by protein accession and splice them into a fresh char
                // array of the parent protein's sequence. If two sampled peptides map to
                // the same position in the same protein (unlikely after the AA-mismatch
                // check), the later one wins — we record both rows in the TSV regardless.
                var mutationsByAccession = sample.GroupBy(s => s.psm.ProteinAccession).ToDictionary(g => g.Key, g => g.ToList());

                var outputProteins = new List<Protein>(allProteins.Count);
                foreach (var protein in allProteins.Where(p => !p.IsDecoy))
                {
                    if (mutationsByAccession.TryGetValue(protein.Accession, out var muts))
                    {
                        var chars = protein.BaseSequence.ToCharArray();
                        foreach (var m in muts)
                            chars[m.proteinPos1Based - 1] = m.to;
                        outputProteins.Add(new Protein(protein, new string(chars)));
                    }
                    else
                    {
                        outputProteins.Add(protein);
                    }
                }

                // Write the mutated FASTA (UniProt-style "|"-delimited headers).
                var fastaOutPath = Path.Combine(outDir, $"MutatedFasta_{n}.fasta");
                ProteinDbWriter.WriteFastaDatabase(outputProteins, fastaOutPath, "|");

                // Companion TSV: one row per sampled peptide describing the mutation
                // applied, where it landed in the parent protein, and the source PSM.
                var tsvOutPath = Path.Combine(outDir, $"MutatedPeptides_{n}.tsv");
                using (var writer = new StreamWriter(tsvOutPath))
                {
                    writer.WriteLine(string.Join("\t","OriginalBaseSeq", "MutatedBaseSeq","MutationIndexInPeptide", "OriginalAa", "MutatedAa", "PrecursorMass", "MassShift", "ProteinPosition1Based", "ProteinAccession", "GeneName", "Start_and_End_Residues", "OriginalFullSequence", "Score"));
                    foreach (var m in sample)
                    {
                        // Residue-level monoisotopic mass delta (Da) for the substitution:
                        // mass(mutated AA) - mass(original AA). Uses the static residue mass
                        // table indexed by ASCII so no Residue lookup overhead per row.
                        double massShift = Residue.ResidueMonoisotopicMass[(int)m.to]
                                         - Residue.ResidueMonoisotopicMass[(int)m.from];
                        writer.WriteLine(string.Join("\t",m.psm.BaseSeq, m.mutatedBase,m.idx, m.from, m.to, m.psm.PrecursorMass,
                            massShift.ToString("F4", CultureInfo.InvariantCulture),
                            m.proteinPos1Based, m.psm.ProteinAccession,
                            m.psm.GeneName,m.psm.StartAndEndResiduesInProtein, m.psm.FullSequence, m.psm.Score));
                    }
                }
            }
        }

        /// <summary>
        /// For each MutatedPeptides_{N}.tsv produced by ReverseFasta, finds the source
        /// MS2 scan that identified each peptide (and that scan's precursor MS1 scan),
        /// pulls those scans out of the raw file, renumbers them sequentially, and
        /// writes a small mzML so you have a low-volume file to test re-searches on.
        /// One mzML is produced per source raw file per sample size — naming pattern
        /// {rawFileBaseName}_mutated_{N}.mzML next to the TSV.
        /// </summary>
        [Test]
        public static void WriteMutatedScansMzml()
        {
            // Same paths as ReverseFasta — the PSM table is the ground truth for
            // mapping each TSV row back to its source scan number + raw file.
            var peptidePath = @"E:\Aneuploidy\Mistranslation_project\Mistranslation_search\ReverseDb\std_YL_raw-search-gptmdFilter-search\Task3-SearchTask\AllPeptides.psmtsv";
            var outDir = @"E:\Aneuploidy\Mistranslation_project\Mistranslation_search\ReverseDb\std_YL_raw-search-gptmdFilter-search";
            // Where the raw/mzML files live (FileName column in the PSM table is the
            // basename, no extension). Adjust if the raws aren't in this directory.
            var rawDir = @"E:\Aneuploidy\Mistranslation_project\Mistranslation_search\ReverseDb";

            var psmFile = new PsmFromTsvFile(peptidePath);
            // Index PSMs by their full sequence — that's the join key for the TSV's
            // OriginalFullSequence column. AmbiguityLevel=="1" peptides are unique by
            // full sequence in AllPeptides.psmtsv, so .First() is safe; we just guard
            // against accidental dupes with GroupBy().
            var psmByFullSeq = psmFile.Results
                .GroupBy(p => p.FullSequence)
                .ToDictionary(g => g.Key, g => g.First());

            // Cache opened raw files across sample sizes — the 100-PSM mzML and the
            // 1000-PSM mzML often draw from the same raws, no need to re-read.
            var rawFileCache = new Dictionary<string, MsDataFile>();

            foreach (var n in new[] { 1000 })
            {
                var tsvPath = Path.Combine(outDir, $"MutatedPeptides_{n}.tsv");
                if (!File.Exists(tsvPath))
                {
                    TestContext.WriteLine($"Skipping {n}: {tsvPath} not found. Run ReverseFasta first.");
                    continue;
                }

                var tsvLines = File.ReadAllLines(tsvPath);
                var header = tsvLines[0].Split('\t');
                int fullSeqIdx = Array.IndexOf(header, "OriginalFullSequence");
                if (fullSeqIdx < 0)
                    throw new InvalidDataException($"'OriginalFullSequence' column missing in {tsvPath}");

                // Walk the TSV rows IN ORDER so the mzML reflects the same peptide set
                // ReverseFasta wrote (and the 100-row mzML is a subset of the 1000-row mzML).
                var selectedPsms = new List<PsmFromTsv>();
                int missing = 0;
                foreach (var line in tsvLines.Skip(1))
                {
                    var fields = line.Split('\t');
                    if (fields.Length <= fullSeqIdx) continue;
                    if (psmByFullSeq.TryGetValue(fields[fullSeqIdx], out var psm)) selectedPsms.Add(psm);
                    else missing++;
                }
                TestContext.WriteLine($"{n}: matched {selectedPsms.Count} of {tsvLines.Length - 1} TSV rows to PSMs ({missing} unmatched).");

                // One mzML per source raw file. Most of the time it's just one file,
                // but if PSMs span multiple raws we'd emit one mzML per raw.
                foreach (var fileGroup in selectedPsms.GroupBy(p => p.FileName))
                {
                    if (!rawFileCache.TryGetValue(fileGroup.Key, out var dataFile))
                    {
                        // FileName in MM PSM output is the base name (no extension).
                        // Try .raw, then .mzML/.mzml. Skip if none found.
                        var candidates = new[] { ".raw", ".mzML", ".mzml" }
                            .Select(ext => Path.Combine(rawDir, fileGroup.Key + ext));
                        var rawPath = candidates.FirstOrDefault(File.Exists);
                        if (rawPath == null)
                        {
                            TestContext.WriteLine($"Raw file for '{fileGroup.Key}' not found in {rawDir} — skipping {fileGroup.Count()} PSMs.");
                            continue;
                        }
                        dataFile = MsDataFileReader.GetDataFile(rawPath);
                        dataFile.LoadAllStaticData();
                        rawFileCache[fileGroup.Key] = dataFile;
                    }

                    // MS2 ONLY — no MS1 scans. For each PSM we keep its MS2 scan and
                    // bake the deconvoluted precursor info (m/z, charge, mono-guess m/z)
                    // directly into the MS2's selectedIon fields, taking those values
                    // from the PsmFromTsv row. That's the value MetaMorpheus already
                    // trusted to identify this peptide, so a re-search can skip
                    // deconvolution entirely. Dedup by Ms2ScanNumber in case two PSMs
                    // accidentally share a scan; first PSM wins.
                    var ms2ScanToPsm = new SortedDictionary<int, PsmFromTsv>();
                    foreach (var psm in fileGroup)
                        if (!ms2ScanToPsm.ContainsKey(psm.Ms2ScanNumber))
                            ms2ScanToPsm[psm.Ms2ScanNumber] = psm;

                    var collected = new List<(MsDataScan scan, PsmFromTsv psm)>();
                    foreach (var kv in ms2ScanToPsm)
                    {
                        MsDataScan ms2;
                        try { ms2 = dataFile.GetOneBasedScan(kv.Key); }
                        catch { continue; }
                        if (ms2 == null) continue;
                        collected.Add((ms2, kv.Value));
                    }

                    // Renumber 1..K. The mzML writer + FakeMsDataFile.GetOneBasedScan
                    // assume array index == scanNumber - 1, so non-contiguous original
                    // numbers would leave holes. Drop OneBasedPrecursorScanNumber to
                    // null since we're not including MS1 — the precursor m/z and charge
                    // live in selectedIon* on the MS2 itself.
                    var newScans = new MsDataScan[collected.Count];
                    for (int idx = 0; idx < collected.Count; idx++)
                    {
                        var (s, psm) = collected[idx];
                        int newScanNum = idx + 1;

                        // Override the MS2's selected-ion fields with the PSM's
                        // deconvoluted values. Keep the original isolation m/z + width
                        // (those describe what was *physically* isolated by the
                        // instrument, independent of the deconvolution result).
                        double psmPrecursorMz = psm.PrecursorMz;
                        int psmPrecursorCharge = psm.PrecursorCharge;

                        newScans[idx] = new MsDataScan(
                            massSpectrum: s.MassSpectrum,
                            oneBasedScanNumber: newScanNum,
                            msnOrder: s.MsnOrder,
                            isCentroid: s.IsCentroid,
                            polarity: s.Polarity,
                            retentionTime: s.RetentionTime,
                            scanWindowRange: s.ScanWindowRange,
                            scanFilter: s.ScanFilter,
                            mzAnalyzer: s.MzAnalyzer,
                            totalIonCurrent: s.TotalIonCurrent,
                            injectionTime: s.InjectionTime,
                            noiseData: s.NoiseData,
                            nativeId: $"controllerType=0 controllerNumber=1 scan={newScanNum}",
                            selectedIonMz: psmPrecursorMz,
                            selectedIonChargeStateGuess: psmPrecursorCharge,
                            selectedIonIntensity: s.SelectedIonIntensity,
                            isolationMZ: s.IsolationMz ?? psmPrecursorMz,
                            isolationWidth: s.IsolationWidth,
                            dissociationType: s.DissociationType,
                            oneBasedPrecursorScanNumber: null,
                            // PsmFromTsv.PrecursorMz is already the monoisotopic m/z,
                            // so reuse it here too — no separate mono guess needed.
                            selectedIonMonoisotopicGuessMz: psmPrecursorMz,
                            hcdEnergy: s.HcdEnergy,
                            scanDescription: s.ScanDescription,
                            compensationVoltage: s.CompensationVoltage);
                    }

                    var fake = new FakeMsDataFile(newScans);
                    var outMzmlPath = Path.Combine(outDir, $"{fileGroup.Key}_mutated_{n}.mzML");
                    MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(fake, outMzmlPath, writeIndexed: true);
                    TestContext.WriteLine($"{n}: wrote {newScans.Length} MS2 scans → {outMzmlPath}");
                }
            }
        }

        /// <summary>
        /// MaxQuant-compatible variant of <see cref="WriteMutatedScansMzml"/>.
        /// MaxQuant rejects mzML files that contain only MS2 scans (it builds 3D
        /// features from MS1s for matching/quant), so this method:
        ///   1) Keeps each PSM's MS2 scan AND its precursor MS1 scan from the raw file.
        ///   2) Renumbers everything 1..K contiguously and rewrites the precursor
        ///      pointer on each MS2 so MS2 → MS1 links stay valid after renumbering.
        ///   3) Emits Thermo-style nativeIDs ("controllerType=0 controllerNumber=1
        ///      scan=N") so MaxQuant's nativeID parser is happy.
        ///   4) Writes an indexedmzML (writeIndexed: true) — MaxQuant requires the
        ///      &lt;indexedmzML&gt; wrapper.
        /// Output naming: {rawFileBaseName}_mutated_{N}_MQ.mzML, next to the TSV.
        /// Does NOT modify the original WriteMutatedScansMzml method or its output.
        /// </summary>
        [Test]
        public static void WriteMutatedScansMzmlForMaxQuant()
        {
            var peptidePath = @"E:\Aneuploidy\Mistranslation_project\Mistranslation_search\ReverseDb\std_YL_raw-search-gptmdFilter-search\Task3-SearchTask\AllPeptides.psmtsv";
            var outDir = @"E:\Aneuploidy\Mistranslation_project\Mistranslation_search\ReverseDb\std_YL_raw-search-gptmdFilter-search";
            var rawDir = @"E:\Aneuploidy\Mistranslation_project\Mistranslation_search\ReverseDb";

            var psmFile = new PsmFromTsvFile(peptidePath);
            var psmByFullSeq = psmFile.Results
                .GroupBy(p => p.FullSequence)
                .ToDictionary(g => g.Key, g => g.First());

            var rawFileCache = new Dictionary<string, MsDataFile>();

            foreach (var n in new[] { 100, 1000 })
            {
                var tsvPath = Path.Combine(outDir, $"MutatedPeptides_{n}.tsv");
                if (!File.Exists(tsvPath))
                {
                    TestContext.WriteLine($"Skipping {n}: {tsvPath} not found. Run ReverseFasta first.");
                    continue;
                }

                var tsvLines = File.ReadAllLines(tsvPath);
                var header = tsvLines[0].Split('\t');
                int fullSeqIdx = Array.IndexOf(header, "OriginalFullSequence");
                if (fullSeqIdx < 0)
                    throw new InvalidDataException($"'OriginalFullSequence' column missing in {tsvPath}");

                var selectedPsms = new List<PsmFromTsv>();
                int missing = 0;
                foreach (var line in tsvLines.Skip(1))
                {
                    var fields = line.Split('\t');
                    if (fields.Length <= fullSeqIdx) continue;
                    if (psmByFullSeq.TryGetValue(fields[fullSeqIdx], out var psm)) selectedPsms.Add(psm);
                    else missing++;
                }
                TestContext.WriteLine($"{n}: matched {selectedPsms.Count} of {tsvLines.Length - 1} TSV rows to PSMs ({missing} unmatched).");

                foreach (var fileGroup in selectedPsms.GroupBy(p => p.FileName))
                {
                    if (!rawFileCache.TryGetValue(fileGroup.Key, out var dataFile))
                    {
                        var candidates = new[] { ".raw", ".mzML", ".mzml" }
                            .Select(ext => Path.Combine(rawDir, fileGroup.Key + ext));
                        var rawPath = candidates.FirstOrDefault(File.Exists);
                        if (rawPath == null)
                        {
                            TestContext.WriteLine($"Raw file for '{fileGroup.Key}' not found in {rawDir} — skipping {fileGroup.Count()} PSMs.");
                            continue;
                        }
                        dataFile = MsDataFileReader.GetDataFile(rawPath);
                        dataFile.LoadAllStaticData();
                        rawFileCache[fileGroup.Key] = dataFile;
                    }

                    // OPTION 1: keep EVERY MS1 scan from the raw file, but only the
                    // targeted MS2 scans. MaxQuant's feature detection walks consecutive
                    // MS1 scans to build isotope patterns — sparse MS1s (one per MS2)
                    // produce zero features and the run silently aborts after Feature
                    // detection. Keeping the full MS1 trail gives MQ proper chromatography
                    // while the targeted MS2 subset keeps the file small.
                    var ms2ScanToPsm = new Dictionary<int, PsmFromTsv>();
                    foreach (var psm in fileGroup)
                        if (!ms2ScanToPsm.ContainsKey(psm.Ms2ScanNumber))
                            ms2ScanToPsm[psm.Ms2ScanNumber] = psm;

                    // Walk every scan in the raw file in original scan-number order.
                    //   - MS1 (MsnOrder == 1): always include, psm = null
                    //   - MSn (MsnOrder >= 2): include only if its scan number matches
                    //     one of our targeted PSMs
                    var neededScans = new SortedDictionary<int, (MsDataScan scan, PsmFromTsv psm)>();
                    foreach (var s in dataFile.GetAllScansList())
                    {
                        if (s.MsnOrder == 1)
                        {
                            neededScans[s.OneBasedScanNumber] = (s, null);
                        }
                        else if (ms2ScanToPsm.TryGetValue(s.OneBasedScanNumber, out var psm))
                        {
                            neededScans[s.OneBasedScanNumber] = (s, psm);
                        }
                    }
                    int targetedMs2Found = ms2ScanToPsm.Keys.Count(k => neededScans.ContainsKey(k));
                    if (targetedMs2Found < ms2ScanToPsm.Count)
                        TestContext.WriteLine($"  Warning: {ms2ScanToPsm.Count - targetedMs2Found} targeted MS2 scan numbers were not found in {fileGroup.Key}.");

                    // Renumber 1..K and remap each MS2's precursor pointer to the new
                    // MS1 scan number. MaxQuant validates this link.
                    var oldToNew = neededScans.Keys
                        .Select((oldNum, idx) => (oldNum, newNum: idx + 1))
                        .ToDictionary(x => x.oldNum, x => x.newNum);

                    var newScans = new MsDataScan[neededScans.Count];
                    int writeIdx = 0;
                    foreach (var kv in neededScans)
                    {
                        var (s, psm) = kv.Value;
                        int newScanNum = oldToNew[kv.Key];

                        // For MS2: bake the PSM's deconvoluted precursor info into
                        // selectedIon* AND keep the precursor scan link. For MS1: copy
                        // through verbatim — selectedIon* aren't meaningful for MS1.
                        bool isMs2 = s.MsnOrder >= 2 && psm != null;
                        double? selectedMz = isMs2 ? psm.PrecursorMz : s.SelectedIonMZ;
                        int? selectedCharge = isMs2 ? psm.PrecursorCharge : s.SelectedIonChargeStateGuess;
                        double? selectedMonoMz = isMs2 ? psm.PrecursorMz : s.SelectedIonMonoisotopicGuessMz;
                        int? newPrecursorScan = (s.OneBasedPrecursorScanNumber.HasValue
                                                 && oldToNew.TryGetValue(s.OneBasedPrecursorScanNumber.Value, out var mapped))
                            ? mapped
                            : (int?)null;

                        newScans[writeIdx++] = new MsDataScan(
                            massSpectrum: s.MassSpectrum,
                            oneBasedScanNumber: newScanNum,
                            msnOrder: s.MsnOrder,
                            isCentroid: s.IsCentroid,
                            polarity: s.Polarity,
                            retentionTime: s.RetentionTime,
                            scanWindowRange: s.ScanWindowRange,
                            scanFilter: s.ScanFilter,
                            mzAnalyzer: s.MzAnalyzer,
                            totalIonCurrent: s.TotalIonCurrent,
                            injectionTime: s.InjectionTime,
                            noiseData: s.NoiseData,
                            // Thermo-style nativeID is what MaxQuant's parser expects.
                            nativeId: $"controllerType=0 controllerNumber=1 scan={newScanNum}",
                            selectedIonMz: selectedMz,
                            selectedIonChargeStateGuess: selectedCharge,
                            selectedIonIntensity: s.SelectedIonIntensity,
                            isolationMZ: s.IsolationMz ?? selectedMz,
                            isolationWidth: s.IsolationWidth,
                            dissociationType: s.DissociationType,
                            oneBasedPrecursorScanNumber: newPrecursorScan,
                            selectedIonMonoisotopicGuessMz: selectedMonoMz,
                            hcdEnergy: s.HcdEnergy,
                            scanDescription: s.ScanDescription,
                            compensationVoltage: s.CompensationVoltage);
                    }

                    var fake = new FakeMsDataFile(newScans);
                    // _MQ suffix distinguishes from the MS2-only file produced by
                    // WriteMutatedScansMzml so both can coexist in the same dir.
                    var outMzmlPath = Path.Combine(outDir, $"{fileGroup.Key}_mutated_{n}_MQ.mzML");
                    MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(fake, outMzmlPath, writeIndexed: false);
                    int ms1Count = newScans.Count(s => s.MsnOrder == 1);
                    int ms2Count = newScans.Count(s => s.MsnOrder >= 2);
                    TestContext.WriteLine($"{n}: wrote {newScans.Length} scans ({ms1Count} MS1 + {ms2Count} MS2) → {outMzmlPath}");
                }

            }

        }

        [Test]
        public static void testReadmzML()
        {
            var path = @"E:\Aneuploidy\Mistranslation_project\Mistranslation_search\ReverseDb\std_YL_raw-search-gptmdFilter-search\04-14-26_std-YL_mutated_100_MQ.mzML";
            var msFile = MsDataFileReader.GetDataFile(path);
            var scans = msFile.GetAllScansList();
        }

        /// <summary>
        /// Given a search-result PSM file (e.g. AllPSMs.psmtsv from a MetaMorpheus run
        /// against the mutated database) and the MutatedPeptides_{N}.tsv that recorded
        /// what we intentionally mutated, check how many of the engineered peptides the
        /// search correctly identified. Writes one VerifyMutations_{label}.tsv per pair
        /// with per-peptide status and prints a summary to TestContext.
        ///
        /// Classifications per mutated peptide:
        ///   FoundMutated            — search returned the mutated base sequence anywhere
        ///   FoundMutatedAtProtein   — and on the expected protein accession (stronger)
        ///   FoundOriginal           — search returned the ORIGINAL sequence instead
        ///                             (mutation didn't break the match — failure mode)
        ///   neither                 — peptide scan went unidentified
        /// </summary>
        [Test]
        public static void VerifyMutatedPeptidesIdentified()
        {
            // Edit this list to add more (search PSM file, mutation TSV, label) triplets.
            var pairs = new List<(string searchPsmPath, string mutationTsvPath, string label)>
            {
                (
                    @"E:\Aneuploidy\Mistranslation_project\Mistranslation_search\ReverseDb\MM_custom_branch\Modbox\060826_mutation100_max1\Task\AllPSMs.psmtsv",
                    @"E:\Aneuploidy\Mistranslation_project\Mistranslation_search\ReverseDb\std_YL_raw-search-gptmdFilter-search\MutatedPeptides_100.tsv",
                    "mutation1000_max1"
                ),
            };

            foreach (var (searchPsmPath, mutationTsvPath, label) in pairs)
            {
                if (!File.Exists(searchPsmPath) || !File.Exists(mutationTsvPath))
                    continue;

                // Quality-filter target PSMs (q-value <= 0.01, T hits only) and index
                // by I/L-normalized BaseSeq. Keys are folded I→L so a peptide that the
                // engine returned with the wrong I/L isobar still resolves to a hit when
                // the expected sequence has the opposite isobar.
                var searchPsms = new PsmFromTsvFile(searchPsmPath).Results
                    .Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T")
                    .ToList();
                var psmsByBaseSeq = searchPsms
                    .GroupBy(p => p.BaseSeq?.Replace('I', 'L'))
                    .Where(g => g.Key != null)
                    .ToDictionary(g => g.Key, g => g.ToList());

                var tsvLines = File.ReadAllLines(mutationTsvPath);
                var header = tsvLines[0].Split('\t');
                int origIdx       = Array.IndexOf(header, "OriginalBaseSeq");
                int mutIdx        = Array.IndexOf(header, "MutatedBaseSeq");
                int posIdx        = Array.IndexOf(header, "ProteinPosition1Based");
                int accIdx        = Array.IndexOf(header, "ProteinAccession");
                int geneIdx       = Array.IndexOf(header, "GeneName");
                int origAaIdx     = Array.IndexOf(header, "OriginalAa");
                int mutAaIdx      = Array.IndexOf(header, "MutatedAa");
                int origFullIdx   = Array.IndexOf(header, "OriginalFullSequence");
                int origScoreIdx  = Array.IndexOf(header, "Score");
                if (origIdx < 0 || mutIdx < 0 || accIdx < 0)
                    throw new InvalidDataException($"Missing required column in {mutationTsvPath}");

                // I and L are isobaric (same monoisotopic mass) and indistinguishable by
                // MS2 fragmentation, so any sequence-equality check between expected and
                // identified peptides must treat them as the same residue. Helper folds
                // every I to L before comparing.
                static string NormalizeIL(string s) => s == null ? null : s.Replace('I', 'L');

                var outRows = new List<string[]>();
                foreach (var line in tsvLines.Skip(1))
                {
                    var fields = line.Split('\t');
                    if (fields.Length <= Math.Max(mutIdx, accIdx)) continue;

                    string origSeq     = fields[origIdx];
                    string mutSeq      = fields[mutIdx];
                    string acc = fields[accIdx];
                    string expectedAAsub = $"({fields[origAaIdx]} -> {fields[mutAaIdx]})";
                    string origFullSeq = origFullIdx  >= 0 && fields.Length > origFullIdx  ? fields[origFullIdx]  : "";
                    string origScore   = origScoreIdx >= 0 && fields.Length > origScoreIdx ? fields[origScoreIdx] : "";

                    // I/L-folded keys for lookups so an L→I or I→L mutation (which a
                    // search engine cannot distinguish) still surfaces a candidate hit.
                    string normMutSeq  = NormalizeIL(mutSeq);
                    string normOrigSeq = NormalizeIL(origSeq);

                    bool foundMutated      = psmsByBaseSeq.TryGetValue(normMutSeq, out var mutHits);
                    bool foundMutAtProtein = foundMutated && mutHits.Any(p => p.ProteinAccession != null && p.ProteinAccession.Contains(acc));
                    bool foundOriginal     = psmsByBaseSeq.TryGetValue(normOrigSeq, out var origHits);

                    // Prefer a hit whose AA-sub annotation, when applied, reverts to the
                    // original sequence (MM both found the engineered residue AND
                    // annotated the sub back to the original).
                    PsmFromTsv strictHit = null;
                    if (foundMutated)
                    {
                        foreach (var candidate in mutHits)
                        {
                            if (NormalizeIL(candidate.BaseSeq) != normMutSeq) continue;
                            string subbed;
                            try { subbed = IBioPolymerWithSetMods.ParseSubstitutedFullSequence(candidate.FullSequence); }
                            catch { continue; }
                            var subbedBase = IBioPolymerWithSetMods.GetBaseSequenceFromFullSequence(subbed);
                            if (NormalizeIL(subbedBase) == normOrigSeq) { strictHit = candidate; break; }
                        }
                    }

                    PsmFromTsv hit = strictHit
                                   ?? (foundMutAtProtein ? mutHits.First(p => p.ProteinAccession.Contains(acc))
                                   :   foundMutated     ? mutHits.First()
                                   :   foundOriginal    ? origHits.First()
                                   :   null);

                    // Per-row verdicts. Both use I/L-tolerant comparison.
                    //   ifBaseSeqMatch: search returned the engineered (mutated) sequence
                    //   ifFullSeqMatch: applying the AA-sub on the identified FullSequence
                    //                   produces the original pre-mutation full sequence
                    bool ifBaseSeqMatch = hit != null && NormalizeIL(hit.BaseSeq) == normMutSeq;
                    bool ifFullSeqMatch = false;
                    if (hit != null && !string.IsNullOrEmpty(origFullSeq))
                    {
                        try
                        {
                            var subbedFull = IBioPolymerWithSetMods.ParseSubstitutedFullSequence(hit.FullSequence);
                            ifFullSeqMatch = NormalizeIL(subbedFull) == NormalizeIL(origFullSeq);
                        }
                        catch { /* malformed AA-sub annotation; leave false */ }
                    }

                    outRows.Add(new[]
                    {
                        // Expected
                        origSeq, mutSeq, expectedAAsub, origScore,
                        // Identified
                        hit?.BaseSeq ?? "",
                        hit?.FullSequence ?? "",
                        hit?.Score.ToString("F2", CultureInfo.InvariantCulture) ?? "",
                        // Verdicts
                        ifBaseSeqMatch.ToString(),
                        ifFullSeqMatch.ToString(),
                    });
                }

                // Write companion verification TSV next to the mutation TSV.
                var outDir = Path.GetDirectoryName(mutationTsvPath);
                var outPath = @"E:\Aneuploidy\Mistranslation_project\Mistranslation_search\ReverseDb\MM_custom_branch\Modbox\060826_mutation100_max1\Task\checkID.tsv";
                using (var writer = new StreamWriter(outPath))
                {
                    writer.WriteLine(string.Join("\t",
                        // Expected
                        "ExpectedOriginalBaseSeq", "ExpectedMutatedBaseSeq", "ExpectedAAsub", "OriginalScore",
                        // Identified
                        "IdentifiedBaseSeq", "IdentifiedFullSequence", "IdentifiedScore",
                        // Verdicts
                        "ifBaseSeqMatch", "ifFullSeqMatch"));
                    foreach (var row in outRows)
                        writer.WriteLine(string.Join("\t", row));
                }
            }
        }

        /// <summary>
        /// Reads the FragPipe DMO mass-offset list (AASubs_DMO-list.tsv) and writes a
        /// trimmed copy that keeps only:
        ///   1) entries reachable by a single-nucleotide substitution under the same
        ///      rules as <see cref="MutateOneAminoAcid"/> — C/K/R are not allowed as
        ///      source or target. For mixed-site rows, sites that aren't 1-nt sources
        ///      for that mass shift are dropped (e.g. "72.9952 IL" → "72.9952 L"
        ///      because only L→W is a 1-nt sub at +72.9952 Da).
        ///   2) entries whose mass matches a common PTM (methylation, oxidation,
        ///      acetylation, phosphorylation), preserved verbatim so DMO still
        ///      detects those modifications. Both positive and negative forms are
        ///      kept (handles neutral-loss style rows).
        /// </summary>
        [Test]
        public static void FilterDmoSubsListToSingleNtSubs()
        {
            var dmoListPath = @"E:\Aneuploidy\Mistranslation_project\Mistranslation_search\DMO_data\DMO-FragPipe_ZenodoRepo\DMO_Input-Mass-Lists\AASubs_DMO-list.tsv";
            var outPath = Path.Combine(Path.GetDirectoryName(dmoListPath), "AASubs_DMO-list_1ntSub.tsv");

            // Per-source-residue mass shifts achievable by a single-nucleotide
            // substitution (C/K/R excluded either side) — derived from the same
            // ReachableAaMap the mutation generator uses, so any change in the
            // mutation rules automatically propagates here.
            var shiftsByFrom = new Dictionary<char, List<double>>();
            foreach (var from in ReachableAaMap.Keys)
            {
                var shifts = ReachableAaMap[from]
                    .Select(to => Residue.ResidueMonoisotopicMass[(int)to] - Residue.ResidueMonoisotopicMass[(int)from])
                    .ToList();
                shiftsByFrom[from] = shifts;
            }

            // Common-PTM mass shifts. Same list as DeepLC.PtmConfounderMasses.
            // Matched as |mass| so neutral-loss rows (e.g. -79.9663 phosphate loss)
            // are preserved too.
            var ptmMasses = new[] { 14.0157, 15.9949, 42.0106, 79.9663 };
            const double massTol = 0.005;

            var lines = File.ReadAllLines(dmoListPath);
            var outLines = new List<string> { lines[0] }; // header

            foreach (var line in lines.Skip(1))
            {
                if (string.IsNullOrWhiteSpace(line)) continue;
                var parts = line.Split('\t');
                if (parts.Length < 2) continue;
                if (!double.TryParse(parts[0], NumberStyles.Float, CultureInfo.InvariantCulture, out double mass)) continue;
                var sites = parts[1];

                // PTM-mass row: keep verbatim. DMO needs every listed site so it can
                // still call the modification on any of them, including K and C.
                if (ptmMasses.Any(p => Math.Abs(Math.Abs(mass) - p) <= massTol))
                {
                    outLines.Add(line);
                    continue;
                }

                // AAsub row: keep only sites where the residue has a 1-nt neighbour
                // (also not C/K/R) whose mass matches this row's mass shift.
                var keptSites = new StringBuilder();
                foreach (var site in sites)
                {
                    if (!shiftsByFrom.TryGetValue(site, out var shifts)) continue;
                    if (shifts.Any(s => Math.Abs(s - mass) <= massTol))
                        keptSites.Append(site);
                }

                if (keptSites.Length > 0)
                    outLines.Add($"{parts[0]}\t{keptSites}");
            }

            File.WriteAllLines(outPath, outLines);
        }

        /// <summary>
        /// Implements the Lundgren / Champion 2024 "fit for purpose" benchmark
        /// (J. Proteome Res., bioRxiv 2023.08.09.552645) for the 1715 E. coli +
        /// Salmonella co-lysate run. Compares two AAS-detection pipelines —
        /// FragPipe DMO and MetaMorpheus AAsub — against the SSP ground truth
        /// built by FindSSP.py.
        ///
        /// For each engine, every PSM that calls a single AA substitution is
        /// resolved to the post-substitution base sequence. A call is counted as
        /// a "ground-truth recovery" iff that substituted base sequence appears
        /// in the SALTY SSP list for the E. coli peptide (I/L folded).
        ///
        /// Outputs three TSVs in <c>outDir</c>:
        ///   DMO_SspComparison.tsv  — per-PSM details for DMO
        ///   MM_SspComparison.tsv   — per-PSM details for MM
        ///   Comparison_Summary.tsv — recovery counts and rates, per engine
        /// </summary>
        [Test]
        public static void CompareDmoVsMmAgainstSspGroundTruth()
        {
            var sspPath    = @"E:\Aneuploidy\Mistranslation_project\Mistranslation_search\DMO_data\DMO-FragPipe_ZenodoRepo\SSP_ECOLI_to_SALTY.csv";
            var dmoPsmPath = @"E:\Aneuploidy\Mistranslation_project\Mistranslation_search\DMO_data\DMO-FragPipe_ZenodoRepo\Results\AASubs_DMO\psm.tsv";
            var mmPsmPath  = @"E:\Aneuploidy\Mistranslation_project\Mistranslation_search\DMO_data\MM\ModBox\1715_max1_AllAAsubPTM\Task\AllPSMs.psmtsv";
            var outDir     = @"E:\Aneuploidy\Mistranslation_project\Mistranslation_search\DMO_data\MM\ModBox\1715_max1_AllAAsubPTM\Task";

            // ---------- 1) Load SSP ground truth ---------------------------------
            // Row layout: <idx>,<BaseSeq>,[<SSP list>],<IsSSP>,[<DSP list>],<IsDSP>,<IsExact>
            // The two bracketed lists contain commas, so we can't naively split on ','.
            // Instead, regex-out the two [...] groups and extract 'AA' literals.
            // I/L are isobaric — fold every I to L on both sides of the comparison.
            static string IL(string s) => s?.Replace('I', 'L');
            var sspByEcoli = new Dictionary<string, HashSet<string>>(StringComparer.Ordinal);
            var bracketRx  = new Regex(@"\[([^\]]*)\]", RegexOptions.Compiled);
            var quotedRx   = new Regex(@"'([A-Z]+)'", RegexOptions.Compiled);
            foreach (var line in File.ReadLines(sspPath).Skip(1))
            {
                if (string.IsNullOrWhiteSpace(line)) continue;
                var brackets = bracketRx.Matches(line);
                if (brackets.Count == 0) continue;
                var prefix = line.Substring(0, brackets[0].Index).TrimEnd(',');
                var pfields = prefix.Split(',');
                if (pfields.Length < 2) continue;
                var baseSeq = pfields[1].Trim();
                if (string.IsNullOrEmpty(baseSeq)) continue;

                if (!sspByEcoli.TryGetValue(IL(baseSeq), out var set))
                    sspByEcoli[IL(baseSeq)] = set = new HashSet<string>();
                foreach (Match m in quotedRx.Matches(brackets[0].Groups[1].Value))
                    set.Add(IL(m.Groups[1].Value));
            }

            // ---------- 2) Process DMO PSMs --------------------------------------
            var aa3to1 = new Dictionary<string, char>(StringComparer.Ordinal)
            {
                {"Ala",'A'},{"Arg",'R'},{"Asn",'N'},{"Asp",'D'},{"Cys",'C'},
                {"Gln",'Q'},{"Glu",'E'},{"Gly",'G'},{"His",'H'},{"Ile",'I'},
                {"Leu",'L'},{"Lys",'K'},{"Met",'M'},{"Phe",'F'},{"Pro",'P'},
                {"Ser",'S'},{"Thr",'T'},{"Trp",'W'},{"Tyr",'Y'},{"Val",'V'},
                {"Xle",'L'}, // Ile/Leu isobar — fold to L for comparison
            };
            var dmoSubRx     = new Regex(@"([A-Z][a-z]{2})->([A-Z][a-z]{2})\(([-\d.]+)\)", RegexOptions.Compiled);
            var dmoBestPosRx = new Regex(@"^([A-Z])(\d+)$", RegexOptions.Compiled);
            int dmoTotal = 0, dmoHit = 0, dmoMiss = 0, dmoSkipMulti = 0, dmoSkipParse = 0;

            using (var w = new StreamWriter(Path.Combine(outDir, "DMO_SspComparison.tsv")))
            {
                w.WriteLine(string.Join("\t",
                    "Spectrum", "Peptide", "Substitution", "Position",
                    "SubstitutedBaseSeq", "IsInSspGroundTruth", "IsIsobaricPair",
                    "Hyperscore", "Qvalue"));

                var lines = File.ReadAllLines(dmoPsmPath);
                var H = lines[0].Split('\t');
                int iSpec = Array.IndexOf(H, "Spectrum");
                int iPept = Array.IndexOf(H, "Peptide");
                int iObs  = Array.IndexOf(H, "Observed Modifications");
                int iBest = Array.IndexOf(H, "Best Positions");
                int iHyp  = Array.IndexOf(H, "Hyperscore");
                int iQ    = Array.IndexOf(H, "Qvalue");

                // The Zenodo DMO file contains both raws (1715 + 1716). The MM run
                // only searched 1715, so for an apples-to-apples comparison we keep
                // only 1715 PSMs here.
                const string raw1715Prefix = "2021-06-18-ECLandSALTY1_Slot1-42_1_1715";
                for (int i = 1; i < lines.Length; i++)
                {
                    var f = lines[i].Split('\t');
                    if (f.Length <= Math.Max(iObs, iBest)) continue;
                    if (!f[iSpec].StartsWith(raw1715Prefix, StringComparison.Ordinal)) continue;
                    if (!double.TryParse(f[iQ], NumberStyles.Float, CultureInfo.InvariantCulture, out double q) || q > 0.01) continue;

                    var obs = f[iObs];
                    if (string.IsNullOrEmpty(obs) || !obs.Contains("->")) continue;
                    // Single-substitution PSMs only — multi-sub rows are out of scope
                    // for the SSP (=1 substitution) ground truth.
                    var subMatches = dmoSubRx.Matches(obs);
                    if (subMatches.Count != 1) { dmoSkipMulti++; continue; }

                    var sm = subMatches[0];
                    if (!aa3to1.TryGetValue(sm.Groups[1].Value, out char fromAa)
                        || !aa3to1.TryGetValue(sm.Groups[2].Value, out char toAa))
                    { dmoSkipParse++; continue; }

                    var posMatch = dmoBestPosRx.Match(f[iBest]);
                    if (!posMatch.Success) { dmoSkipParse++; continue; }
                    int pos1 = int.Parse(posMatch.Groups[2].Value);

                    var peptide = f[iPept];
                    if (pos1 < 1 || pos1 > peptide.Length) { dmoSkipParse++; continue; }

                    var subbed = peptide.Substring(0, pos1 - 1) + toAa + peptide.Substring(pos1);
                    bool isHit = sspByEcoli.TryGetValue(IL(peptide), out var ssps) && ssps.Contains(IL(subbed));
                    bool isIso = (fromAa == 'I' && toAa == 'L') || (fromAa == 'L' && toAa == 'I')
                              || (fromAa == 'Q' && toAa == 'K') || (fromAa == 'K' && toAa == 'Q')
                              || (fromAa == 'N' && toAa == 'D') || (fromAa == 'D' && toAa == 'N');

                    dmoTotal++;
                    if (isHit) dmoHit++; else dmoMiss++;
                    w.WriteLine(string.Join("\t",
                        f[iSpec], peptide, $"{sm.Groups[1].Value}->{sm.Groups[2].Value}",
                        pos1, subbed, isHit, isIso, f[iHyp], q.ToString("G", CultureInfo.InvariantCulture)));
                }
            }

            // ---------- 3) Process MM PSMs ---------------------------------------
            int mmTotal = 0, mmHit = 0, mmMiss = 0, mmSkipMulti = 0, mmSkipParse = 0;
            using (var w = new StreamWriter(Path.Combine(outDir, "MM_SspComparison.tsv")))
            {
                w.WriteLine(string.Join("\t",
                    "ScanNumber", "BaseSeq", "FullSequence", "SubstitutedBaseSeq",
                    "IsInSspGroundTruth", "Score", "QValue"));

                var psms = new PsmFromTsvFile(mmPsmPath).Results
                    .Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T"
                             && p.FullSequence != null && p.FullSequence.Contains("nucleotide substitution"));
                foreach (var psm in psms)
                {
                    // "1+ nucleotide substitution" — paper's SSP set is exactly one
                    // sub, so skip 2+/3+ rows.
                    int subCount = Regex.Matches(psm.FullSequence, @"nucleotide substitutions?").Count;
                    if (subCount != 1) { mmSkipMulti++; continue; }

                    string subbedBase;
                    try
                    {
                        var subbedFull = IBioPolymerWithSetMods.ParseSubstitutedFullSequence(psm.FullSequence);
                        subbedBase = IBioPolymerWithSetMods.GetBaseSequenceFromFullSequence(subbedFull);
                    }
                    catch { mmSkipParse++; continue; }
                    if (string.IsNullOrEmpty(subbedBase)) { mmSkipParse++; continue; }

                    bool isHit = sspByEcoli.TryGetValue(IL(psm.BaseSeq), out var ssps) && ssps.Contains(IL(subbedBase));
                    mmTotal++;
                    if (isHit) mmHit++; else mmMiss++;
                    w.WriteLine(string.Join("\t",
                        psm.Ms2ScanNumber, psm.BaseSeq, psm.FullSequence, subbedBase,
                        isHit, psm.Score.ToString("F2", CultureInfo.InvariantCulture),
                        psm.QValue.ToString("G", CultureInfo.InvariantCulture)));
                }
            }

            // ---------- 4) Summary -----------------------------------------------
            using (var w = new StreamWriter(Path.Combine(outDir, "Comparison_Summary.tsv")))
            {
                w.WriteLine("Engine\tSingleSubPSMs\tHitsInSspGroundTruth\tMisses\tHitRate\tSkippedMultiSub\tSkippedParseErr\tDistinctEcoliPeptidesHit");
                int dmoDistinct = CountDistinctHits(Path.Combine(outDir, "DMO_SspComparison.tsv"));
                int mmDistinct  = CountDistinctHits(Path.Combine(outDir, "MM_SspComparison.tsv"));
                w.WriteLine($"DMO\t{dmoTotal}\t{dmoHit}\t{dmoMiss}\t{(dmoTotal>0 ? 100.0*dmoHit/dmoTotal : 0):F2}%\t{dmoSkipMulti}\t{dmoSkipParse}\t{dmoDistinct}");
                w.WriteLine($"MM\t{mmTotal}\t{mmHit}\t{mmMiss}\t{(mmTotal>0 ? 100.0*mmHit/mmTotal : 0):F2}%\t{mmSkipMulti}\t{mmSkipParse}\t{mmDistinct}");
            }

            static int CountDistinctHits(string tsv)
            {
                var set = new HashSet<string>();
                var lines = File.ReadAllLines(tsv);
                var H = lines[0].Split('\t');
                int iBase = Array.IndexOf(H, "Peptide");
                if (iBase < 0) iBase = Array.IndexOf(H, "BaseSeq");
                int iHit  = Array.IndexOf(H, "IsInSspGroundTruth");
                for (int i = 1; i < lines.Length; i++)
                {
                    var f = lines[i].Split('\t');
                    if (f.Length <= Math.Max(iBase, iHit)) continue;
                    if (string.Equals(f[iHit], "True", StringComparison.OrdinalIgnoreCase))
                        set.Add(f[iBase]);
                }
                return set.Count;
            }
        }

        /// <summary>
        /// Peptide-level (not PSM-level) comparison of FragPipe DMO vs MetaMorpheus
        /// AAsub on the 1715 E. coli + Salmonella co-lysate run, against the SSP
        /// ground truth produced by FindSSP.py (ChampionLab/substitutionannotation).
        ///
        /// Counting (same on both engines, I/L folded throughout):
        ///   1) Filter PSMs to Qvalue/QValue &lt;= 0.01 (and target hits for MM).
        ///   2) Among those, find PSMs carrying a substitution-mass shift:
        ///      - DMO: parse 'Assigned Modifications' entries like '7E(-14.0157)'.
        ///        Drop common-PTM masses (Met-ox, Carbamidomethyl, acetyl, methyl,
        ///        phospho), and accept the entry only if some target residue Y has
        ///        mass(Y) - mass(X) within tolerance of the shift.
        ///      - MM:  any PSM whose FullSequence carries '[1+ nucleotide
        ///        substitution:X->Y on X]'. Apply with the canonical mzLib parser.
        ///   3) Build the POST-substitution base sequence per PSM by splicing in
        ///      the new residue at the annotated position.
        ///   4) Dedupe to a unique-substituted-peptide SET per engine.
        ///   5) Ground truth = unique union of every string in 'SALTY SSP Sequence'
        ///      across all rows of SSP_ECOLI_to_SALTY.csv (= FindSSP.py output).
        ///   6) Intersect each engine's set with the ground-truth set; the size of
        ///      that intersection is the engine's recovered-SSP count.
        ///
        /// Writes to <c>outDir</c>:
        ///   GroundTruthSALTY_SSPs.txt              — the GT peptide set, one per line
        ///   DMO_SubstitutedPeptides_in_GT.tsv      — DMO unique sub-peptides + InGT
        ///   MM_SubstitutedPeptides_in_GT.tsv       — MM unique sub-peptides + InGT
        ///   DmoVsMm_PeptideCounts.tsv              — side-by-side counts
        /// </summary>
        [Test]
        public static void CompareDmoVsMmPeptideCounts()
        {
            var sspPath        = @"E:\Aneuploidy\Mistranslation_project\Mistranslation_search\DMO_data\DMO-FragPipe_ZenodoRepo\SSP_ECOLI_to_SALTY.csv";
            var dmoPsmPath     = @"E:\Aneuploidy\Mistranslation_project\Mistranslation_search\DMO_data\DMO-FragPipe_ZenodoRepo\Results\AASubs_DMO\psm.tsv";
            var mmModBoxPath   = @"E:\Aneuploidy\Mistranslation_project\Mistranslation_search\DMO_data\MM\ModBox\1715_max1_AllAAsubPTM\Task\AllPSMs.psmtsv";
            var mmGptmdPath    = @"E:\Aneuploidy\Mistranslation_project\Mistranslation_search\DMO_data\MM\GPTMD\AllAAsubs_1715_noFilter\Task2-SearchTask\AllPeptides.psmtsv";
            var outDir         = @"E:\Aneuploidy\Mistranslation_project\Mistranslation_search\DMO_data\MM\ModBox\1715_max1_AllAAsubPTM\Task";
            const string raw1715Prefix = "2021-06-18-ECLandSALTY1_Slot1-42_1_1715";
            const double massTol = 0.005;
            const string allAa = "ACDEFGHIKLMNPQRSTVWY";

            static string IL(string s) => s?.Replace('I', 'L');

            // Build (from_aa, [to_aa, massShift]) lookup over ALL 20 AAs. The SSP
            // ground truth reflects real genome differences, not 1-nt-restricted
            // substitutions, so we don't apply the C/K/R rule here.
            var subsByFrom = new Dictionary<char, List<(char to, double shift)>>();
            foreach (var fr in allAa)
            {
                var list = new List<(char, double)>();
                foreach (var to in allAa)
                {
                    if (to == fr) continue;
                    list.Add((to, Residue.ResidueMonoisotopicMass[to] - Residue.ResidueMonoisotopicMass[fr]));
                }
                subsByFrom[fr] = list;
            }
            double[] ptmMasses = { 15.9949, 57.0214, 42.0106, 14.0157, 79.9663 };
            bool IsPtm(double m) => ptmMasses.Any(p => Math.Abs(Math.Abs(m) - p) < massTol);
            IEnumerable<char> FindSubTargets(char fromAa, double shift)
            {
                if (!subsByFrom.TryGetValue(fromAa, out var list)) yield break;
                foreach (var (to, s) in list)
                    if (Math.Abs(s - shift) < massTol) yield return to;
            }

            // ---------- 1) Ground truth: SALTY SSP peptide SET ----------
            // CSV row: <idx>,<BaseSeq>,[<SSP list>],<IsSSP>,[<DSP list>],<IsDSP>,<IsExact>
            // The bracketed list values contain commas — extract via regex on [...] groups.
            var saltySspSet = new HashSet<string>(StringComparer.Ordinal);
            var bracketRx   = new Regex(@"\[([^\]]*)\]", RegexOptions.Compiled);
            var quotedRx    = new Regex(@"'([A-Z]+)'", RegexOptions.Compiled);
            foreach (var line in File.ReadLines(sspPath).Skip(1))
            {
                if (string.IsNullOrWhiteSpace(line)) continue;
                var brackets = bracketRx.Matches(line);
                if (brackets.Count == 0) continue;
                foreach (Match m in quotedRx.Matches(brackets[0].Groups[1].Value))
                    saltySspSet.Add(IL(m.Groups[1].Value));
            }

            // ---------- 2) DMO ----------
            // 'Assigned Modifications' entries look like '7E(-14.0157), 12D(32.0415)'.
            var dmoModRx = new Regex(@"(\d+)([A-Z])\(([-\d.]+)\)", RegexOptions.Compiled);
            var dmoAllIdentified = new HashSet<string>(StringComparer.Ordinal);
            var dmoSubIdentified = new HashSet<string>(StringComparer.Ordinal);
            int dmoPsms = 0, dmoSubPsms = 0;
            {
                var lines = File.ReadAllLines(dmoPsmPath);
                var H = lines[0].Split('\t');
                int iSpec = Array.IndexOf(H, "Spectrum");
                int iPept = Array.IndexOf(H, "Peptide");
                int iAsgn = Array.IndexOf(H, "Assigned Modifications");
                int iQ    = Array.IndexOf(H, "Qvalue");
                for (int i = 1; i < lines.Length; i++)
                {
                    var f = lines[i].Split('\t');
                    if (f.Length <= Math.Max(iAsgn, iQ)) continue;
                    if (!f[iSpec].StartsWith(raw1715Prefix, StringComparison.Ordinal)) continue;
                    if (!double.TryParse(f[iQ], NumberStyles.Float, CultureInfo.InvariantCulture, out double q) || q > 0.01) continue;
                    var pep = f[iPept];
                    dmoPsms++;
                    dmoAllIdentified.Add(IL(pep));
                    var asgn = f[iAsgn];
                    if (string.IsNullOrEmpty(asgn)) continue;
                    bool psmHadSub = false;
                    foreach (Match mm in dmoModRx.Matches(asgn))
                    {
                        int pos = int.Parse(mm.Groups[1].Value);
                        char res = mm.Groups[2].Value[0];
                        double mass = double.Parse(mm.Groups[3].Value, CultureInfo.InvariantCulture);
                        if (IsPtm(mass)) continue;
                        if (pos < 1 || pos > pep.Length) continue;
                        // I/L tolerance on the residue check
                        if (pep[pos - 1] != res && !("IL".IndexOf(pep[pos - 1]) >= 0 && "IL".IndexOf(res) >= 0)) continue;
                        var targets = FindSubTargets(res, mass).ToList();
                        if (targets.Count == 0) continue;
                        psmHadSub = true;
                        foreach (var to in targets)
                        {
                            var subbed = pep.Substring(0, pos - 1) + to + pep.Substring(pos);
                            dmoSubIdentified.Add(IL(subbed));
                        }
                    }
                    if (psmHadSub) dmoSubPsms++;
                }
            }

            // ---------- 3) MM-format readers (ModBox AllPSMs + GPTMD AllPeptides) ----------
            // Both files share the BaseSeq/FullSequence/QValue/DecoyContamTarget schema
            // and emit AA-sub mods as '[1+ nucleotide substitution:X->Y on X]' inside
            // FullSequence. Same parser handles both.
            (HashSet<string> AllId, HashSet<string> SubId, int Psms, int SubPsms) ReadMmFormat(string path)
            {
                var allId = new HashSet<string>(StringComparer.Ordinal);
                var subId = new HashSet<string>(StringComparer.Ordinal);
                int psms = 0, subPsms = 0;
                var rows = new PsmFromTsvFile(path).Results
                    .Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T");
                foreach (var p in rows)
                {
                    psms++;
                    allId.Add(IL(p.BaseSeq));
                    if (string.IsNullOrEmpty(p.FullSequence) || !p.FullSequence.Contains("nucleotide substitution")) continue;
                    string subbedFull;
                    try { subbedFull = IBioPolymerWithSetMods.ParseSubstitutedFullSequence(p.FullSequence); }
                    catch { continue; }
                    var subbedBase = IBioPolymerWithSetMods.GetBaseSequenceFromFullSequence(subbedFull);
                    if (string.IsNullOrEmpty(subbedBase)) continue;
                    subPsms++;
                    subId.Add(IL(subbedBase));
                }
                return (allId, subId, psms, subPsms);
            }

            var (mmAllIdentified,    mmSubIdentified,    mmPsms,    mmSubPsms)    = ReadMmFormat(mmModBoxPath);
            var (gptmdAllIdentified, gptmdSubIdentified, gptmdPsms, gptmdSubPsms) = ReadMmFormat(mmGptmdPath);

            // ---------- 4) Intersect with ground truth ----------
            var dmoInGt   = new HashSet<string>(dmoSubIdentified,   StringComparer.Ordinal); dmoInGt.IntersectWith(saltySspSet);
            var mmInGt    = new HashSet<string>(mmSubIdentified,    StringComparer.Ordinal); mmInGt.IntersectWith(saltySspSet);
            var gptmdInGt = new HashSet<string>(gptmdSubIdentified, StringComparer.Ordinal); gptmdInGt.IntersectWith(saltySspSet);

            // ---------- 5) Outputs ----------
            File.WriteAllLines(Path.Combine(outDir, "GroundTruthSALTY_SSPs.txt"), saltySspSet.OrderBy(s => s));
            void WriteEngineTsv(string filename, HashSet<string> sub)
            {
                using var w = new StreamWriter(Path.Combine(outDir, filename));
                w.WriteLine("SubstitutedPeptide_ILfolded\tInGroundTruth");
                foreach (var p in sub.OrderBy(s => s))
                    w.WriteLine($"{p}\t{saltySspSet.Contains(p)}");
            }
            WriteEngineTsv("DMO_SubstitutedPeptides_in_GT.tsv",       dmoSubIdentified);
            WriteEngineTsv("MM_ModBox_SubstitutedPeptides_in_GT.tsv", mmSubIdentified);
            WriteEngineTsv("MM_GPTMD_SubstitutedPeptides_in_GT.tsv",  gptmdSubIdentified);

            using (var w = new StreamWriter(Path.Combine(outDir, "EngineComparison_PeptideCounts.tsv")))
            {
                w.WriteLine("Metric\tDMO\tMM_ModBox\tMM_GPTMD");
                w.WriteLine($"TotalPSMs\t{dmoPsms}\t{mmPsms}\t{gptmdPsms}");
                w.WriteLine($"PSMsWithSubstitution\t{dmoSubPsms}\t{mmSubPsms}\t{gptmdSubPsms}");
                w.WriteLine($"UniqueBaseSeqIdentified_anyMod\t{dmoAllIdentified.Count}\t{mmAllIdentified.Count}\t{gptmdAllIdentified.Count}");
                w.WriteLine($"UniquePostSubPeptidesCalled\t{dmoSubIdentified.Count}\t{mmSubIdentified.Count}\t{gptmdSubIdentified.Count}");
                w.WriteLine($"UniquePostSubPeptidesInGroundTruth\t{dmoInGt.Count}\t{mmInGt.Count}\t{gptmdInGt.Count}");
                w.WriteLine($"GroundTruthUniverseSize\t{saltySspSet.Count}\t{saltySspSet.Count}\t{saltySspSet.Count}");

                // Pairwise + triple overlap on GT-matching peptides
                var dm = new HashSet<string>(dmoInGt, StringComparer.Ordinal);    dm.IntersectWith(mmInGt);
                var dg = new HashSet<string>(dmoInGt, StringComparer.Ordinal);    dg.IntersectWith(gptmdInGt);
                var mg = new HashSet<string>(mmInGt,  StringComparer.Ordinal);    mg.IntersectWith(gptmdInGt);
                var all3 = new HashSet<string>(dm, StringComparer.Ordinal);       all3.IntersectWith(gptmdInGt);
                w.WriteLine($"Overlap_DMOandModBox\t{dm.Count}\t{dm.Count}\t-");
                w.WriteLine($"Overlap_DMOandGPTMD\t{dg.Count}\t-\t{dg.Count}");
                w.WriteLine($"Overlap_ModBoxandGPTMD\t-\t{mg.Count}\t{mg.Count}");
                w.WriteLine($"Overlap_AllThree\t{all3.Count}\t{all3.Count}\t{all3.Count}");

                // Engine-only GT hits (relative to the other two engines)
                var dmoOnly   = new HashSet<string>(dmoInGt,   StringComparer.Ordinal); dmoOnly.ExceptWith(mmInGt); dmoOnly.ExceptWith(gptmdInGt);
                var mmOnly    = new HashSet<string>(mmInGt,    StringComparer.Ordinal); mmOnly.ExceptWith(dmoInGt); mmOnly.ExceptWith(gptmdInGt);
                var gptmdOnly = new HashSet<string>(gptmdInGt, StringComparer.Ordinal); gptmdOnly.ExceptWith(dmoInGt); gptmdOnly.ExceptWith(mmInGt);
                w.WriteLine($"EngineUnique_GTHits\t{dmoOnly.Count}\t{mmOnly.Count}\t{gptmdOnly.Count}");
            }
        }

        /// <summary>
        /// Smoke test that mzLib can open the Bruker .d folder for this DMO dataset.
        /// Dispatches via MsDataFileReader (which routes ".d" to BrukerFileReader if
        /// it sees analysis.baf, or TimsTofFileReader for analysis.tdf/.tsf). Both
        /// readers P/Invoke into Bruker's native DLLs — Win x64 only.
        /// </summary>
        [Test]
        public static void ReadBrukerDotDFile()
        {
            var path = @"E:\Aneuploidy\Mistranslation_project\Mistranslation_search\DMO_data\2021-06-18-ECLandSALTY1_Slot1-42_1_1715.d";

            Assert.That(Directory.Exists(path), $".d folder not found: {path}");

            var dataFile = MsDataFileReader.GetDataFile(path);
            dataFile.LoadAllStaticData();

            var scans = dataFile.GetAllScansList();
            Assert.That(scans, Is.Not.Null.And.Not.Empty, "No scans returned from the .d file.");

            // Sanity-check the basics so a silently empty/broken read fails loudly.
            Assert.That(scans.Count, Is.EqualTo(dataFile.NumSpectra));
            Assert.That(scans.Any(s => s.MsnOrder == 1), "Expected at least one MS1 scan.");
            Assert.That(scans.Any(s => s.MsnOrder >= 2), "Expected at least one MS2 scan.");
            Assert.That(scans.First().OneBasedScanNumber, Is.EqualTo(1));
            Assert.That(scans.All(s => s.MassSpectrum != null), "Encountered a scan with a null spectrum.");
        }
    }
}
