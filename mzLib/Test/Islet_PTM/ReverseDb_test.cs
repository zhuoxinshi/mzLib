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
using Omics.Modifications;

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
            var sspPath        = @"E:\Aneuploidy\Mistranslation_project\Mistranslation_search\DMO_data\DMO-FragPipe_ZenodoRepo\SSP_ECOLI_to_SALTY.csv";                                              // FindSSP.py output — defines the ground-truth SALTY SSP peptides
            var dmoPsmPath     = @"E:\Aneuploidy\Mistranslation_project\Mistranslation_search\DMO_data\DMO-FragPipe_ZenodoRepo\Results\AASubs_DMO\psm.tsv";                                          // FragPipe Detailed Mass Offset PSM table (covers both 1715 and 1716)
            var mmModBoxPath   = @"E:\Aneuploidy\Mistranslation_project\Mistranslation_search\DMO_data\MM\ModBox\1715_max1_AllAAsubPTM\Task\AllPSMs.psmtsv";                                          // MetaMorpheus ModBox search of 1715 only (AllPSMs = one row per PSM)
            var mmGptmdPath    = @"E:\Aneuploidy\Mistranslation_project\Mistranslation_search\DMO_data\MM\GPTMD\AllAAsubs_1715_noFilter\Task2-SearchTask\AllPeptides.psmtsv";                          // MetaMorpheus GPTMD post-search (AllPeptides = one row per unique peptide)
            var outDir         = @"E:\Aneuploidy\Mistranslation_project\Mistranslation_search\DMO_data\MM\ModBox\1715_max1_AllAAsubPTM\Task";                                                          // Where all comparison TSVs land
            const string raw1715Prefix = "2021-06-18-ECLandSALTY1_Slot1-42_1_1715";                                                                                                                    // Spectrum-ID prefix used to filter the DMO file down to 1715 only (MM only has 1715)
            const double massTol = 0.005;                                                                                                                                                              // ±5 mDa: tight enough to discriminate, loose enough for Orbitrap rounding
            const string allAa = "ACDEFGHIKLMNPQRSTVWY";                                                                                                                                               // 20 canonical amino acids — enumerate every from/to substitution pair

            static string IL(string s) => s?.Replace('I', 'L');                                                                                                                                        // I/L are isobaric and MS2-indistinguishable; fold both to L for set comparison

            // Build (from_aa, [to_aa, massShift]) lookup over ALL 20 AAs. The SSP
            // ground truth reflects real genome differences, not 1-nt-restricted
            // substitutions, so we don't apply the C/K/R rule here.
            var subsByFrom = new Dictionary<char, List<(char to, double shift)>>();                                                                                                                    // For each source residue, the list of (target residue, mass-shift) it could become
            foreach (var fr in allAa)                                                                                                                                                                  // 20 possible source residues
            {
                var list = new List<(char, double)>();                                                                                                                                                 // Buffer for this source residue's targets
                foreach (var to in allAa)                                                                                                                                                              // 20 possible target residues
                {
                    if (to == fr) continue;                                                                                                                                                            // Skip identity (not a substitution)
                    list.Add((to, Residue.ResidueMonoisotopicMass[to] - Residue.ResidueMonoisotopicMass[fr]));                                                                                         // Δmass = mass(to) − mass(from), can be positive or negative
                }
                subsByFrom[fr] = list;                                                                                                                                                                 // 19 entries per source residue → 380 total pairs
            }
            double[] ptmMasses = { 15.9949, 57.0214, 42.0106, 14.0157, 79.9663 };                                                                                                                      // Met-ox / Carbamidomethyl / acetyl / methyl / phospho — common PTMs to exclude from sub-shift interpretation
            bool IsPtm(double m) => ptmMasses.Any(p => Math.Abs(Math.Abs(m) - p) < massTol);                                                                                                           // |mass| match handles both +PTM and -PTM (neutral-loss) forms
            IEnumerable<char> FindSubTargets(char fromAa, double shift)                                                                                                                                 // Given source residue + mass shift, yield every target AA whose Δmass matches within tolerance
            {
                if (!subsByFrom.TryGetValue(fromAa, out var list)) yield break;                                                                                                                        // Source residue not in lookup → no candidates
                foreach (var (to, s) in list)                                                                                                                                                          // Walk all 19 candidates for this source
                    if (Math.Abs(s - shift) < massTol) yield return to;                                                                                                                                // Return every target whose mass shift matches (may be 0, 1, or several — isobaric subs are rare but happen)
            }

            // ---------- 1) Ground truth: SALTY SSP peptide SET ----------
            // CSV row: <idx>,<BaseSeq>,[<SSP list>],<IsSSP>,[<DSP list>],<IsDSP>,<IsExact>
            // The bracketed list values contain commas — extract via regex on [...] groups.
            var saltySspSet = new HashSet<string>(StringComparer.Ordinal);                                                                                                                              // Union of every SALTY peptide that's an SSP of some E. coli peptide
            var bracketRx   = new Regex(@"\[([^\]]*)\]", RegexOptions.Compiled);                                                                                                                       // Matches the two Python-list columns in the SSP CSV
            var quotedRx    = new Regex(@"'([A-Z]+)'", RegexOptions.Compiled);                                                                                                                         // Each quoted entry inside a list is one peptide string
            foreach (var line in File.ReadLines(sspPath).Skip(1))                                                                                                                                       // Stream the CSV, skip the header row
            {
                if (string.IsNullOrWhiteSpace(line)) continue;                                                                                                                                         // Skip blank lines
                var brackets = bracketRx.Matches(line);                                                                                                                                                // First bracket = SSP list, second = DSP list
                if (brackets.Count == 0) continue;                                                                                                                                                     // Malformed row
                foreach (Match m in quotedRx.Matches(brackets[0].Groups[1].Value))                                                                                                                     // Only iterate the SSP bracket (index 0); DSP is out of scope
                    saltySspSet.Add(IL(m.Groups[1].Value));                                                                                                                                            // Add the SALTY peptide string (I/L folded) to the GT set
            }

            // ---------- 2) DMO ----------
            // 'Assigned Modifications' entries look like '7E(-14.0157), 12D(32.0415)'.
            var dmoModRx = new Regex(@"(\d+)([A-Z])\(([-\d.]+)\)", RegexOptions.Compiled);                                                                                                              // Captures (position, residue, mass-shift) for each assigned mod entry
            var dmoAllIdentified = new HashSet<string>(StringComparer.Ordinal);                                                                                                                        // Set of all unique BaseSeq DMO returned (any mod state)
            var dmoSubIdentified = new HashSet<string>(StringComparer.Ordinal);                                                                                                                        // Set of unique POST-substitution peptide strings DMO called
            int dmoPsms = 0, dmoSubPsms = 0;                                                                                                                                                           // Total PSM count after Q-filter, plus subset carrying ≥1 substitution
            {
                var lines = File.ReadAllLines(dmoPsmPath);                                                                                                                                             // Tab-separated PSM table; small enough to slurp
                var H = lines[0].Split('\t');                                                                                                                                                          // Header row
                int iSpec = Array.IndexOf(H, "Spectrum");                                                                                                                                              // Column with raw-file-prefixed scan identifier
                int iPept = Array.IndexOf(H, "Peptide");                                                                                                                                               // Base peptide (no mods)
                int iAsgn = Array.IndexOf(H, "Assigned Modifications");                                                                                                                                // Where substitution mass shifts live
                int iQ    = Array.IndexOf(H, "Qvalue");                                                                                                                                                // Philosopher-reported FDR
                for (int i = 1; i < lines.Length; i++)                                                                                                                                                 // Data rows start at index 1
                {
                    var f = lines[i].Split('\t');                                                                                                                                                      // Tokenize this row
                    if (f.Length <= Math.Max(iAsgn, iQ)) continue;                                                                                                                                     // Truncated row — skip
                    if (!f[iSpec].StartsWith(raw1715Prefix, StringComparison.Ordinal)) continue;                                                                                                       // Drop 1716 PSMs — MM only searched 1715, keep comparison apples-to-apples
                    if (!double.TryParse(f[iQ], NumberStyles.Float, CultureInfo.InvariantCulture, out double q) || q > 0.01) continue;                                                                 // 1% FDR cutoff (same as MM)
                    var pep = f[iPept];                                                                                                                                                                // E. coli base sequence as MM identified it
                    dmoPsms++;                                                                                                                                                                         // PSM passed all filters
                    dmoAllIdentified.Add(IL(pep));                                                                                                                                                     // Add to "any mod state" set (used as a coverage denominator)
                    var asgn = f[iAsgn];                                                                                                                                                               // The assigned-mods field, possibly empty
                    if (string.IsNullOrEmpty(asgn)) continue;                                                                                                                                          // No mods → not a substitution PSM
                    bool psmHadSub = false;                                                                                                                                                            // Did at least one entry decode to a real substitution?
                    foreach (Match mm in dmoModRx.Matches(asgn))                                                                                                                                       // Could be 1 or more "<pos><res>(<mass>)" entries
                    {
                        int pos = int.Parse(mm.Groups[1].Value);                                                                                                                                       // 1-based residue position
                        char res = mm.Groups[2].Value[0];                                                                                                                                              // Residue letter at that position
                        double mass = double.Parse(mm.Groups[3].Value, CultureInfo.InvariantCulture);                                                                                                  // Mass shift in Da
                        if (IsPtm(mass)) continue;                                                                                                                                                     // Skip common-PTM shifts so we don't mis-call Met-ox/Carbam/etc. as substitutions
                        if (pos < 1 || pos > pep.Length) continue;                                                                                                                                     // Out-of-bounds position — bad parse
                        // I/L tolerance on the residue check
                        if (pep[pos - 1] != res && !("IL".IndexOf(pep[pos - 1]) >= 0 && "IL".IndexOf(res) >= 0)) continue;                                                                              // Bail unless the peptide residue at that position matches the mod's residue (or they're both I/L)
                        var targets = FindSubTargets(res, mass).ToList();                                                                                                                              // What target residues does this (residue, mass-shift) decode to?
                        if (targets.Count == 0) continue;                                                                                                                                              // No matching AA pair → not a substitution shift (probably unusual PTM)
                        psmHadSub = true;                                                                                                                                                              // We commit this PSM to the substitution count
                        foreach (var to in targets)                                                                                                                                                    // For each candidate target (usually exactly 1, occasionally several)
                        {
                            var subbed = pep.Substring(0, pos - 1) + to + pep.Substring(pos);                                                                                                          // Splice in the new residue → post-substitution base sequence
                            dmoSubIdentified.Add(IL(subbed));                                                                                                                                          // Add to the unique-substituted-peptide set (I/L folded)
                        }
                    }
                    if (psmHadSub) dmoSubPsms++;                                                                                                                                                       // Count the PSM once if any of its mod entries was a substitution
                }
            }

            // ---------- 3) MM-format readers (ModBox AllPSMs + GPTMD AllPeptides) ----------
            // Both files share the BaseSeq/FullSequence/QValue/DecoyContamTarget schema
            // and emit AA-sub mods as '[1+ nucleotide substitution:X->Y on X]' inside
            // FullSequence. Same parser handles both.
            (HashSet<string> AllId, HashSet<string> SubId, int Psms, int SubPsms) ReadMmFormat(string path)                                                                                            // Local function returns four parallel measures per engine
            {
                var allId = new HashSet<string>(StringComparer.Ordinal);                                                                                                                               // Unique BaseSeq across all PSMs (denominator)
                var subId = new HashSet<string>(StringComparer.Ordinal);                                                                                                                               // Unique post-substitution peptide strings (the numerator candidates)
                int psms = 0, subPsms = 0;                                                                                                                                                             // PSM counts
                var rows = new PsmFromTsvFile(path).Results                                                                                                                                            // mzLib's MM PSM reader — handles either AllPSMs or AllPeptides schema
                    .Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T");                                                                                                                       // Same 1% FDR + target-only filter as DMO
                foreach (var p in rows)
                {
                    psms++;                                                                                                                                                                            // PSM count after filter
                    allId.Add(IL(p.BaseSeq));                                                                                                                                                          // Coverage denominator
                    if (string.IsNullOrEmpty(p.FullSequence) || !p.FullSequence.Contains("nucleotide substitution")) continue;                                                                         // No sub annotation → skip
                    string subbedFull;                                                                                                                                                                 // Will hold the FullSequence with all AA-sub annotations applied
                    try { subbedFull = IBioPolymerWithSetMods.ParseSubstitutedFullSequence(p.FullSequence); }                                                                                          // Canonical mzLib parser — handles MM's "1+ nucleotide substitution" form
                    catch { continue; }                                                                                                                                                                // Malformed annotation — skip rather than crash
                    var subbedBase = IBioPolymerWithSetMods.GetBaseSequenceFromFullSequence(subbedFull);                                                                                               // Strip mod brackets → bare AA string
                    if (string.IsNullOrEmpty(subbedBase)) continue;                                                                                                                                    // Empty result — skip
                    subPsms++;                                                                                                                                                                         // PSM count of substitution-bearing rows
                    subId.Add(IL(subbedBase));                                                                                                                                                         // Add I/L-folded post-sub peptide to the unique set
                }
                return (allId, subId, psms, subPsms);                                                                                                                                                  // Hand back the four measures
            }

            var (mmAllIdentified,    mmSubIdentified,    mmPsms,    mmSubPsms)    = ReadMmFormat(mmModBoxPath);                                                                                        // Apply the reader to MM ModBox
            var (gptmdAllIdentified, gptmdSubIdentified, gptmdPsms, gptmdSubPsms) = ReadMmFormat(mmGptmdPath);                                                                                         // Apply the reader to MM GPTMD

            // ---------- 4) Intersect with ground truth ----------
            var dmoInGt   = new HashSet<string>(dmoSubIdentified,   StringComparer.Ordinal); dmoInGt.IntersectWith(saltySspSet);                                                                       // DMO's post-sub peptides that the GT set contains
            var mmInGt    = new HashSet<string>(mmSubIdentified,    StringComparer.Ordinal); mmInGt.IntersectWith(saltySspSet);                                                                        // Same for MM ModBox
            var gptmdInGt = new HashSet<string>(gptmdSubIdentified, StringComparer.Ordinal); gptmdInGt.IntersectWith(saltySspSet);                                                                     // Same for MM GPTMD

            // ---------- 5) Outputs ----------
            File.WriteAllLines(Path.Combine(outDir, "GroundTruthSALTY_SSPs.txt"), saltySspSet.OrderBy(s => s));                                                                                        // Dump the GT set, one peptide per line, for reference
            void WriteEngineTsv(string filename, HashSet<string> sub)                                                                                                                                  // Per-engine "called set" TSV with InGT True/False flag
            {
                using var w = new StreamWriter(Path.Combine(outDir, filename));                                                                                                                        // Auto-flush + dispose
                w.WriteLine("SubstitutedPeptide_ILfolded\tInGroundTruth");                                                                                                                             // Header
                foreach (var p in sub.OrderBy(s => s))                                                                                                                                                 // Sorted alphabetically for easy diff/inspection
                    w.WriteLine($"{p}\t{saltySspSet.Contains(p)}");                                                                                                                                    // Each unique called sub peptide + whether GT contains it
            }
            WriteEngineTsv("DMO_SubstitutedPeptides_in_GT.tsv",       dmoSubIdentified);                                                                                                               // DMO output
            WriteEngineTsv("MM_ModBox_SubstitutedPeptides_in_GT.tsv", mmSubIdentified);                                                                                                                // MM ModBox output
            WriteEngineTsv("MM_GPTMD_SubstitutedPeptides_in_GT.tsv",  gptmdSubIdentified);                                                                                                             // MM GPTMD output

            using (var w = new StreamWriter(Path.Combine(outDir, "EngineComparison_PeptideCounts.tsv")))                                                                                                // Summary table; columns = engines, rows = metrics
            {
                w.WriteLine("Metric\tDMO\tMM_ModBox\tMM_GPTMD");                                                                                                                                       // Header
                w.WriteLine($"TotalPSMs\t{dmoPsms}\t{mmPsms}\t{gptmdPsms}");                                                                                                                           // PSMs after Q-filter (any mod state)
                w.WriteLine($"PSMsWithSubstitution\t{dmoSubPsms}\t{mmSubPsms}\t{gptmdSubPsms}");                                                                                                       // PSMs that carry a substitution call
                w.WriteLine($"UniqueBaseSeqIdentified_anyMod\t{dmoAllIdentified.Count}\t{mmAllIdentified.Count}\t{gptmdAllIdentified.Count}");                                                         // Coverage of the proteome (denominator-like)
                w.WriteLine($"UniquePostSubPeptidesCalled\t{dmoSubIdentified.Count}\t{mmSubIdentified.Count}\t{gptmdSubIdentified.Count}");                                                            // Unique substituted-peptide STRINGS — the main per-engine "called" set
                w.WriteLine($"UniquePostSubPeptidesInGroundTruth\t{dmoInGt.Count}\t{mmInGt.Count}\t{gptmdInGt.Count}");                                                                                // The headline metric — recovered GT SSPs
                w.WriteLine($"GroundTruthUniverseSize\t{saltySspSet.Count}\t{saltySspSet.Count}\t{saltySspSet.Count}");                                                                                // Total achievable — same on every column

                // Pairwise + triple overlap on GT-matching peptides
                var dm = new HashSet<string>(dmoInGt, StringComparer.Ordinal);    dm.IntersectWith(mmInGt);                                                                                            // DMO ∩ MM(ModBox)
                var dg = new HashSet<string>(dmoInGt, StringComparer.Ordinal);    dg.IntersectWith(gptmdInGt);                                                                                         // DMO ∩ MM(GPTMD)
                var mg = new HashSet<string>(mmInGt,  StringComparer.Ordinal);    mg.IntersectWith(gptmdInGt);                                                                                         // MM(ModBox) ∩ MM(GPTMD)
                var all3 = new HashSet<string>(dm, StringComparer.Ordinal);       all3.IntersectWith(gptmdInGt);                                                                                       // Triple overlap = high-confidence consensus call set
                w.WriteLine($"Overlap_DMOandModBox\t{dm.Count}\t{dm.Count}\t-");                                                                                                                       // Pair counts written twice so the column under each engine matches
                w.WriteLine($"Overlap_DMOandGPTMD\t{dg.Count}\t-\t{dg.Count}");
                w.WriteLine($"Overlap_ModBoxandGPTMD\t-\t{mg.Count}\t{mg.Count}");
                w.WriteLine($"Overlap_AllThree\t{all3.Count}\t{all3.Count}\t{all3.Count}");                                                                                                            // Triple overlap broadcast across all three columns

                // Engine-only GT hits (relative to the other two engines)
                var dmoOnly   = new HashSet<string>(dmoInGt,   StringComparer.Ordinal); dmoOnly.ExceptWith(mmInGt); dmoOnly.ExceptWith(gptmdInGt);                                                     // GT hits only DMO found
                var mmOnly    = new HashSet<string>(mmInGt,    StringComparer.Ordinal); mmOnly.ExceptWith(dmoInGt); mmOnly.ExceptWith(gptmdInGt);                                                      // GT hits only MM ModBox found
                var gptmdOnly = new HashSet<string>(gptmdInGt, StringComparer.Ordinal); gptmdOnly.ExceptWith(dmoInGt); gptmdOnly.ExceptWith(mmInGt);                                                   // GT hits only MM GPTMD found
                w.WriteLine($"EngineUnique_GTHits\t{dmoOnly.Count}\t{mmOnly.Count}\t{gptmdOnly.Count}");                                                                                               // Three engine-unique counts on one row
            }
        }

        /// <summary>
        /// For each MetaMorpheus PSM carrying an AA-substitution annotation
        /// (<c>[1+ nucleotide substitution:X->Y on X]</c> inside FullSequence), check
        /// whether the substitution's residue mass shift Δm = mass(Y) − mass(X) could
        /// instead be explained by a user-specified common PTM applied to residue X.
        ///
        /// Following Mordret et al. (Mol. Cell 2019, "Systematic Detection of Amino
        /// Acid Substitutions in Proteomes ...", https://www.sciencedirect.com/science/article/pii/S1097276519304988),
        /// a PTM "explains" the shift only when BOTH:
        ///   1) the PTM's monoisotopic mass matches Δm within tolerance, AND
        ///   2) the PTM's motif/site (Target uppercase residue) and LocationRestriction
        ///      are compatible with the substitution site — e.g. Phosphorylation
        ///      (+79.9663) on S/T only flags S→? and T→? subs whose mass shift lands
        ///      at +79.9663; an N-terminal-restricted PTM only flags subs at peptide
        ///      (or protein) N-terminus.
        ///
        /// PTMs are pulled from mzLib's bundled Mods.txt by ID (caller-supplied list).
        ///
        /// When an SSP ground-truth CSV is provided (FindSSP.py output —
        /// SSP_ECOLI_to_SALTY.csv for the 1715 ECLandSALTY co-lysate), the test also
        /// reports per-substitution TP/FP counts BEFORE and AFTER the filter. The
        /// "true substitution" set is the union of every SALTY SSP sequence; a
        /// substitution annotation is a TP if its singly-substituted peptide
        /// (I/L folded) is in that set, else a FP ("nonSSP" call, in the paper's
        /// terms). This mirrors the false-positive comparison in Mordret et al.
        ///
        /// Writes next to the PSM file:
        ///   <c>{stem}_AaSubFiltered_kept.tsv</c>      — substitutions that survived
        ///   <c>{stem}_AaSubFiltered_dropped.tsv</c>   — substitutions explained by a common PTM
        ///   <c>{stem}_AaSubFiltered_summary.tsv</c>   — TP/FP counts before vs after
        /// </summary>
        [Test]
        public static void FilterMmAaSubsByCommonPtmMassAndMotif()
        {
            // ---- Inputs --------------------------------------------------------
            var mmPsmPath = @"E:\Aneuploidy\Mistranslation_project\Mistranslation_search\DMO_data\MM\ModBox\1715_max1_AllAAsubPTM\Task\AllPSMs.psmtsv";
            // Optional — set to "" to skip the TP/FP comparison block. FindSSP.py
            // output defines the SALTY SSP peptides we treat as the "true" sub set.
            var sspGroundTruthCsv = @"E:\Aneuploidy\Mistranslation_project\Mistranslation_search\DMO_data\DMO-FragPipe_ZenodoRepo\SSP_ECOLI_to_SALTY.csv";
            // Caller-supplied "common PTMs" by Mods.txt ID. Anything here, when
            // matched on BOTH mass AND motif, is treated as a confounding PTM and
            // the corresponding AA-sub call is removed.
            var commonPtmIds = new HashSet<string>(StringComparer.OrdinalIgnoreCase)
            {
                "Phosphorylation",
                "Oxidation on M",
                "Acetylation",
                "Methylation",
                "Carbamidomethyl",
            };
            // mzLib ships a Mods.txt with Anywhere/N-term/C-term/Protein-term entries already annotated — reuse instead of hand-rolling a PTM table.
            var modsTxtPath = Path.Combine(
                TestContext.CurrentContext.TestDirectory, "..", "..", "..", "..",
                "Omics", "Resources", "Mods.txt");
            modsTxtPath = Path.GetFullPath(modsTxtPath);
            const double massTol = 0.005;

            // ---- Load the common PTM rows (one Modification per site/term entry) - A single PTM ID (e.g. "Phosphorylation") may expand to several rows
            // because Mods.txt splits by Target and LocationRestriction (S/T vs Y, Anywhere vs N-terminal, ...). Keep them all so each combination is tested independently.
            var allMods = PtmListLoader.ReadModsFromFile(modsTxtPath, out _).ToList();
            var commonPtms = allMods
                .Where(m => m.OriginalId != null && commonPtmIds.Contains(m.OriginalId)
                            && m.MonoisotopicMass.HasValue && m.Target != null)
                .ToList();
            TestContext.WriteLine($"Loaded {commonPtms.Count} common-PTM rows covering {commonPtms.Select(m => m.OriginalId).Distinct().Count()} distinct PTM IDs.");
            Assert.That(commonPtms, Is.Not.Empty, "No common PTMs matched the supplied IDs in Mods.txt.");

            // ---- Walk a FullSequence and extract every embedded sub annotation --
            // Returns (position1Based, originalAa, substitutedAa) for each [N+ nucleotide substitution:X->Y on Z] tag. The position is the index
            // of the most-recently-seen base-sequence residue at the time the bracket opens — that's the residue carrying the sub annotation.
            var subRx = new Regex(@"^\d+\+?\s*nucleotide substitution:\s*([A-Z])->([A-Z]) on ([A-Z])$", RegexOptions.Compiled);
            static IEnumerable<(int pos1, char from, char to)> ExtractSubs(string fullSeq, Regex rx)
            {
                if (string.IsNullOrEmpty(fullSeq)) yield break;
                int basePos = 0;
                for (int i = 0; i < fullSeq.Length; i++)
                {
                    char c = fullSeq[i];
                    if (c == '[')
                    {
                        int close = fullSeq.IndexOf(']', i + 1);
                        if (close < 0) yield break;
                        var content = fullSeq.Substring(i + 1, close - i - 1);
                        var m = rx.Match(content);
                        if (m.Success && basePos >= 1)
                            yield return (basePos, m.Groups[1].Value[0], m.Groups[2].Value[0]);
                        i = close;
                    }
                    else if (c >= 'A' && c <= 'Z')
                    {
                        basePos++;
                    }
                    // dashes, digits, lowercase tags inside brackets are skipped
                }
            }

            // ---- PTM compatibility predicate (mass + motif + location) ----------
            // peptideLength is needed to test C-terminal restrictions; pepIsProtNterm
            // / pepIsProtCterm distinguish "Protein N-terminal" from peptide N-term.
            static bool PtmCompatible(Modification ptm, double deltaMass, char originalAa,
                int subPos1Based, int peptideLength,
                bool pepIsProtNterm, bool pepIsProtCterm, double tol)
            {
                if (!ptm.MonoisotopicMass.HasValue) return false;
                if (Math.Abs(ptm.MonoisotopicMass.Value - deltaMass) > tol) return false;

                // Motif: the only uppercase letter in the Target string is the AA
                // the PTM targets. The substitution's original residue must equal it.
                var motif = ptm.Target?.ToString();
                if (string.IsNullOrEmpty(motif)) return false;
                char target = motif.First(char.IsUpper);
                // 'X' in Mods.txt motifs means "any residue" (used for terminal mods
                // like N-terminal acetylation that don't care about identity).
                if (target != 'X' && target != originalAa) return false;

                // Location restriction. Mods.txt uses trailing period — normalize.
                var loc = (ptm.LocationRestriction ?? "Anywhere.").TrimEnd('.');
                switch (loc)
                {
                    case "Anywhere": return true;
                    case "N-terminal": return subPos1Based == 1;
                    case "C-terminal": return subPos1Based == peptideLength;
                    case "Peptide N-terminal": return subPos1Based == 1;
                    case "Peptide C-terminal": return subPos1Based == peptideLength;
                    case "Protein N-terminal": return subPos1Based == 1 && pepIsProtNterm;
                    case "Protein C-terminal": return subPos1Based == peptideLength && pepIsProtCterm;
                    default: return false; // Unassigned. or unknown → don't filter
                }
            }

            // ---- Iterate PSMs --------------------------------------------------
            var psms = new PsmFromTsvFile(mmPsmPath).Results
                .Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T"
                         && p.FullSequence != null && p.FullSequence.Contains("nucleotide substitution"))
                .ToList();

            var outDir = Path.GetDirectoryName(mmPsmPath) ?? ".";
            var stem = Path.GetFileNameWithoutExtension(mmPsmPath);
            var keptPath = Path.Combine(outDir, $"{stem}_AaSubFiltered_kept.tsv");
            var droppedPath = Path.Combine(outDir, $"{stem}_AaSubFiltered_dropped.tsv");
            var summaryPath = Path.Combine(outDir, $"{stem}_AaSubFiltered_summary.tsv");

            // ---- Optional: load SSP ground truth for TP/FP comparison ---------
            // Row layout: <idx>,<EcoliBaseSeq>,[<SSP list>],<IsSSP>,[<DSP list>],<IsDSP>,<IsExact>
            // The bracket lists contain commas, so we regex out the [...] groups.
            // I/L is folded everywhere — same convention as the existing SSP tests.
            static string IL(string s) => s?.Replace('I', 'L');
            HashSet<string> saltySspSet = null;
            if (!string.IsNullOrEmpty(sspGroundTruthCsv) && File.Exists(sspGroundTruthCsv))
            {
                saltySspSet = new HashSet<string>(StringComparer.Ordinal);
                var bracketRx = new Regex(@"\[([^\]]*)\]", RegexOptions.Compiled);
                var quotedRx = new Regex(@"'([A-Z]+)'", RegexOptions.Compiled);
                foreach (var line in File.ReadLines(sspGroundTruthCsv).Skip(1))
                {
                    if (string.IsNullOrWhiteSpace(line)) continue;
                    var brackets = bracketRx.Matches(line);
                    if (brackets.Count == 0) continue;
                    foreach (Match m in quotedRx.Matches(brackets[0].Groups[1].Value))
                        saltySspSet.Add(IL(m.Groups[1].Value));
                }
                TestContext.WriteLine($"Loaded {saltySspSet.Count} unique SALTY SSP peptides from {sspGroundTruthCsv}");
            }
            else
            {
                TestContext.WriteLine("No SSP ground-truth CSV provided — skipping TP/FP block.");
            }

            // For protein-terminus checks, parse "1 to N" out of StartAndEndResiduesInProtein.
            var posRx = new Regex(@"(\d+)\s+to\s+(\d+)", RegexOptions.Compiled);

            // Counters split four ways: kept vs dropped, GT-hit vs GT-miss. Only the
            // last two are meaningful when the SSP set is loaded. "TP" here means
            // the singly-substituted peptide for this annotation is in the SALTY
            // SSP set; "FP" (= the paper's "nonSSP") means it isn't.
            int kept = 0, dropped = 0;
            int keptTp = 0, keptFp = 0, droppedTp = 0, droppedFp = 0;
            using var keptWriter = new StreamWriter(keptPath);
            using var dropWriter = new StreamWriter(droppedPath);
            string tpCol = saltySspSet != null ? "\tInSspGroundTruth" : "";
            keptWriter.WriteLine("BaseSeq\tFullSequence\tSubPosInPeptide\tFromAa\tToAa\tDeltaMass\tSubstitutedPeptide\tScore\tQValue\tProteinAccession\tStartAndEnd" + tpCol);
            dropWriter.WriteLine("BaseSeq\tFullSequence\tSubPosInPeptide\tFromAa\tToAa\tDeltaMass\tSubstitutedPeptide\tMatchedPtms\tScore\tQValue\tProteinAccession\tStartAndEnd" + tpCol);

            foreach (var psm in psms)
            {
                var subs = ExtractSubs(psm.FullSequence, subRx).ToList();
                if (subs.Count == 0) continue;

                // Pre-compute peptide context for terminus restrictions.
                int pepLen = psm.BaseSeq?.Length ?? 0;
                bool pepIsProtNterm = false, pepIsProtCterm = false;
                var pm = posRx.Match(psm.StartAndEndResiduesInProtein ?? "");
                if (pm.Success)
                {
                    pepIsProtNterm = int.Parse(pm.Groups[1].Value) == 1;
                    // We don't have protein length here, so leave protein-C-term false
                    // unless explicitly stated; conservative — protein-C-term-restricted
                    // PTMs simply won't fire, which is the safer default.
                    pepIsProtCterm = false;
                }

                foreach (var (pos1, from, to) in subs)
                {
                    double deltaMass = Residue.ResidueMonoisotopicMass[(int)to]
                                     - Residue.ResidueMonoisotopicMass[(int)from];

                    // Singly-substituted peptide: splice ONLY this sub into the
                    // original base seq. Used both for the ground-truth check and
                    // for the output column (so the user can audit calls by eye).
                    string subbedPep = "";
                    if (!string.IsNullOrEmpty(psm.BaseSeq) && pos1 >= 1 && pos1 <= psm.BaseSeq.Length)
                        subbedPep = psm.BaseSeq.Substring(0, pos1 - 1) + to + psm.BaseSeq.Substring(pos1);

                    bool? inGt = saltySspSet != null && !string.IsNullOrEmpty(subbedPep)
                        ? saltySspSet.Contains(IL(subbedPep))
                        : (bool?)null;
                    string inGtCol = inGt.HasValue ? $"\t{inGt.Value}" : "";

                    var matched = commonPtms
                        .Where(m => PtmCompatible(m, deltaMass, from, pos1, pepLen,
                                                  pepIsProtNterm, pepIsProtCterm, massTol))
                        .Select(m => $"{m.OriginalId}({m.LocationRestriction}/{m.Target})")
                        .Distinct()
                        .ToList();

                    if (matched.Count > 0)
                    {
                        dropWriter.WriteLine(string.Join("\t",
                            psm.BaseSeq, psm.FullSequence, pos1, from, to,
                            deltaMass.ToString("F4", CultureInfo.InvariantCulture),
                            subbedPep,
                            string.Join(";", matched),
                            psm.Score.ToString("F2", CultureInfo.InvariantCulture),
                            psm.QValue.ToString("G", CultureInfo.InvariantCulture),
                            psm.ProteinAccession ?? "",
                            psm.StartAndEndResiduesInProtein ?? "") + inGtCol);
                        dropped++;
                        if (inGt == true) droppedTp++;
                        else if (inGt == false) droppedFp++;
                    }
                    else
                    {
                        keptWriter.WriteLine(string.Join("\t",
                            psm.BaseSeq, psm.FullSequence, pos1, from, to,
                            deltaMass.ToString("F4", CultureInfo.InvariantCulture),
                            subbedPep,
                            psm.Score.ToString("F2", CultureInfo.InvariantCulture),
                            psm.QValue.ToString("G", CultureInfo.InvariantCulture),
                            psm.ProteinAccession ?? "",
                            psm.StartAndEndResiduesInProtein ?? "") + inGtCol);
                        kept++;
                        if (inGt == true) keptTp++;
                        else if (inGt == false) keptFp++;
                    }
                }
            }

            TestContext.WriteLine($"Substitution annotations evaluated: {kept + dropped} (kept={kept}, dropped={dropped})");
            TestContext.WriteLine($"  kept    → {keptPath}");
            TestContext.WriteLine($"  dropped → {droppedPath}");

            // ---- TP/FP summary (only when ground-truth was loaded) -------------
            if (saltySspSet != null)
            {
                int preTp = keptTp + droppedTp;
                int preFp = keptFp + droppedFp;
                int preTotal = preTp + preFp;
                int postTp = keptTp;
                int postFp = keptFp;
                int postTotal = postTp + postFp;

                static double Pct(int num, int den) => den == 0 ? 0 : 100.0 * num / den;

                using (var w = new StreamWriter(summaryPath))
                {
                    w.WriteLine("Metric\tPreFilter\tPostFilter\tDelta");
                    w.WriteLine($"TotalSubstitutions\t{preTotal}\t{postTotal}\t{postTotal - preTotal}");
                    w.WriteLine($"TPs_inSSP\t{preTp}\t{postTp}\t{postTp - preTp}");
                    w.WriteLine($"FPs_nonSSP\t{preFp}\t{postFp}\t{postFp - preFp}");
                    w.WriteLine($"FPRate_pct\t{Pct(preFp, preTotal):F2}\t{Pct(postFp, postTotal):F2}\t{Pct(postFp, postTotal) - Pct(preFp, preTotal):F2}");
                    w.WriteLine($"TPRate_pct\t{Pct(preTp, preTotal):F2}\t{Pct(postTp, postTotal):F2}\t{Pct(postTp, postTotal) - Pct(preTp, preTotal):F2}");
                    w.WriteLine();
                    w.WriteLine("FilterAccountingMetric\tCount\tPctOfDropped");
                    w.WriteLine($"Dropped_TPs(true subs lost)\t{droppedTp}\t{Pct(droppedTp, dropped):F2}");
                    w.WriteLine($"Dropped_FPs(common-PTM masquerades removed)\t{droppedFp}\t{Pct(droppedFp, dropped):F2}");
                    w.WriteLine($"FilterPrecision_pct(Dropped_FPs/Dropped_All)\t{Pct(droppedFp, dropped):F2}\t");
                }

                TestContext.WriteLine("TP/FP comparison vs SALTY SSP ground truth:");
                TestContext.WriteLine($"  Pre-filter : total={preTotal}  TP={preTp}  FP={preFp}  FP%={Pct(preFp, preTotal):F2}");
                TestContext.WriteLine($"  Post-filter: total={postTotal}  TP={postTp}  FP={postFp}  FP%={Pct(postFp, postTotal):F2}");
                TestContext.WriteLine($"  Filter accounting: dropped {droppedFp} FPs and {droppedTp} TPs (precision {Pct(droppedFp, dropped):F2}%)");
                TestContext.WriteLine($"  summary → {summaryPath}");
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
