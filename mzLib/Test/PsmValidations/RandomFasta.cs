using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra.Factorization;
using NUnit.Framework;
using Omics.Modifications;
using Proteomics;
using Readers;
using UsefulProteomicsDatabases;

namespace Test
{
    public class RandomFasta
    {
        [Test]
        public static void WriteNewFasta()
        {
            var xmlDb = @"E:\Aneuploidy\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16.xml";
            var fastaDb = @"E:\Aneuploidy\uniprotkb_taxonomy_id_559292_AND_review_2024_10_02.fasta";
            var allProteins = ProteinDbLoader.LoadProteinXML(xmlDb, true, DecoyType.None, null, false,
                new List<string>(), out Dictionary<string, Modification> un);
            var psmPath = @"E:\DIA\Data\DIA_bu_250114\fasta\Task1-SearchTask\AllPeptides.psmtsv";
            var allPeptides = SpectrumMatchTsvReader.ReadPsmTsv(psmPath, out List<string> _).Where(p => p.DecoyContamTarget == "T" && p.QValue <= 0.01).ToList();
            var allPeptidesNoMod = allPeptides.Select(p => !p.FullSequence.Contains("[")).ToList();

            Random rand = new Random();
            int numberOfPeptideSeqToMutate = 100;
            int[] randomIndices = Enumerable.Range(0, allPeptidesNoMod.Count).OrderBy(_ => rand.Next()).Take(numberOfPeptideSeqToMutate).ToArray();
            var randomPeptides = randomIndices.Select(i => allPeptides[i]).ToList();
            
            foreach(var peptide in randomPeptides)
            {
                var protein = allProteins.FirstOrDefault(p => p.Accession == peptide.ProteinAccession);
                var mutatedSequence = RandomlyMutateAminoAcid(peptide.FullSequence);
                var pos = peptide.StartAndEndResiduesInProtein.Trim('[', ']').Split(new[] { " to " }, StringSplitOptions.None);
                var newProteinSequence = protein.BaseSequence.Substring(0, int.Parse(pos[0]) - 1) + mutatedSequence + protein.BaseSequence.Substring(int.Parse(pos[1]));
                var modifiedProtein = new Protein(protein, newProteinSequence);
            }
        }

        public static string RandomlyMutateAminoAcid(string sequence)
        {
            if (string.IsNullOrEmpty(sequence) || sequence.Length == 1) return sequence;

            const string aminoAcids = "ADEFGHIKLNPQRSTVWY";
            var rand = new Random();
            int pos = rand.Next(sequence.Length);

            char original = sequence[pos];
            char newAmino;
            do
            {
                newAmino = aminoAcids[rand.Next(aminoAcids.Length)];
            } while (newAmino == original);

            char[] arr = sequence.ToCharArray();
            arr[pos] = newAmino;
            return new string(arr);
        }

        [Test]
        public static void TestMakeNewProteins()
        {
            var xmlDb = @"E:\Aneuploidy\uniprotkb_taxonomy_id_559292_AND_review_2024_08_16.xml";
            var fastaDb = @"E:\Aneuploidy\uniprotkb_taxonomy_id_559292_AND_review_2024_10_02.fasta";
            var allProteins = ProteinDbLoader.LoadProteinXML(xmlDb, true, DecoyType.None, null, false,
                new List<string>(), out Dictionary<string, Modification> un);
            var psmPath = @"E:\DIA\Data\DIA_bu_250114\fasta\Task1-SearchTask\AllPeptides.psmtsv";
            var allPeptides = SpectrumMatchTsvReader.ReadPsmTsv(psmPath, out List<string> _).Where(p => p.DecoyContamTarget == "T" && p.QValue <= 0.01).ToList();
            var peptide = allPeptides.First(p => !p.FullSequence.Contains("["));

            Random rand = new Random();
            var protein = allProteins.FirstOrDefault(p => p.Accession == peptide.ProteinAccession);
            var mutatedSequence = RandomlyMutateAminoAcid(peptide.FullSequence);
            var pos = peptide.StartAndEndResiduesInProtein.Trim('[', ']').Split(new[] { " to " }, StringSplitOptions.None);
            var newProteinSequence = protein.BaseSequence.Substring(0, int.Parse(pos[0]) - 1) + mutatedSequence + protein.BaseSequence.Substring(int.Parse(pos[1]));
            var modifiedProtein = new Protein(protein, newProteinSequence);
            Assert.That(newProteinSequence.Length, Is.EqualTo(protein.BaseSequence.Length));
        }
    }
}
