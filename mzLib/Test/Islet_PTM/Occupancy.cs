using NUnit.Framework;
using Readers;
using Readers.Generated;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
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
            var allPsms = allPSMsTMT_file.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T" && !p.FullSequence.Contains("|"));
            var allPeptides_path = @"E:\Islets\Brian_data\PTM\MS3_all_LFgptmdFilterPrunedDb\Task1-SearchTask\AllPeptides.psmtsv";
            var allPeptidesTMT_file = new PsmFromTsvFile(allPeptides_path);
            var allPeptides = allPeptidesTMT_file.Results.Where(p => p.QValue <= 0.01 && p.DecoyContamTarget == "T" && !p.FullSequence.Contains("|"));
            var allPeptidesWithMod = allPeptides.Where(p => SpectrumMatchFromTsv.ParseModifications(p.FullSequence).Values.Any(v => !notInteresting.Any(key => v.Contains(key))));
            
            var allPtms = new List<PtmRecord>();
            foreach (var peptide in allPeptidesWithMod)
            {
                var allMods = SpectrumMatchFromTsv.ParseModifications(peptide.FullSequence).Where(kvp => !notInteresting.Any(key => kvp.Value.Contains(key)));
                int startAA = int.Parse(peptide.StartAndEndResiduesInParentSequence.Split()[0].Split("[")[1]);
                foreach (var mod in allMods)
                {
                    int modSite = mod.Key + startAA;
                    if (mod.Key != 0) modSite = modSite - 1;
                    allPtms.Add(new PtmRecord { Accession = peptide.ProteinAccession, ModSite = modSite, Mod = mod.Value });
                }
            }
            var allPtmGroups = allPtms.GroupBy(ptm => new {a = ptm.Accession, m = ptm.ModSite}).ToList();
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
