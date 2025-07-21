using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NUnit.Framework;
using Readers;

namespace Test
{
    public class TestPsmValidations
    {
        [Test]
        public static void Test()
        {
            var psmTsvPath = @"E:\Aneuploidy\DDA\071525\1614_E1-8_calied-generalGPTMD+1NAsub_noTrunc\Task2-SearchTask\Individual File Results\07-15-25_1614-R1-Q_E1+5-calib_Peptides.psmtsv";
            var allPsmTsv = SpectrumMatchTsvReader.ReadTsv(psmTsvPath, out List<string> warnings).Where(p => p.DecoyContamTarget == "T" && p.QValue <= 0.01).ToList();
        }
    }
}
