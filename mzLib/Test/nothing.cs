using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Readers;
using NUnit.Framework;

namespace Test
{
    public class nothing
    {
        [Test]
        public static void TestReading()
        {

            var osmPath = @"E:\Temp\AllOSMs.osmtsv";
            var allOsms = SpectrumMatchTsvReader.ReadOsmTsv(osmPath, out List<string> warnings).ToList();
        }
    }
}
