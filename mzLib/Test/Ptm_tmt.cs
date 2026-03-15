using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Readers;
using NUnit.Framework;
using System.Windows.Markup;

namespace Test
{
    public class Ptm_tmt
    {
        [Test]
        public static void TransferIds()
        {
            string fullSequence = "AAGGAGAQVGGSISSGSSASSVTVTR";
            string converted = TmtPair.ConvertSequence(fullSequence);
            string fullSeq2 = "EVIITPNSAWGGEGSLGC[Common Fixed:Carbamidomethyl on C]GIGYGYLHR";
            string converted2 = TmtPair.ConvertSequence(fullSeq2);
            string fullSeq3 = "FGIVTSSAGTGTTEDTEAK[]";
            string converted3 = TmtPair.ConvertSequence(fullSeq3);
        }
    }

    public class TmtPair
    {
        public string FullSequence { get; set; }
        public string BaseSequence { get; set; }
        public int Fraction_LF { get; set; }
        public int Fraction_TMT { get; set; }

        public PsmFromTsv BestPsm_LF { get; set; }
        public PsmFromTsv BestPsm_TMT { get; set; }

        public static string ConvertSequence(string fullSequence)
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
    }
}
