using CsvHelper.Configuration;
using CsvHelper.Configuration.Attributes;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Readers
{
    public class TmtProteinResult
    {
        public static CsvConfiguration CsvConfiguration = new CsvConfiguration(System.Globalization.CultureInfo.InvariantCulture)
        {
            Delimiter = "\t",
            HasHeaderRecord = true,
            IgnoreBlankLines = true,
            TrimOptions = CsvHelper.Configuration.TrimOptions.Trim,
        };

        [Name("Protein Accession")]
        public string Accession { get; set; }

        [Name("Number of Proteins in Group")]
        public string NumProteinsInGroup { get; set; }

        [Name("Unique Peptides")]
        public string UniquePeptides { get; set; }

        [Name("Shared Peptides")]
        public string SharedPeptides { get; set; }

        [Name("Protein QValue")]
        public double ProteinQValue { get; set; }


    }
}
