using System.IO.Compression;
using Sequence.Position;

namespace PolyploidQtlSeqCore.QtlAnalysis.Preprocess.IO
{
    /// <summary>
    /// VCFファイルパーサ
    /// </summary>
    internal static class VcfFileParser
    {
        private const string COMMENT = "#";
        private static readonly char[] _delimiter = new[] { '\t' };
        private static readonly char[] _altAlleleDelimiter = new[] { ',' };
        private static readonly char[] _infoDelimiter = new[] { ';' };
        private static readonly char[] _infoKeyValueDelimiter = new[] { '=' };

        private const string _annotationKey = "ANN=";

        private const int _chrIndex = 0;
        private const int _startIndex = 1;
        private const int _refAlleleIndex = 3;
        private const int _altAllelesIndex = 4;
        private const int _infoIndex = 7;
        private const int _formatIndex = 8;
        private const int _basicInfoCount = 9;

        /// <summary>
        /// VCF変異情報を取得する。
        /// </summary>
        /// <param name="filePath">VCFファイルのPath</param>
        /// <returns>Vcf変異情報リスト</returns>
        public static IReadOnlyList<VcfVariant> Parse(string filePath)
        {
            var variantList = new List<VcfVariant>();

            using var readStream = new FileStream(filePath, FileMode.Open, FileAccess.Read);
            using var gzStream = new GZipStream(readStream, CompressionMode.Decompress);
            using var reader = new StreamReader(gzStream);

            while (!reader.EndOfStream)
            {
                var line = reader.ReadLine();
                if (string.IsNullOrEmpty(line)) continue;
                if (line.StartsWith(COMMENT)) continue;

                var variant = ConvertToVcfVariant(line);
                variantList.Add(variant);
            }

            return variantList;
        }

        /// <summary>
        /// VCF行テキストを変異情報に変換する。
        /// </summary>
        /// <param name="line">VCF行テキスト</param>
        /// <returns>VCF変異</returns>
        private static VcfVariant ConvertToVcfVariant(string line)
        {
            var items = line.Split(_delimiter);

            var chr = items[_chrIndex];
            var start = int.Parse(items[_startIndex]);
            var refAllele = items[_refAlleleIndex];
            var altAlleles = items[_altAllelesIndex].Split(_altAlleleDelimiter).ToArray();
            var format = items[_formatIndex];

            var allAlleles = new[] { refAllele }.Concat(altAlleles).ToArray();
            var pos = GetGenomePosition(chr, start, refAllele);
            var type = GetType(allAlleles);

            var sampleValues = items.Skip(_basicInfoCount).ToArray();
            var samples = new RawVcfSamples(altAlleles, format, sampleValues);
            var snpEffAnns = GetSnpEffAnnotations(items[_infoIndex]);

            return new VcfVariant(
                pos,
                refAllele,
                type,
                samples.IsMultiAllele,
                samples.GetVcfParent1(allAlleles),
                samples.GetVcfParent2(allAlleles),
                samples.GetVcfBulk1(allAlleles),
                samples.GetVcfBulk2(allAlleles),
                snpEffAnns);
        }

        /// <summary>
        /// ゲノム位置を取得する。
        /// </summary>
        /// <param name="chrName">染色多名</param>
        /// <param name="start">Start</param>
        /// <param name="refAllele">refアレル</param>
        /// <returns>ゲノム位置</returns>
        private static GenomePosition GetGenomePosition(string chrName, int start, string refAllele)
        {
            var end = start + refAllele.Length - 1;
            return new GenomePosition(chrName, start, end);
        }

        /// <summary>
        /// 変異の種類を取得する。
        /// snpかそれ以外かに分ける。
        /// </summary>
        /// <param name="allAlleles">全アレル</param>
        /// <returns>種類</returns>
        private static VariantType GetType(string[] allAlleles)
        {
            return allAlleles.All(x => x.Length == 1)
                ? VariantType.snp
                : VariantType.other;
        }

        /// <summary>
        /// SnpEffアノテーションコレクションを作成する。
        /// </summary>
        /// <param name="info">INFO情報</param>
        /// <returns>SnpEffアノテーションコレクション</returns>
        private static SnpEffAnnotations GetSnpEffAnnotations(string info)
        {
            if (!info.Contains(_annotationKey)) return new SnpEffAnnotations();

            var keyValuePairs = info.Split(_infoDelimiter);
            var annLine = keyValuePairs.First(x => x.StartsWith(_annotationKey));
            var items = annLine.Split(_infoKeyValueDelimiter);

            return SnpEffAnnotationParser.Parse(items[1]);
        }
    }
}
