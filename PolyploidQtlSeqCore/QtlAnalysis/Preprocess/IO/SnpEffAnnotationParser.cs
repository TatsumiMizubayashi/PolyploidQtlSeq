using System.Web;

namespace PolyploidQtlSeqCore.QtlAnalysis.Preprocess.IO
{
    /// <summary>
    /// SnpEffアノテーション パーサー
    /// </summary>
    internal static class SnpEffAnnotationParser
    {
        /// <summary>
        /// アノテーションの種類(Allele, Annotationなど）を区切るためのデリミタ
        /// </summary>
        private static readonly char[] _annTypeDelimiter = new[] { '|' };

        /// <summary>
        /// 複数のアノテーション情報を区切るためのデリミタ
        /// </summary>
        private static readonly char[] _lineDelimiter = new[] { ',' };

        private const int _annotationIndex = 1;
        private const int _impactIndex = 2;
        private const int _geneIdIndex = 4;
        private const int _hgvsCIndex = 9;
        private const int _hgvsPIndex = 10;

        /// <summary>
        /// SnpEffアノテーションコレクションに変換する。
        /// </summary>
        /// <param name="annValue">ANN値</param>
        /// <returns>SnpEffアノテーションコレクション</returns>
        public static SnpEffAnnotations Parse(string annValue)
        {
            if (string.IsNullOrWhiteSpace(annValue)) throw new ArgumentException(null, nameof(annValue));

            // URLエンコードされている場合があるのでデコードしておく
            var annLines = annValue.Split(_lineDelimiter)
                .Select(x => HttpUtility.UrlDecode(x))
                .ToArray();

            var snpEffAnnotations = new SnpEffAnnotation[annLines.Length];
            for (var i = 0; i < annLines.Length; i++)
            {
                var items = annLines[i].Split(_annTypeDelimiter);
                var id = items[_geneIdIndex];
                var impact = items[_impactIndex];
                var ann = items[_annotationIndex];
                var hgvsc = items[_hgvsCIndex];
                var hgvsp = items[_hgvsPIndex];

                snpEffAnnotations[i] = new SnpEffAnnotation(id, impact, ann, hgvsc, hgvsp);
            }

            return new SnpEffAnnotations(snpEffAnnotations);
        }
    }
}
