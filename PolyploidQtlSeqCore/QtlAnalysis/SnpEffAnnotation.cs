using PolyploidQtlSeqCore.QtlAnalysis.IO;

namespace PolyploidQtlSeqCore.QtlAnalysis
{
    /// <summary>
    /// SnpEff アノテーション情報
    /// </summary>
    internal class SnpEffAnnotation
    {
        private const string DELIMITER = "=";

        /// <summary>
        /// SnpEffアノテーション情報を作成する。
        /// </summary>
        /// <param name="id">Gene ID</param>
        /// <param name="impact">Impact</param>
        /// <param name="ann">アノテーション</param>
        /// <param name="hgvsc">HGVS.C</param>
        /// <param name="hgvsp">HGVS.P</param>
        public SnpEffAnnotation(string id, string impact, string ann, string hgvsc, string hgvsp)
        {
            if (string.IsNullOrWhiteSpace(id)) throw new ArgumentException(null, nameof(id));
            if (string.IsNullOrWhiteSpace(impact)) throw new ArgumentException(null, nameof(impact));
            if (string.IsNullOrWhiteSpace(ann)) throw new ArgumentException(null, nameof(ann));

            GeneId = id;
            Impact = (Impact)Enum.Parse(typeof(Impact), impact);
            if (Impact == Impact.None) throw new ArgumentException("Impact is None.");

            Annotation = ann;
            HgvsC = hgvsc;
            HgvsP = hgvsp;
        }

        /// <summary>
        /// Gene Idを取得する。
        /// </summary>
        public string GeneId { get; }

        /// <summary>
        /// Impactを取得する。
        /// </summary>
        public Impact Impact { get; }

        /// <summary>
        /// Annotationを取得する。
        /// </summary>
        public string Annotation { get; }

        /// <summary>
        /// HGVS.Cを取得する。
        /// </summary>
        public string HgvsC { get; }

        /// <summary>
        /// HGVS.Pを取得する。
        /// </summary>
        public string HgvsP { get; }

        /// <summary>
        /// アノテーションを表示するかどうかを判断する。
        /// </summary>
        /// <param name="displayImpacts">表示するImapactフラグ</param>
        /// <returns>表示するならtrue</returns>
        public bool Display(DisplayAnnotationImpacts displayImpacts)
        {
            return displayImpacts.DisplayFlag.HasFlag(Impact);
        }

        /// <summary>
        /// 出力用Impact情報に変換する。
        /// </summary>
        /// <returns>Impact情報</returns>
        public string ToImpactInfo()
        {
            return $"{GeneId}{DELIMITER}{Impact}";
        }

        /// <summary>
        /// 出力用Annotation情報に変換する。
        /// </summary>
        /// <returns>Annotation情報</returns>
        public string ToAnnotationInfo()
        {
            return $"{GeneId}{DELIMITER}{Annotation}";
        }

        /// <summary>
        /// 出力用HGVS.C情報に変換する。
        /// </summary>
        /// <returns>HGVS.C情報</returns>
        public string ToHgvsCInfo()
        {
            return $"{GeneId}{DELIMITER}{HgvsC}";
        }

        /// <summary>
        /// 出力用HGVS.P情報に変換する。
        /// </summary>
        /// <returns>HGVS.P情報</returns>
        public string ToHgvsPInfo()
        {
            return string.IsNullOrWhiteSpace(HgvsP)
                ? ""
                : $"{GeneId}{DELIMITER}{HgvsP}";
        }

    }
}
